#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"
#include "jmt.h"

using namespace std;
#define NUM_TOTAL_STATES 3

const double SPEED_DIFF = 0.06;
const double MAX_SPEED = 21.6;
const double STOP_COST = 0.8;
const double BUFFER_SPEED = 2.0;
const double LC_COST = 0.3;
const double LC_BUFFER = 25.0;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

void addPointsToSplineFit(vector<double>& ptsx, vector<double>& ptsy, double& car_x, double& car_y, double end_s, int lane, double& ref_yaw, const vector<double>& previous_path_x, const vector<double>& previous_path_y,
						  const vector<double>& map_waypoints_s, const vector<double>& map_waypoints_x, const vector<double>& map_waypoints_y)
{
	int prev_path_size = previous_path_x.size();
	if (prev_path_size < 2)
	{
		ptsx.push_back(car_x - cos(ref_yaw));
		ptsx.push_back(car_x);

		ptsy.push_back(car_y - sin(ref_yaw));
		ptsy.push_back(car_y);

		//std::cout << "X points 1 and 2: " << ptsx[0] << " " << ptsx[1] << std::endl;
		//std::cout << "Y points 1 and 2: " << ptsy[0] << " " << ptsy[1] << std::endl;
	}
	else
	{
		car_x = previous_path_x[prev_path_size - 1];
		car_y = previous_path_y[prev_path_size - 1];
		// compute car_s from car_x and car_y
		//vector<double> frenet = getFrenet(car_x, car_y, ref_yaw, map_waypoints_x, map_waypoints_y);
		//car_s = frenet[0];
		ptsx.push_back(previous_path_x[prev_path_size - 2]);
		ptsx.push_back(car_x);

		ptsy.push_back(previous_path_y[prev_path_size - 2]);
		ptsy.push_back(car_y);
		//std::cout << "X points 1 and 2: " << ptsx[0] << " " << ptsx[1] << std::endl;
		//std::cout << "Y points 1 and 2: " << ptsy[0] << " " << ptsy[1] << std::endl;

		ref_yaw = atan2((car_y - previous_path_y[prev_path_size - 2]), (car_x - previous_path_x[prev_path_size - 2]));
	}
	/* Add anchor points in intended lane from end of previous path */
	//vector<double> next_wp00 = getXY(end_s + 20, (2+4*lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_wp0 = getXY(end_s + 30, (2+4*lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_wp1 = getXY(end_s + 50, (2+4*lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_wp2 = getXY(end_s + 70, (2+4*lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
	vector<double> next_wp3 = getXY(end_s + 90, (2+4*lane),map_waypoints_s, map_waypoints_x, map_waypoints_y);
	//std::cout << "Frenet value: " << end_s << std::endl;

	//ptsx.push_back(next_wp00[0]);
	ptsx.push_back(next_wp0[0]);
	ptsx.push_back(next_wp1[0]);
	ptsx.push_back(next_wp2[0]);
	ptsx.push_back(next_wp3[0]);

	//std::cout << "X points 3, 4 and 5: " << ptsx[2] << " " << ptsx[3] << " " << ptsx[4] << std::endl;
	//ptsy.push_back(next_wp00[1]);
	ptsy.push_back(next_wp0[1]);
	ptsy.push_back(next_wp1[1]);
	ptsy.push_back(next_wp2[1]);
	ptsy.push_back(next_wp3[1]);
	//std::cout << "Y points 3, 4 and 5: " << ptsy[2] << " " << ptsy[3] << " " << ptsy[4] << std::endl;
}

/* Find the lane in which car is based on its 'd' value*/
int getCarLane(double d)
{
	if (d >= 0 && d < 4.0)
		return 0;
	else if (d >= 4.0 && d < 8.0)
		return 1;
	else if (d >= 8.0 && d<=12.0)
		return 2;
	else
		return -1;
}

struct State {
	int lane;
	double velocity;
	bool lc_allowed;
};

inline double getSpeedLimitCost(double vel) {
	if (vel > MAX_SPEED) return 1.0;
	double target = MAX_SPEED - BUFFER_SPEED;
	return (STOP_COST * (target - vel) / target);
}

State getLowestCostState(int current_lane, bool car_in_front, bool car_on_left, bool car_on_right, double ref_vel,
		double front_car_s, double left_car_s, double right_car_s, double car_s, double ref_left_speed, double ref_right_speed,
		double left_car_behind_s, double right_car_behind_s, bool lane_change_allowed, int intended_lane)
{
  	std::cout << "Lane Change Allowed: " << lane_change_allowed << std::endl;
	State nextState;
	if (!car_in_front) // no need to change lanes
	{
		if (ref_vel < MAX_SPEED) ref_vel += SPEED_DIFF;
		nextState.lane = current_lane;
		nextState.velocity = ref_vel;
		nextState.lc_allowed = true;
		std::cout << "Next intended lane (keep): " << nextState.lane << "current lane: " << current_lane << std::endl;
		return nextState;
	}
	/*
	 * Cost fn1: Speed limit- since car cannot exceed speed limit (v > MAX_SPEED) => cost = 1
	 *           if (v < MAX_SPEED - BUFFER) cost = STOP_COST * (MAX_SPEED - BUFFER - v) / (MAX_SPEED - BUFFER)
	 * Cost fn2: Penalize a lane change anyway - cost = 0.8
	 * Cost fn3: Penalize getting too close to car in front - cost = SDIFF_COST * (check_car_s - car_s)
	 * Cost fn4: Penalize getting in way of car behind during lane change
	 */
	double cost[NUM_TOTAL_STATES] = {1,1,1}; //costs for KL, LCL, LCR
	//for keeping current lane
	{
		double target_speed = ref_vel - SPEED_DIFF;
		double cost1 = getSpeedLimitCost(target_speed);
		double cost2 = 0.0; //lane change cost
		double cost3 = 1.0 - exp(-1.0 / (double)(front_car_s - car_s));
		cost[0] = 0.25 * cost1 + 0.25 * cost2 + 0.25 * cost3;
		//std::cout << "Cost of KL is: " << cost[0] << std::endl;
	}
	if (lane_change_allowed && !car_on_left && (current_lane != 0)) // if car_on_left or already in left lane then can't perform LCL
	{
		double target_speed = ref_left_speed - SPEED_DIFF;
		double cost1 = getSpeedLimitCost(target_speed);
		double cost2 = LC_COST;
		double cost3 = 1.0 - exp(-1.0 / (double)(left_car_s - car_s));
		double cost4 = 1.0 - exp(-1.0 / (double)(left_car_behind_s - car_s));
		cost[1] = 0.25 * cost1 + 0.15 * cost2 + 0.3 * cost3 + 0.3*cost4;
		//std::cout << "Cost of LCL is: " << cost[0] << std::endl;
	}
	if (lane_change_allowed && !car_on_right && (current_lane != 2))
	{
		double target_speed = ref_right_speed - SPEED_DIFF;
		double cost1 = getSpeedLimitCost(target_speed);
		double cost2 = LC_COST;
		double cost3 = 1.0 - exp(-1.0 / (double)(right_car_s - car_s));
		double cost4 = 1.0 - exp(-1.0 / (double)(right_car_behind_s - car_s));
		cost[2] = 0.25 * cost1 + 0.15 * cost2 + 0.3 * cost3 + 0.3*cost4;
		//std::cout << "Cost of LCR is: " << cost[0] << std::endl;
	}
	int minCostIndex = 0;
	double minCost = 1.0;
	for (int i = 0; i < 3; i++)
	{
		if (cost[i] < minCost) minCostIndex = i;
	}
	switch(minCostIndex) {
	case 0:
		nextState.lane = intended_lane;
		nextState.velocity = ref_vel - SPEED_DIFF;
		nextState.lc_allowed = true;
		break;
	case 1:
		nextState.lane = current_lane - 1;
		nextState.velocity = ref_vel;
		nextState.lc_allowed = false;
		break;
	case 2:
		nextState.lane = current_lane + 1;
		nextState.velocity = ref_vel;
		nextState.lc_allowed = false;
		break;
	}
	std::cout << "Next intended lane: " << nextState.lane << "current lane: " << current_lane << std::endl;
	return nextState;
}

int computeNextLaneAndSpeed(const vector<vector<double> >& sensor_fusion, int current_lane, double& ref_vel, int prev_path_size, double car_s, bool& lane_change_allowed, int intended_lane)
{
	bool car_in_front = false, car_on_left = false, car_on_right = false;
	double front_car_s = 0.0, left_car_s = 0.0, right_car_s = 0.0;
	double left_car_behind_s = 0.0, right_car_behind_s = 0.0;
	double ref_left_speed = ref_vel, ref_right_speed = ref_vel;
	for (int i = 0; i < sensor_fusion.size(); i++)
	{
		// for each car, check:
		int lane = getCarLane(sensor_fusion[i][6]);
		if (lane < 0) continue;

		// compute car speed and position
		double vx = sensor_fusion[i][3];
		double vy = sensor_fusion[i][4];
		double check_speed = sqrt(vx*vx + vy*vy);

		double check_car_s = sensor_fusion[i][5];
		//estimated position of car at the end of previous path
		check_car_s += ((double)prev_path_size*0.02*check_speed);

		if (lane == current_lane)
		{
			// if car is within LC_BUFFER of our car
			if (check_car_s > car_s && (check_car_s - car_s < LC_BUFFER))
			{
				if (front_car_s == 0.0 || check_car_s < front_car_s) front_car_s = check_car_s;
				car_in_front = true;
			}
		}
		else if (current_lane - lane == 1) // car in left lane
		{
			if (check_car_s < car_s+LC_BUFFER && check_car_s > car_s - LC_BUFFER) {
				car_on_left = true;
			}
			if (check_car_s > car_s && (left_car_s == 0.0 || check_car_s < left_car_s)) {
				left_car_s = check_car_s;
				ref_left_speed = check_speed;
			}
			if (check_car_s < car_s && (left_car_behind_s == 0.0 || check_car_s > left_car_behind_s)) left_car_behind_s = check_car_s;
		}
		else if (current_lane - lane == -1)
		{
			if (check_car_s < car_s+LC_BUFFER && check_car_s > car_s - LC_BUFFER) {
				car_on_right = true;
			}
			if (check_car_s > car_s && (right_car_s == 0.0 || check_car_s < right_car_s)) {
				right_car_s = check_car_s;
				ref_right_speed = check_speed;
			}
			if (check_car_s < car_s && (right_car_behind_s == 0.0 || check_car_s > right_car_behind_s)) right_car_behind_s = check_car_s;
		}
	}

	State nextState = getLowestCostState(current_lane, car_in_front, car_on_left, car_on_right, ref_vel, front_car_s, left_car_s, right_car_s, car_s,
			ref_left_speed, ref_right_speed, left_car_behind_s, right_car_behind_s, lane_change_allowed, intended_lane);
	ref_vel = nextState.velocity;
	lane_change_allowed = nextState.lc_allowed;
	return nextState.lane;
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  double ref_velocity = 0.0; // in m/s
  bool lane_change_allowed = true;
  int intended_lane = -1;

  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy, &ref_velocity, &lane_change_allowed, &intended_lane](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

          	int prev_path_size = previous_path_x.size();
          	// Since we will continue on unprocessed previous_path co-ordinates
          	if (prev_path_size > 0) {
          		car_s = end_path_s;
          		car_d = end_path_d;
          	}
          	int current_lane = getCarLane(car_d); // 0, 1 or 2;
          	if (current_lane < 0) std::cout << "ERROR!! Invalid lane!!" << std::endl;
          	if (intended_lane < 0) intended_lane = current_lane;

          	int next_lane = computeNextLaneAndSpeed(sensor_fusion, current_lane, ref_velocity, prev_path_size, car_s, lane_change_allowed, intended_lane);

          	vector<double> next_x_vals;
          	vector<double> next_y_vals;

          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
          	// first get the car to drive in its own lane smoothly

          	/* Add leftover previous waypoints */
          	for (int i= 0; i < prev_path_size; i++)
          	{
          		next_x_vals.push_back(previous_path_x[i]);
          		next_y_vals.push_back(previous_path_y[i]);
          	}

          	vector<double> ptsx;
          	vector<double> ptsy;
          	/* Add 2 points from previous path for spline fit */
          	double ref_x = car_x, ref_y = car_y, ref_yaw = deg2rad(car_yaw);
          	double ref_s = car_s;
          	addPointsToSplineFit(ptsx, ptsy, ref_x, ref_y, ref_s, next_lane, ref_yaw, previous_path_x, previous_path_y, map_waypoints_s, map_waypoints_x, map_waypoints_y);

          	/* Shift waypoints to car's co-ordinate system */
          	for (int i = 0; i < ptsx.size(); i++)
          	{
          		double mod_x = ptsx[i] - ref_x, mod_y = ptsy[i] - ref_y;
          		ptsx[i] = mod_x * cos(0-ref_yaw) - mod_y * sin(0-ref_yaw);
          		ptsy[i] = mod_x * sin(0-ref_yaw) + mod_y * cos(0-ref_yaw);

          		//std::cout << "Spline X: " << ptsx[i] << std::endl;
          		//std::cout << "Spline Y: " << ptsy[i] << std::endl;
          	}

          	/* Generate smooth path use a spline curve fit */
          	tk::spline s;
          	s.set_points(ptsx, ptsy);


        	// trajectory generated using JMT
          	/*
        	vector<double> startS = {car_s, car_speed, 0};
        	vector<double> endS = {car_s+30, ref_velocity, 0};
        	vector<double> startD = {(double)2+4*intended_lane, 0.0, 0.0};
        	vector<double> endD = {(double)2+4*next_lane, 0.0, 0.0};
        	vector<double> s_coeffs = getJMT(startS, endS, 1.5);
        	vector<double> d_coeffs = getJMT(startD, endD, 1.5);
        	*/

          	/* Calculate how to break up spline points to travel at desired velocity */
          	auto targetX = 30.0; //project spline out to horizon at 30 m
          	auto targetY = s(targetX);
          	auto targetDist = sqrt((targetX * targetX) + (targetY * targetY));
          	double x_add_on = 0;
          	for (int i = 0; i < 50 - prev_path_size; i++)
          	{
          		// Get points on the trajectory spaced such that car travels at ref_velocity
          		double x_point = x_add_on + (targetX * (0.02 * ref_velocity)) / targetDist;
          		double y_point = s(x_point);

          		x_add_on = x_point;

          		// Shift points back to global co-ordinate system
          		double x_wp = x_point*cos(ref_yaw) - y_point*sin(ref_yaw) + ref_x;
          		double y_wp = x_point*sin(ref_yaw) + y_point*cos(ref_yaw) + ref_y;

          		next_x_vals.push_back(x_wp);
          		next_y_vals.push_back(y_wp);

          	}
          	intended_lane = next_lane;

          	msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
