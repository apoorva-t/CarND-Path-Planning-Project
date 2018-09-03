/*
 * jmt.cpp
 *
 *  Created on: Sep 3, 2018
 *      Author: apoorva
 */

#include "jmt.h"
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "Eigen-3.3/Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

vector<double> getJMT(const vector<double>& start, const vector<double>& end, double T)
{
	MatrixXd A = MatrixXd(3, 3);
	A << T*T*T, T*T*T*T, T*T*T*T*T,
				    3*T*T, 4*T*T*T,5*T*T*T*T,
				    6*T, 12*T*T, 20*T*T*T;

	MatrixXd B = MatrixXd(3,1);
	B << end[0]-(start[0]+start[1]*T+.5*start[2]*T*T),
				    end[1]-(start[1]+start[2]*T),
				    end[2]-start[2];

	MatrixXd Ai = A.inverse();

	MatrixXd C = Ai*B;

	vector <double> result = {start[0], start[1], .5*start[2]};
	for(int i = 0; i < C.size(); i++)
	{
		result.push_back(C.data()[i]);
	}

	return result;
}
