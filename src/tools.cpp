#include <iostream>
#include "tools.h"
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
    * Calculate the RMSE here.
  */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  if(estimations.size() == 0 || estimations.size() != ground_truth.size()) {
      cout << "Invalid estimation data, size is 0 or not equal." << endl;
      return rmse;
  }

  for(int i = 0; i < estimations.size(); i++) {
      VectorXd err = estimations[i] - ground_truth[i];
      err = err.array() * err.array();
      rmse += err;
  }

  rmse = rmse.array() / estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
    * Calculate a Jacobian here.
  */
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  float power_xy = px*px + py*py;

  MatrixXd Hj(3,4);

  if(fabs(power_xy) < 0.0001) {
      cout << "cannot divided by zero." << endl;
      return Hj;
  }

  float s1 = sqrt(power_xy);
  float s2 = power_xy * s1;

  Hj(0, 0) = px/s1;
  Hj(0, 1) = py/s1;
  Hj(0, 2) = 0;
  Hj(0, 3) = 0;
  Hj(1, 0) = -py/power_xy;
  Hj(1, 1) = px/power_xy;
  Hj(1, 2) = 0;
  Hj(1, 3) = 0;
  Hj(2, 0) = py*(vx*py-vy*px)/s2;
  Hj(2, 1) = px*(vy*px-vx*py)/s2;
  Hj(2, 2) = px/s1;
  Hj(2, 3) = py/s1;

  return Hj;

}
