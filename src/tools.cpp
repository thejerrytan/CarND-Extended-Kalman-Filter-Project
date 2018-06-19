#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

const float EPS = 0.0001;

Tools::Tools() {
  rmse = VectorXd(4);
  rmse << 0,0,0,0;
}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  
  // Check the validity of inputs
  assert(estimations.size() > 0);
  assert(ground_truth.size() == estimations.size());

  // Smarter way, we keep a running RMSE and add the contribution from the last reading
  int lastIdx = estimations.size() - 1;
  VectorXd sumOfSquares = (rmse.array().square() * (estimations.size() - 1)) + ((estimations[lastIdx] - ground_truth[lastIdx]).array() * (estimations[lastIdx] - ground_truth[lastIdx]).array());
  rmse = sumOfSquares / estimations.size();
  rmse = rmse.array().sqrt();

  // Naive way of calculating RSME, a lot of repeated calculations
  // rmse << 0,0,0,0;
  // for (unsigned int i = 0; i < estimations.size(); ++i) {
  //   VectorXd rse = (estimations[i] - ground_truth[i]).array() * (estimations[i] - ground_truth[i]).array();
  //   rmse += rse;
  // }
  // rmse /= estimations.size();
  // rmse = rmse.array().sqrt();

  cout << "RMSE = " << rmse << endl;

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);

  float norm = sqrt(px*px + py*py);
  float norm2 = px*px + py*py;
  float norm3 = pow(norm, 3);
  
  //check division by zero
  if (norm == 0 || norm2 == 0 || norm3 == 0) {
    std::cout << "Error: divisionByZero" << std::endl;
    const double c2 = std::max(EPS, px*px + py*py); // equivalent of norm2
    const double c1 = sqrt(c2); // equivalent of norm
    const double c3 = pow(c1, 3); // equivalent of norm3
    Hj(0,0) = px / c1;
    Hj(0,1) = py / c1;
    Hj(0,2) = 0;
    Hj(0,3) = 0;
    Hj(1,0) = -py / c2;
    Hj(1,1) = px / c2;
    Hj(1,2) = 0;
    Hj(1,3) = 0;
    Hj(2,0) = py * (vx * py - vy * px) / c3;
    Hj(2,1) = px * (vy * px - vx * py) / c3;
    Hj(2,2) = px / c1;
    Hj(2,3) = py / c1;
  } else {
    //compute the Jacobian matrix
    Hj(0,0) = px / norm;
    Hj(0,1) = py / norm;
    Hj(0,2) = 0;
    Hj(0,3) = 0;
    Hj(1,0) = -py / norm2;
    Hj(1,1) = px / norm2;
    Hj(1,2) = 0;
    Hj(1,3) = 0;
    Hj(2,0) = py * (vx * py - vy * px) / norm3;
    Hj(2,1) = px * (vy * px - vx * py) / norm3;
    Hj(2,2) = px / norm;
    Hj(2,3) = py / norm;
  }

  return Hj;
}
