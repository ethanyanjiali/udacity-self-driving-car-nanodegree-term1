#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::endl;
using std::cout;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  if(estimations.size() != ground_truth.size()){
    cout << "ERROR - CalculateRMSE() - estimations and ground truth should be in same size." << endl;
    return rmse;
  }
  
  // square error
  for (unsigned int i = 0; i < estimations.size(); i++) {
    VectorXd d = estimations[i] - ground_truth[i];
    d = d.array() * d.array();
    rmse += d;
  }
  // mean
  rmse = rmse / estimations.size();
  // root
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3,4);
  if (x_state.size() != 4) {
    cout << "ERROR - CalculateJacobian() - The state vector should have size 4." << endl;
    return Hj;
  }
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  double c1 = px * px + py * py;
  double c2 = sqrt(c1);
  double c3 = c1 * c2;
  
  if (fabs(c1) < 0.0001) {
    cout << "CalculateJacobian() - Error - Division by Zero" << endl;
    return Hj;
  }
  Hj << (px/c2), (py/c2), 0, 0,
        -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  return Hj;
}
