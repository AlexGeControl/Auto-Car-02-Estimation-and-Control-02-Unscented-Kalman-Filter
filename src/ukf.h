#ifndef UKF_H
#define UKF_H

#include <vector>
#include <string>
#include <fstream>
#include "Eigen/Dense"
#include "measurement_package.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;
  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* state dimension
  static const int n_x_;
  ///* augmented state dimension
  static const int n_aug_;
  ///* sigma point spreading parameter
  double lambda_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;
  ///* state covariance matrix
  MatrixXd P_xx_;

  ///* process noise -- standard deviation velocity magnitude in m/s
  double std_v_;
  ///* process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;
  ///* process noise -- standard deviation yaw angle in rad
  double std_yaw_;
  ///* process noise -- standard deviation yaw velocity in rad/s
  double std_yawd_;
  ///* process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* mean of sigma points:
  VectorXd x_aug_;
  ///* covariance of sigma points:
  MatrixXd P_aug_;
  ///* weights of sigma points
  VectorXd weights_;
  ///* sigma points matrix:
  MatrixXd Xsig_;
  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;
  ///* deviation between predicted sigma points and their mean
  MatrixXd Xsig_pred_dev_;

  ///* measurement dimension, laser:
  static const int n_z_laser_;
  ///* measurement dimension, radar:
  static const int n_z_radar_;

  ///* predicted observation sigma points matrix:
  MatrixXd Zsig_laser_pred_;
  MatrixXd Zsig_radar_pred_;
  // predicted mean observation:
  VectorXd z_laser_;
  VectorXd z_radar_;
  ///* deviation between predicted observations and their mean:
  MatrixXd Zsig_laser_pred_dev_;
  MatrixXd Zsig_radar_pred_dev_;
  // predicted observation covariance matrix:
  MatrixXd P_zz_laser_;
  MatrixXd P_zz_radar_;
  // predicted state-observation correlation matrix:
  MatrixXd P_xz_laser_;
  MatrixXd P_xz_radar_;

  ///* Kalman gains:
  MatrixXd K_laser_;
  MatrixXd K_radar_;

  ///* laser measurement noise standard deviation position1 in m
  double std_laspx_;
  ///* laser measurement noise standard deviation position2 in m
  double std_laspy_;
  ///* laser measurement noise covariance matrix:
  MatrixXd R_laser_;

  ///* radar measurement noise standard deviation radius in m
  double std_radr_;
  ///* radar measurement noise standard deviation angle in rad
  double std_radphi_;
  ///* radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;
  ///* radar measurement noise covariance matrix:
  MatrixXd R_radar_;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Generate sigma points
   */
  void GenerateSigmaPoints(void);

  /**
   * Make state prediction using discretized Constant Turing Rate and Velocity(CTRV) model
   * @param {const VectorXd &} x current state
   * @param {const double} delta_t time elapsed between current and last state
   */
  VectorXd DoPredict(const VectorXd &x, const double delta_t);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Predict(double delta_t);

  /**
   * Make observation prediction using lidar measurement model
   * @param {const VectorXd &} predicted k + 1 state using k state
   */
  VectorXd DoUpdateLidar(const VectorXd &x);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param z The measurement at k+1
   */
  void UpdateLidar(const VectorXd &z);

  /**
   * Make observation prediction using radar measurement model
   * @param {const VectorXd &} predicted k + 1 state using k state
   */
  VectorXd DoUpdateRadar(const VectorXd &x);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param z The measurement at k+1
   */
  void UpdateRadar(const VectorXd &z);
};

#endif /* UKF_H */
