#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using std::vector;

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 *
 */
const int UKF::n_x_ = 5;
const int UKF::n_aug_ = 7;

const int UKF::n_z_laser_ = 2;
const int UKF::n_z_radar_ = 3;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initially set to false, set to true in first call of ProcessMeasurement
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;
  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  /**
   * Process
   */
  // sigma point spreading parameter
  lambda_ = 3.0 - n_aug_;
  // initial state vector
  x_ = VectorXd::Zero(n_x_);
  // initial covariance matrix
  P_xx_ = MatrixXd::Zero(n_x_, n_x_);

  // process noise -- standard deviation velocity magnitude in m/s
  std_v_ = 5;
  // process noise -- standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;
  // process noise -- standard deviation yaw angle in rad
  std_yaw_ = M_PI / 4;
  // process noise -- standard deviation yaw velocity in rad/s
  std_yawd_ = M_PI / 12;
  // process noise -- standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI / 12;

  // weights of sigma points
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  weights_.tail(2 * n_aug_) = VectorXd::Constant(2 * n_aug_, 0.5 / (lambda_ + n_aug_));
  // sigma points mean:
  x_aug_ = VectorXd::Zero(n_aug_);
  // sigma points covariance matrix:
  P_aug_ = MatrixXd::Zero(n_aug_, n_aug_);
  P_aug_(n_aug_ - 2, n_aug_ - 2) = pow(std_a_, 2);
  P_aug_(n_aug_ - 1, n_aug_ - 1) = pow(std_yawdd_, 2);
  // sigma points matrix:
  Xsig_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  // deviation between predicted sigma points and their mean:
  Xsig_pred_dev_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  /**
   * Measurement:
   */
  // predicted observation sigma points matrix:
  Zsig_laser_pred_ = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);
  Zsig_radar_pred_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  // predicted mean observation:
  z_laser_ = VectorXd(n_z_laser_);
  z_radar_ = VectorXd(n_z_radar_);
  // deviation between predicted observations and their mean
  Zsig_laser_pred_dev_ = MatrixXd(n_z_laser_, 2 * n_aug_ + 1);
  Zsig_radar_pred_dev_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  // predicted observation covariance matrix:
  P_zz_laser_ = MatrixXd(n_z_laser_, n_z_laser_);
  P_zz_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  // predicted state-observation correlation matrix:
  P_xz_laser_ = MatrixXd(n_x_, n_z_laser_);
  P_xz_radar_ = MatrixXd(n_x_, n_z_radar_);

  // Kalman gains:
  K_laser_ = MatrixXd(n_x_, n_z_laser_);
  K_radar_ = MatrixXd(n_x_, n_z_radar_);

  // measurement noise -- laser, standard deviation position x in m
  std_laspx_ = 0.15;
  // measurement noise -- laser, standard deviation position y in m
  std_laspy_ = 0.15;
  R_laser_ = MatrixXd::Zero(n_z_laser_, n_z_laser_);
  R_laser_(0, 0) = pow(std_laspx_, 2);
  R_laser_(1, 1) = pow(std_laspy_, 2);

  // measurement noise -- radar, standard deviation radius in m
  std_radr_ = 0.3;
  // measurement noise -- radar, standard deviation angle in rad
  std_radphi_ = 0.03;
  // measurement noise -- radar, standard deviation radius change in m/s
  std_radrd_ = 0.3;
  R_radar_ = MatrixXd(n_z_radar_, n_z_radar_);
  R_radar_(0, 0) = pow(std_radr_, 2);
  R_radar_(1, 1) = pow(std_radphi_, 2);
  R_radar_(2, 2) = pow(std_radrd_, 2);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  if (!is_initialized_) {
    /**
     * Initialize the state ekf_.x_ with the first measurement.
     * Create the covariance matrix.
     * Remember: you'll need to convert radar from polar to cartesian coordinates.
     */
    cout << "[UKF]: Initialization..." << endl;

    // first timestamp:
    time_us_ = meas_package.timestamp_;

    // first measurement:
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       * Convert radar from polar to cartesian coordinates and initialize state.
       */
      // parse radar measurements:
      double ro, theta, ro_dot;

      ro = meas_package.raw_measurements_(0);
      theta = meas_package.raw_measurements_(1);
      ro_dot = meas_package.raw_measurements_(2);

      // position:
      x_(0) =  ro * cos(theta);
      x_(1) =  ro * sin(theta);

      // set initial state covariance matrix according to radar measurement noise:
      P_xx_(0, 0) = R_radar_(0, 0);
      P_xx_(1, 1) = R_radar_(0, 0);
      P_xx_(2, 2) = pow(std_v_, 2);
      P_xx_(3, 3) = R_radar_(1, 1);
      P_xx_(4, 4) = pow(std_yawd_, 2);
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
       * Initialize state from lidar measurements:
       */
      double px, py;

      px = meas_package.raw_measurements_(0);
      py = meas_package.raw_measurements_(1);

      // position:
      x_(0) = px;
      x_(1) = py;

      // set initial state covariance matrix according to lidar measurement noise:
      P_xx_(0, 0) = R_laser_(0, 0);
      P_xx_(1, 1) = R_laser_(1, 1);
      P_xx_(2, 2) = pow(std_v_, 2);
      P_xx_(3, 3) = pow(std_yaw_, 2);
      P_xx_(4, 4) = pow(std_yawd_, 2);
    }

    cout << "[UKF]: Initialization done." << endl;

    // done initializing, no need to predict or update
    is_initialized_ = true;

    return;
  }

  /*****************************************************************************
  *  Time elapsed
  ****************************************************************************/
  long long current_time_us = meas_package.timestamp_;
  long long delta_t = current_time_us - time_us_;
  time_us_ = current_time_us;

  /*****************************************************************************
  *  Sampling
  ****************************************************************************/
  GenerateSigmaPoints();

  /*****************************************************************************
  *  Prediction
  ****************************************************************************/
  Predict(delta_t / 1e6);

  /*****************************************************************************
  *  Update
  ****************************************************************************/
  // parse measurements
  const VectorXd &z = meas_package.raw_measurements_;

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(z);
  } else {
    // Laser updates
    UpdateLidar(z);
  }
}

/**
 * Generate sigma points
 */
void UKF::GenerateSigmaPoints(void) {
  // set up distribution params for sigma sampling:
  x_aug_.head(n_x_) = x_;
  P_aug_.topLeftCorner(n_x_, n_x_) = P_xx_;

  // covariance matrix square root through Cholesky decomposition:
  MatrixXd A = P_aug_.llt().matrixL();

  // get sigma points:
  Xsig_.col(0) = x_aug_;
  Xsig_.block<n_aug_, n_aug_>(0,          1) = (+sqrt(lambda_ + n_aug_) * A).colwise() + x_aug_;
  Xsig_.block<n_aug_, n_aug_>(0, n_aug_ + 1) = (-sqrt(lambda_ + n_aug_) * A).colwise() + x_aug_;
}

/**
 * Make state prediction using discretized Constant Turing Rate and Velocity(CTRV) model
 * @param {const VectorXd &} x current state
 * @param {const double} delta_t time elapsed between current and last state
 */
VectorXd UKF::DoPredict(const VectorXd &x, const double delta_t){
    // division-by-zero check:
    const double epsilon = 1e-7;

    // parse augmented state:
    const double p_x = x(0);
    const double p_y = x(1);
    const double v = x(2);
    const double psi = x(3);
    const double psi_dot = x(4);
    const double v_dot = x(5);
    const double psi_ddot = x(6);

    // allocate:
    VectorXd x_pred(5);
    // common terms:
    x_pred << p_x + 0.5 * cos(psi) * v_dot * pow(delta_t, 2),
              p_y + 0.5 * sin(psi) * v_dot * pow(delta_t, 2),
              v + v_dot * delta_t,
              psi + psi_dot * delta_t + 0.50 * psi_ddot * pow(delta_t, 2),
              psi_dot + psi_ddot * delta_t;

    if (abs(psi_dot) < epsilon) {
        // use straight line model:
        x_pred(0) += v * cos(psi) * delta_t;
        x_pred(1) += v * sin(psi) * delta_t;
    } else {
        // use curve line model:
        x_pred(0) += v / psi_dot * (
            sin(psi + psi_dot * delta_t) - sin(psi)
        );
        x_pred(1) += v / psi_dot * (
            -cos(psi + psi_dot * delta_t) + cos(psi)
        );
    }

    return x_pred;
}

/**
* Predicts sigma points, the state, and the state covariance matrix.
* @param {double} delta_t the change in time (in seconds) between the last
* measurement and this one.
*/
void UKF::Predict(const double delta_t) {
  // get predicted states using sigma points:
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    Xsig_pred_.col(i) = DoPredict(Xsig_.col(i), delta_t);
  }

  // predicted state mean:
  x_ = Xsig_pred_ * weights_;

  // deviations from mean:
  Xsig_pred_dev_ = Xsig_pred_.colwise() - x_;

  // predicted state covariance matrix:
  P_xx_ = Xsig_pred_dev_ * (Xsig_pred_dev_.transpose().array().colwise() * weights_.array()).matrix();
}

/**
 * Make observation prediction using lidar measurement model
 * @param {const VectorXd &} predicted k + 1 state using k state
 */
VectorXd UKF::DoUpdateLidar(const VectorXd &x) {
  // allocate:
  VectorXd z_pred(2);

  z_pred = x.head(2);

  return z_pred;
}

/**
* Updates the state and the state covariance matrix using a laser measurement.
* @param {const VectorXd &} z
*/
void UKF::UpdateLidar(const VectorXd &z) {
  // Get predicted observations:
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    Zsig_laser_pred_.col(i) = DoUpdateLidar(Xsig_pred_.col(i));
  }

  // predicted observation mean:
  z_laser_ = Zsig_laser_pred_ * weights_;

  // deviation from mean:
  Zsig_laser_pred_dev_ = Zsig_laser_pred_.colwise() - z_laser_;

  // predicted observation covariance matrix:
  P_zz_laser_ = Zsig_laser_pred_dev_ * (Zsig_laser_pred_dev_.transpose().array().colwise() * weights_.array()).matrix() + R_laser_;

  // predicted state-observation correlation matrix:
  P_xz_laser_ = Xsig_pred_dev_ * (Zsig_laser_pred_dev_.transpose().array().colwise() * weights_.array()).matrix();

  // kalman gain:
  K_laser_ = P_xz_laser_ * P_zz_laser_.inverse();

  // update state:
  x_ = x_ + K_laser_ * (z - z_laser_);
  // update state covariance matrix:
  P_xx_ = P_xx_ - K_laser_ * P_zz_laser_ * K_laser_.transpose();
}

/**
 * Make observation prediction using radar measurement model
 * @param {const VectorXd &} predicted k + 1 state using k state
 */
VectorXd UKF::DoUpdateRadar(const VectorXd &x) {
  // parse state:
  const double p_x = x(0);
  const double p_y = x(1);
  const double v = x(2);
  const double psi = x(3);

  // allocate predicted observation:
  VectorXd z_pred(3);

  // predict observation:
  double r = sqrt(pow(p_x, 2) + pow(p_y, 2));
  z_pred << r,
            atan2(p_y, p_x),
            ((p_x*cos(psi) + p_y*sin(psi)) * v) / r;

  return z_pred;
}

/**
* Updates the state and the state covariance matrix using a radar measurement.
* @param {const VectorXd &} z
*/
void UKF::UpdateRadar(const VectorXd &z) {
  // Get predicted observations:
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {
    Zsig_radar_pred_.col(i) = DoUpdateRadar(Xsig_pred_.col(i));
  }

  // predicted observation mean:
  z_radar_ = Zsig_radar_pred_ * weights_;

  // deviation from mean:
  Zsig_radar_pred_dev_ = Zsig_radar_pred_.colwise() - z_radar_;

  // predicted observation covariance matrix:
  P_zz_radar_ = Zsig_radar_pred_dev_ * (Zsig_radar_pred_dev_.transpose().array().colwise() * weights_.array()).matrix() + R_radar_;

  // predicted state-observation correlation matrix:
  P_xz_radar_ = Xsig_pred_dev_ * (Zsig_radar_pred_dev_.transpose().array().colwise() * weights_.array()).matrix();

  // kalman gain:
  K_radar_ = P_xz_radar_ * P_zz_radar_.inverse();

  // update state:
  x_ = x_ + K_radar_ * (z - z_radar_);
  // update state covariance matrix:
  P_xx_ = P_xx_ - K_radar_ * P_zz_radar_ * K_radar_.transpose();
}
