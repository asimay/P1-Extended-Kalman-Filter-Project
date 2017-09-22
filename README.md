# Extended Kalman Filter Project 
Project Purpose:

Utilize a kalman filter to estimate the state of a moving object of interest with noisy lidar and radar measurements. Passing the project requires obtaining RMSE values that are lower that the tolerance outlined in the project rubric. 

---

## 1. Read the Radar and Lidar measurement data into program.

This is provided in main.cpp file.

## 2. Kalman Filter or Extended Kalman Filter Process:

1. Lidar measurement module is linear, so we use Kalman Filter to handle the data.

2. Radar measurement module is nolinear, so we need to transform some nonlinear module or matrix to linear to handle the data, which is Extended Kalman Filter.

3. Kalman Filter is "predict --> update --> predict --> update" loop to handle the data.

4. Initializing the State Vector, use lidar measurement to initialize the state variable locations px, py, use the radar measurements ρ and ϕ to initialize the state variable locations px, py, vx, vy.

5. Calculating KF/EKF equations:

`y = z - H * x' or y = z - h(x);`

S = H_ * P_ * H_.transpose() + R_;

K = P_ * H_.transpose() * S.inverse();`

new estimate:

`x_ = x_ + (K * y);

P_ = (I - K * H_) * P_;`

This is for measurement update. For EKF, I calculated the Hj(Jacobian Matrix) to do measurement update.

predict the state:

x_ = F_ * x_ ;

P_ = F_ * P_ * F_.transpose() + Q_;

This is for predict.

6. Normalizing Angles:

  float theta_offset = y(1);
  // normalize theta
  if(theta_offset < -PI) {
      y(1) = theta_offset + 2 * PI;
  }
  else if (theta_offset > PI) {
      y(1) = theta_offset - 2 * PI;
  }

## Test Implementation:

1. Fusion both Radar and Lidar data:
RMSE screenshot as below:


2. Fusion only Radar data:
RMSE screenshot as below:


3. Fusion only Lidar data:
RMSE screenshot as below:


4. Fusion both Radar and Lidar data, but tune a good Q process noise, I increase the noise_ax and noise_ay  up to 13, and got good result:
RMSE screenshot as below:

