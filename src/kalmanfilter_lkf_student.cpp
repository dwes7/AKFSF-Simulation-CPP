// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Linear Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr bool INIT_ON_FIRST_PREDICTION = true;
constexpr double INIT_POS_STD = 10;
constexpr double INIT_VEL_STD = 10;
constexpr double ACCEL_STD = 0.0;
constexpr double GPS_POS_STD = 3.0;
// -------------------------------------------------- //

void KalmanFilter::predictionStep(double dt)
{
    if (!isInitialised() && INIT_ON_FIRST_PREDICTION)
    {
        // Implement the State Vector and Covariance Matrix Initialisation in the
        // section below if you want to initialise the filter WITHOUT waiting for
        // the first measurement to occur. Make sure you call the setState() /
        // setCovariance() functions once you have generated the initial conditions.
        // Hint: Assume the state vector has the form [X,Y,VX,VY].
        // Hint: You can use the constants: INIT_POS_STD, INIT_VEL_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
            VectorXd state = Vector4d::Zero();
            MatrixXd cov = Matrix4d::Zero();



            // Assume the initial position is (X,Y) = (0,0) m
            // Assume the initial velocity is 5 m/s at 45 degrees (VX,VY) = (5*cos(45deg),5*sin(45deg)) m/s
            
            /**
             * x = [0, 0, 3.536, 3.536]^T
            */
           
            state(0) = 0;
            state(1) = 0;
            // Test initial state
            // state(2) = 5.0*cos(M_PI/4);
            // state(3) = 5.0*sin(M_PI/4);
            // Test covariance estimate
            state(2) = 0;
            state(3) = 0;
            cov(0,0) = INIT_POS_STD*INIT_POS_STD;
            cov(1,1) = INIT_POS_STD*INIT_POS_STD;
            cov(2,2) = INIT_VEL_STD*INIT_VEL_STD;
            cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;


            setState(state.transpose());
            setCovariance(cov);
        // ----------------------------------------------------------------------- //
    }

    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Prediction Step for the system in the  
        // section below.
        // Hint: You can use the constants: ACCEL_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        /**
         *        | 1 0 dt 0 |
         *        | 0 1 0 dt |
         *   F =  | 0 0 1  0 |
         *        | 0 0 0  1 |
        */
        MatrixXd F_process = Matrix4d::Zero();
        F_process(0, 0) = 1.0;
        F_process(0, 2) = dt;
        F_process(1, 1) = 1.0;
        F_process(1, 3) = dt;
        F_process(2, 2) = 1.0;
        F_process(3, 3) = 1.0;

        state = F_process * state;


        MatrixXd Q = Matrix2d::Zero();
        Q(0, 0) = ACCEL_STD*ACCEL_STD;
        Q(1, 1) = ACCEL_STD*ACCEL_STD;


        MatrixXd L = MatrixXd::Zero(4, 2);
        double hdt2 = 0.5*dt*dt;
        L(0, 0) = hdt2;
        L(1, 1) = hdt2;
        L(2, 0) = dt;
        L(3, 1) = dt;

        cov = F_process*cov *F_process.transpose() + L*Q*L.transpose();

        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    }
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Update Step for the GPS Measurements in the 
        // section below.
        // Hint: Assume that the GPS sensor has a 3m (1 sigma) position uncertainty.
        // Hint: You can use the constants: GPS_POS_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE 
        // Measurement model sensitivity matrix
        MatrixXd H = MatrixXd(2, 4);
        H << 1,0,0,0,0,1,0,0;


        // Measurement covarience matrix
        Matrix2d R = Matrix2d::Zero();
        R(0, 0) = GPS_POS_STD*GPS_POS_STD;
        R(1, 1) = GPS_POS_STD*GPS_POS_STD;

        // measurement state
        VectorXd zk = Vector2d::Zero();
        zk(0) = meas.x;
        zk(1) = meas.y;

        // innovation difference between current measurement and predicted sensor measurement
        VectorXd innovation = zk - H * state;;


        // innovation covariance
        // This tells us how much uncertainty in the 
        //2x4*2x2*2*4
        Matrix2d Sk = H * cov * H.transpose() + R;

        // kalman gain matrix
        MatrixXd kalman_gain = cov * H.transpose() * Sk.inverse();

        // state update step
        state = state + kalman_gain*innovation;

        // covariance update step
        cov = (Matrix4d::Identity(4, 4) - kalman_gain * H)*cov;
        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    }
    else
    {
        // Implement the State Vector and Covariance Matrix Initialisation in the
        // section below. Make sure you call the setState/setCovariance functions
        // once you have generated the initial conditions.
        // Hint: Assume the state vector has the form [X,Y,VX,VY].
        // Hint: You can use the constants: GPS_POS_STD, INIT_VEL_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
            VectorXd state = Vector4d::Zero();
            MatrixXd cov = Matrix4d::Zero();

            // state << 0,0,0,0;
            // cov << 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0;

            state(0) = meas.x;
            state(1) = meas.y;
            cov(0,0) = GPS_POS_STD*GPS_POS_STD;
            cov(1,1) = GPS_POS_STD*GPS_POS_STD;
            cov(2,2) = INIT_VEL_STD*INIT_VEL_STD;
            cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;


            setState(state);
            setCovariance(cov);
        // ----------------------------------------------------------------------- //
    }        
}

Matrix2d KalmanFilter::getVehicleStatePositionCovariance()
{
    Matrix2d pos_cov = Matrix2d::Zero();
    MatrixXd cov = getCovariance();
    if (isInitialised() && cov.size() != 0){pos_cov << cov(0,0), cov(0,1), cov(1,0), cov(1,1);}
    return pos_cov;
}

VehicleState KalmanFilter::getVehicleState()
{
    if (isInitialised())
    {
        VectorXd state = getState(); // STATE VECTOR [X,Y,VX,VY]
        double psi = std::atan2(state[3],state[2]);
        double V = std::sqrt(state[2]*state[2] + state[3]*state[3]);
        return VehicleState(state[0],state[1],psi,V);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt){predictionStep(dt);}
void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map){}
void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map){}

