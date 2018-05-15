%
%
%  Process GPS-INS from flight sensor data
%
%  Sensors - 
%       GPS - OEMV1
%       IMU - SiIMU04
%       MAG - Nil
%
%



%datapath = 'data/Flight_Data/dsa_31Aug10/flight1/';

%load(datapath);

% set number formatting for display
format long g;

% load GPS constants
GPSConstants;
LOCAL_GRAVITY = 9.80;


%% initial conditions - close enough for calculating earth params
InitialPosition =   [-0.461899368467055          2.64237565630361          429.047781012021];
InitialVelocity = [0; 0; 0 ];
InitialAttitude = [0; 0; 0];


RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));
RMh = RM + InitialPosition(3);
RPh = RP + InitialPosition(3);

ElevationMaskAngle = 7.5;

%% organise data to required format
magdata = 0;
MAG_RATE = 0;

Xsf = 1.0;
Ysf = 1.0;
Xoff = 0;
Yoff = 0;  
   
CompassCalibration = [Xsf,Ysf,Xoff,Yoff];


IMU_RATE = 200;
imudata(:,1) = imu(:,2)*1e-9;
imudata(:,5:7) = imu(:,6:8)*IMU_RATE;  % convert rates
imudata(:,8:10) = imu(:,9:11)*IMU_RATE;
imudata(:,11:12) = imu(:,12:13);


gpsstart = 1;

for i=gpsstart:size(gps.BestXYZData,2)
    gpsdata(i-gpsstart+1,1) = gps.BestXYZData(i).rtTimeStamp*1e-9;
    gpsdata(i-gpsstart+1,2) = gps.BestXYZData(i).GPSWeek;
    gpsdata(i-gpsstart+1,3) = gps.BestXYZData(i).GPSSec;
    gpsdata(i-gpsstart+1,6:8) = [gps.BestXYZData(i).dPosX,gps.BestXYZData(i).dPosY,gps.BestXYZData(i).dPosZ];
    gpsdata(i-gpsstart+1,9:11) = [gps.BestXYZData(i).fPosXSigma,gps.BestXYZData(i).fPosYSigma,gps.BestXYZData(i).fPosZSigma];
    
    gpsdata(i-gpsstart+1,14:16) = [gps.BestXYZData(i).dVelX,gps.BestXYZData(i).dVelY,gps.BestXYZData(i).dVelZ];
    gpsdata(i-gpsstart+1,17:19) = [gps.BestXYZData(i).fVelXSigma,gps.BestXYZData(i).fVelYSigma,gps.BestXYZData(i).fVelZSigma];
end

GPS_RATE = 1;


barodata = 0;
BARO_RATE = 0;

StartTime = 00;
StopTime = size(gps.BestXYZData);


%% IMU paramaters - SiIMU04
Sensor_Params.gyro_beta(1) = 1/300;
Sensor_Params.gyro_beta(2) = 1/300;
Sensor_Params.gyro_beta(3) = 1/300;
Sensor_Params.accel_beta(1) = 1/300;
Sensor_Params.accel_beta(2) = 1/300;
Sensor_Params.accel_beta(3) = 1/300;
Sensor_Params.gyro_Q(1) = (0.00125)^2;
Sensor_Params.gyro_Q(2) = (0.00125)^2;
Sensor_Params.gyro_Q(3) = (0.00125)^2;
Sensor_Params.accel_Q(1) = (0.0219)^2;
Sensor_Params.accel_Q(2) = (0.0219)^2;
Sensor_Params.accel_Q(3) = (0.0219)^2;

% IMU to GPS Antenna in meters; IMU body coordinates
%GPSLeverArm = [-.311;-2.334;-0.233];
GPSLeverArm = [0.085;0.805;0.960];
%GPSLeverArm = [0;0;0];

%% do loose EKF solution
results_EKF_Loose = GARD_GPSINS_Loose_EKF(imudata,gpsdata,barodata,magdata,GPSLeverArm, IMU_RATE,GPS_RATE,MAG_RATE,BARO_RATE,Sensor_Params,CompassCalibration);


%% do stand-alone GPS solution
%GARDSim_LSQSolution;

if ~exist('index_i','var')
    for i = 1:size(span,1)
        index = find(results_EKF_Loose.insupdatetime > span(i,2)*1e-9 - 74.018,1) - 1;
        if isempty(index)
            break;
        end
        % save indeces for next time
        index_i(i,1) = index;
    end
end

roll_error_LEKF = zeros(size(span,1),1);
pitch_error_LEKF = zeros(size(span,1),1);
yaw_error_LEKF = zeros(size(span,1),1);

for i = 1:size(span,1)

    if index_i(i) == 0
        continue;
    end
    
    
    roll_error_LEKF(i) = results_EKF_Loose.phi(index_i(i)) - span(i,11)*pi/180;
    pitch_error_LEKF(i) = results_EKF_Loose.theta(index_i(i)) - span(i,12)*pi/180;
    yaw_error_LEKF(i) = results_EKF_Loose.psi(index_i(i)) - span(i,13)*pi/180;
    
    if yaw_error_LEKF(i) > pi
        yaw_error_LEKF(i) = yaw_error_LEKF(i) - 2*pi;
    end
    
    if yaw_error_LEKF(i) < -pi
        yaw_error_LEKF(i) = yaw_error_LEKF(i) + 2*pi;
    end
    
end


