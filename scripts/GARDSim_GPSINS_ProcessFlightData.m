%
%
%  Process GPS-INS from flight sensor data
%
%  Sensors - 
%       GPS - OEMV1
%       IMU - SiIMU04
%       MAG - HMR2300
%       truth - Novatel SPAN (Honeywell HG58)
%
%

datapath = 'data/Flight_Data/log_22Jul09/';
datafile = 'log_22Jul09-clean.mat';
datafile_full_path = strcat(datapath,datafile);

%% UNCOMMENT THIS SECTION TO RE-READ the DATA
% % imufile = 'log_20090722_imu.out';
% % gpsfile = 'log_20090722_gps.out';
% % magfile = 'log_20090722_mag.out';
% % spanfile = 'log_20090722_span.out';
% % 
% % 

% % 
% % imu = dlmread(strcat(datapath,imufile));
% % mag = dlmread(strcat(datapath,magfile));
% % span = dlmread(strcat(datapath,spanfile));
% % gps = GARD_ReadNovatelLogData(strcat(datapath,gpsfile));
% % 
% % % required because there double range data records due to a logger bug
% % for i=1:size(gps.RangeData,2) / 2
% %    gps.RangeData2(i) = gps.RangeData(2*i); 
% % end
% % gps.RangeData = gps.RangeData2;
% % clear gps.RangeData2;
% % 
% % 
% % 
% % save(datafile_full_path);

load(datafile_full_path);

% set number formatting for display
format long g;

% load GPS constants
GPSConstants;
LOCAL_GRAVITY = 9.80;


%% initial conditions
InitialPosition =  [ -0.481198373257095 ;        2.67055752492006  ;        58.7557216892019];
InitialVelocity = [0; 0; 0 ];
InitialAttitude = [0; 0; 100*pi/180];

GPSLeverArm = [0.085;0.805;0.960]; % in IMU body coorindates, meters

RM = MeridianRadius(InitialPosition(1));
RP = PrimeRadius(InitialPosition(1));
RMh = RM + InitialPosition(3);
RPh = RP + InitialPosition(3);

ElevationMaskAngle = 7.5;

%% organise data to required format

magdata(:,1) = mag(:,2)*1e-9;
magdata(:,2:4) = mag(:,5:7);
MAG_RATE = 20;


mafilt = ones(40,1) / 40;

magdata2(:,1) = magdata(:,1);
magdata2(:,2) = filter(mafilt,1,magdata(:,2));
magdata2(:,3) = filter(mafilt,1,magdata(:,3));
magdata2(:,4) = filter(mafilt,1,magdata(:,4));



Xsf = 1.0;
Ysf = 1.0;
Xoff = 0;
Yoff = 0;  
   
CompassCalibration = [Xsf,Ysf,Xoff,Yoff];


IMU_RATE = 200;
imudata(:,1) = imu(:,2)*1e-9;
imudata(:,5:7) = imu(:,6:8)*IMU_RATE;
imudata(:,8:10) = imu(:,9:11)*IMU_RATE;
imudata(:,11:12) = imu(:,12:13);



for i=1:size(gps.BestXYZData,2)
    gpsdata(i,1) = gps.BestXYZData(i).rtTimeStamp*1e-9;
    gpsdata(i,2) = gps.BestXYZData(i).GPSWeek;
    gpsdata(i,3) = gps.BestXYZData(i).GPSSec;
    gpsdata(i,6:8) = [gps.BestXYZData(i).dPosX,gps.BestXYZData(i).dPosY,gps.BestXYZData(i).dPosZ];
    gpsdata(i,9:11) = [gps.BestXYZData(i).fPosXSigma,gps.BestXYZData(i).fPosYSigma,gps.BestXYZData(i).fPosZSigma];
    
    gpsdata(i,14:16) = [gps.BestXYZData(i).dVelX,gps.BestXYZData(i).dVelY,gps.BestXYZData(i).dVelZ];
    gpsdata(i,17:19) = [gps.BestXYZData(i).fVelXSigma,gps.BestXYZData(i).fVelYSigma,gps.BestXYZData(i).fVelZSigma];
end
GPS_RATE = 1;

SPAN_RATE = 50;
barodata(:,2) = span(:,7);
barodata(:,1) = span(:,2)*1e-9;
BARO_RATE = SPAN_RATE;

gps.SV_Ephemeris = GARD_GPSEphemStruct_to_Table(gps.GPSEphem);


StartTime = 00;
StopTime = 1500;

%% do loose EKF solution
% Sensor_Params.gyro_beta(1) = 1/500;
% Sensor_Params.gyro_beta(2) = 1/500;
% Sensor_Params.gyro_beta(3) = 1/500;
% Sensor_Params.accel_beta(1) = 1/500;
% Sensor_Params.accel_beta(2) = 1/500;
% Sensor_Params.accel_beta(3) = 1/500;
% Sensor_Params.gyro_Q(1) = (1e-2)^2;
% Sensor_Params.gyro_Q(2) = (1e-2)^2;
% Sensor_Params.gyro_Q(3) = (1e-2)^2;
% Sensor_Params.accel_Q(1) = (1e-2)^2;
% Sensor_Params.accel_Q(2) = (1e-2)^2;
% Sensor_Params.accel_Q(3) = (1e-2)^2;
Sensor_Params.gyro_beta(1) = 1/500;
Sensor_Params.gyro_beta(2) = 1/500;
Sensor_Params.gyro_beta(3) = 1/500;
Sensor_Params.accel_beta(1) = 1/500;
Sensor_Params.accel_beta(2) = 1/500;
Sensor_Params.accel_beta(3) = 1/500;
Sensor_Params.gyro_Q(1) = (0.00125)^2;
Sensor_Params.gyro_Q(2) = (0.00125)^2;
Sensor_Params.gyro_Q(3) = (0.00125)^2;
Sensor_Params.accel_Q(1) = (0.0219)^2;
Sensor_Params.accel_Q(2) = (0.0219)^2;
Sensor_Params.accel_Q(3) = (0.0219)^2;

%results_EKF_Loose = GARD_GPSINS_Loose_EKF(imudata,gpsdata,barodata,magdata,GPSLeverArm,IMU_RATE,GPS_RATE,MAG_RATE,BARO_RATE,Sensor_Params,CompassCalibration);

%% do tight EKF solution
%results_EKF_Tight = GARD_GPSINS_TightlyCoupledEKF(imudata, gps,barodata,magdata, ...
%    IMU_RATE,GPS_RATE,MAG_RATE,BARO_RATE, Sensor_Params, CompassCalibration, StartTime, StopTime, InitialPosition, InitialVelocity,InitialAttitude,GPSLeverArm );

%results_EKF_Tight_error = GARD_GPSINS_TightlyCoupledEKF(imudata, gps,barodata,magdata, ...
%    IMU_RATE,GPS_RATE,MAG_RATE,BARO_RATE, Sensor_Params, CompassCalibration, StartTime, StopTime, InitialPosition, InitialVelocity,InitialAttitude );


%% do stand-alone GPS solution
%GARDSim_LSQSolution;





%% do UKF solution
results_UKF_Tight = GARD_GPSINS_UKF(imudata, gps,barodata,magdata, IMU_RATE,GPS_RATE,MAG_RATE,BARO_RATE, Sensor_Params, CompassCalibration, StartTime, StopTime, InitialPosition, InitialVelocity,InitialAttitude );

% %% plot results
figure(); grid on; hold on;
%plot(imudata(:,1),results_EKF_Loose.phi*180/pi,'r');
plot([StartTime:1/IMU_RATE:StopTime],results_EKF_Tight.phi_q * 180/pi,'b');
plot(span(:,2)*1e-9-87.2,span(:,11),'g');


figure(); grid on; hold on;
%plot(imudata(:,1),results.theta*180/pi,'r');

plot([StartTime:1/IMU_RATE:StopTime],results_EKF_Tight.theta_q * 180/pi,'b');
plot(span(:,2)*1e-9-87.2,span(:,12),'g');

figure(); grid on; hold on;
plot([0:1/IMU_RATE:StopTime],results_EKF_Tight.psi_q * 180/pi,'b');
plot(span(:,2)*1e-9-87.2,span(:,13),'g');


figure(); grid on; hold on;
plot([StartTime:1/IMU_RATE:StopTime],results_EKF_Tight.Vel_NED(:,1),'b');
plot(span(:,2)*1e-9-87.2,span(:,8),'g');



% 
% %% calculate roll error
span_offset = [-1.35 * pi/180; +1.60 * pi/180; 0.0 ];


roll_error_LEKF = zeros(1,length(results_EKF_Loose.GPSTime));
pitch_error_LEKF = zeros(1,length(results_EKF_Loose.GPSTime));
yaw_error_LEKF = zeros(1,length(results_EKF_Loose.GPSTime));

roll_error_TEKF = zeros(1,76380);
pitch_error_TEKF = zeros(1,76380);
yaw_error_TEKF = zeros(1,76380);


disp('Calculating errors...');

if exist('index_i','var')
    CalculateIndex = 0;
else
    CalculateIndex = 1;
end


%% new method to calculate errors

for i = 1:length(results_EKF_Loose.GPSTime)
   
    span_i = find(span(:,4) > results_EKF_Loose.GPSTime(i),1);
    
    diff_time(i) = span(span_i,4) - results_EKF_Loose.GPSTime(i);
    
    roll_error_LEKF(i) = results_EKF_Loose.phi(i) - span(span_i,11) * pi/180;
    pitch_error_LEKF(i) = results_EKF_Loose.theta(i) - span(span_i,12) * pi/180;
    yaw_error_LEKF(i) = results_EKF_Loose.psi(i) - span(span_i,13) * pi/180;
    
    if yaw_error_LEKF(i) > pi
        yaw_error_LEKF(i) = yaw_error_LEKF(i) - 2*pi;
    end
    
    if yaw_error_LEKF(i) < -pi
        yaw_error_LEKF(i) = yaw_error_LEKF(i) + 2*pi;
    end
end

%% old method
for i = 40:75600

    if CalculateIndex == 1
        index = find(imudata(:,1) > span(i,2)*1e-9 - 74.018,1) - 1;
%        index2 = find(results_EKF_Tight.imudatatime(:,1) > span(i,2)*1e-9-74.018,1);
        
        % save indeces for next time
        index_i(i,2) = index2;
        index_i(i,1) = index;
    else
        index = index_i(i,1)-20;
        index2 = index_i(i,2)-20;
        
        if index < 1 index = 1; end
        if index2 < 1 index2 = 1; end
    end
    
    roll_error_LEKF(i) = results_EKF_Loose.phi(index) - span(i,11)*pi/180;
    pitch_error_LEKF(i) = results_EKF_Loose.theta(index) - span(i,12)*pi/180;
    yaw_error_LEKF(i) = results_EKF_Loose.psi(index) - span(i,13)*pi/180;
    
    if yaw_error_LEKF(i) > pi
        yaw_error_LEKF(i) = yaw_error_LEKF(i) - 2*pi;
    end
    
    if yaw_error_LEKF(i) < -pi
        yaw_error_LEKF(i) = yaw_error_LEKF(i) + 2*pi;
    end
    
    
    

    
    roll_error_TEKF(i) = results_EKF_Tight.phi_q(index2) - span(i,11)*pi/180;
    pitch_error_TEKF(i) = results_EKF_Tight.theta_q(index2) - span(i,12)*pi/180;
    yaw_error_TEKF(i) = results_EKF_Tight.psi_q(index2) - span(i,13)*pi/180;
    
    if yaw_error_TEKF(i) > pi
        yaw_error_TEKF(i) = yaw_error_TEKF(i) - 2*pi;
    end
    
    if yaw_error_TEKF(i) < -pi
        yaw_error_TEKF(i) = yaw_error_TEKF(i) + 2*pi;
    end
    
    PosLLH_error_TEKF(i,:) = results_EKF_Tight.Pos_LLH(index2,:) - span(i,5:7) .* [pi/180 pi/180 1];
    VelNED_error_TEKF(i,:) = results_EKF_Tight.Vel_NED(index2,:) - span(i,8:10) ./ [1 1 -1];
    
    EKF_HPL(i) = results_EKF_Tight.EKF_HPL(ceil(index2/IMU_RATE/GPS_RATE));
    EKF_VPL(i) = results_EKF_Tight.EKF_VPL(ceil(index2/IMU_RATE/GPS_RATE));
    UKF_HPL(i) = results_UKF_Tight.UKF_HPL(ceil(index2/IMU_RATE/GPS_RATE));
    UKF_VPL(i) = results_UKF_Tight.UKF_VPL(ceil(index2/IMU_RATE/GPS_RATE));
    
    
    %% UKF index should be same as EKF
    roll_error_TUKF(i) = results_UKF_Tight.phi_q(index2) - span(i,11)*pi/180;
    pitch_error_TUKF(i) = results_UKF_Tight.theta_q(index2) - span(i,12)*pi/180;
    yaw_error_TUKF(i) = results_UKF_Tight.psi_q(index2) - span(i,13)*pi/180;
    
    if yaw_error_TUKF(i) > pi
        yaw_error_TUKF(i) = yaw_error_TUKF(i) - 2*pi;
    end
    
    if yaw_error_TUKF(i) < -pi
        yaw_error_TUKF(i) = yaw_error_TUKF(i) + 2*pi;
    end
    
    PosLLH_error_TUKF(i,:) = results_UKF_Tight.UKF_x_hat_save(index2,1:3) - span(i,5:7) .* [pi/180 pi/180 1];
    VelNED_error_TUKF(i,:) = results_UKF_Tight.UKF_x_hat_save(index2,4:6) - span(i,8:10) ./ [1 1 -1];
    
    
     
end




for k=60:1558
    index3 = find(results_UKF_Tight.imudatatime(:,1) > gpstime(k,1),1,'First');
    PosLLH_error_TUKF_DGPS(k,:) = (results_UKF_Tight.UKF_x_hat_save(index3,1:3) - LSQ_Solution_LLH(k,:)) .* [RM RP*cos(-0.47) 1];
end


cptruthoff = 75;
for k=76:1558
    if gpstime(k,3) ~= cptruth(k-cptruthoff,2) 
        warn('discrepancy');
    end
    index4 = find(results_UKF_Tight.imudatatime(:,1) > gpstime(k,1),1,'First');
    
    PosLLH_error_TUKF_CPDGPS(k,:) = (results_UKF_Tight.UKF_x_hat_save(index4,1:3) - [cptruth(k-cptruthoff,15), cptruth(k-cptruthoff,17), cptruth(k-cptruthoff,19)] ) .* [RMh RPh*cos(InitialPosition(1)) 1];
    PosLLH_error_TEKF_CPDGPS(k,:) = (results_EKF_Tight.Pos_LLH(index4,1:3) - [cptruth(k-cptruthoff,15), cptruth(k-cptruthoff,17), cptruth(k-cptruthoff,19)] ) .* [RMh RPh*cos(InitialPosition(1)) 1];

    
    PosLLH_error_T(k,1) = gpstime(k,1);
    PosLLH_error_T(k,2) = gpstime(k,3);
end


% compare span truth with Grafnav (cptruth)
for k = 1:1499
   % find corresponding span data point
   
end

roll_error_LEKF = roll_error_LEKF - span_offset(1);
pitch_error_LEKF = pitch_error_LEKF - span_offset(2);
yaw_error_LEKF = yaw_error_LEKF - span_offset(3);

roll_error_TEKF = roll_error_TEKF - span_offset(1);
pitch_error_TEKF = pitch_error_TEKF + span_offset(2);
yaw_error_TEKF = yaw_error_TEKF - span_offset(3);

roll_error_TUKF = roll_error_TUKF - span_offset(1);
pitch_error_TUKF = pitch_error_TUKF + span_offset(2);
yaw_error_TUKF = yaw_error_TUKF - span_offset(3);


span_time = span(1:76380,2) * 1e-9 - span_toff;

figure(); grid on; hold on;
plot(span_time(1:75600)/60, roll_error_LEKF(1:75600)*180/pi,'b');
plot(span_time(1:75600)/60, roll_error_TEKF(1:75600)*180/pi,'g');
plot(span_time(1:75600)/60, roll_error_TUKF(1:75600)*180/pi,'r');
legend('Loose EKF','Tight EKF','Tight UKF');
xlabel('Test Time (mins)')
ylabel('Error (deg)');
title('Roll Estimate Errors');


figure(); grid on; hold on;
%plot(span_time(1:75600)/60,pitch_error_LEKF(1:75600)*180/pi,'b');
plot(span_time(1:75600)/60,pitch_error_TEKF(1:75600)*180/pi,'g');
plot(span_time(1:75600)/60,pitch_error_TUKF(1:75600)*180/pi,'r');
legend('Loose EKF','Tight EKF','Tight UKF');
xlabel('Test Time (mins)')
ylabel('Error (deg)');
title('Pitch Estimate Errors');

figure(); grid on; hold on;
plot(span_time(1:75600)/60,yaw_error_LEKF(1:75600)*180/pi,'b');
plot(span_time(1:75600)/60,yaw_error_TEKF(1:75600)*180/pi,'g');
plot(span_time(1:75600)/60,yaw_error_TUKF(1:75600)*180/pi,'r');
legend('Loose EKF','Tight EKF','Tight UKF');
xlabel('Test Time (mins)')
ylabel('Error (deg)');
title('Yaw Estimate Errors');

mean(roll_error_LEKF)*180/pi
std(roll_error_LEKF)*180/pi
mean(pitch_error_LEKF)*180/pi
std(pitch_error_LEKF)*180/pi
mean(yaw_error_LEKF)*180/pi
std(yaw_error_LEKF)*180/pi

mean(roll_error_TEKF)*180/pi
std(roll_error_TEKF)*180/pi
mean(pitch_error_TEKF)*180/pi
std(pitch_error_TEKF)*180/pi
mean(yaw_error_TEKF)*180/pi
std(yaw_error_TEKF)*180/pi





mean(roll_error_TUKF)*180/pi
std(roll_error_TUKF)*180/pi
mean(pitch_error_TUKF)*180/pi
std(pitch_error_TUKF)*180/pi
mean(yaw_error_TUKF)*180/pi
std(yaw_error_TUKF)*180/pi


figure(); grid on; hold on;
plot(span_time(1:75600),sqrt(PosLLH_error_TEKF(1:75600,1).^2 + PosLLH_error_TEKF(1:75600,2).^2)*RM,'g');
plot(span_time(1:75600),sqrt(PosLLH_error_TUKF(1:75600,1).^2 + PosLLH_error_TUKF(1:75600,2).^2)*RM,'r');
legend('EKF','UKF');
xlabel('Test Time (mins)')
ylabel('Error (m)');
title('Horizontal Errors');

figure(); grid on; hold on;
plot(span_time(1:75600)/60,PosLLH_error_TEKF(1:75600,3),'g');
plot(span_time(1:75600)/60,PosLLH_error_TUKF(1:75600,3),'r');
legend('EKF','UKF');
xlabel('Test Time (mins)')
ylabel('Error (m)');
title('Vertical Errors');


figure(); hold on; grid on;
plot(results_EKF_Tight.EKF_HPL*RM,'g');
plot(results_UKF_Tight.UKF_HPL*RM,'r');
legend('EKF','UKF');
xlabel('Test Time (mins)')
ylabel('Protection Level (m)');
title('Horizontal Protection Level');


figure(); hold on; grid on;
plot(results_EKF_Tight.EKF_VPL,'g');
plot(results_UKF_Tight.UKF_VPL,'r');
legend('EKF','UKF');
xlabel('Test Time (mins)')
ylabel('Protection Level (m)');
title('Horizontal Protection Level');


%base = span(1,5:7) .* [pi/180 pi/180 1];


% load background images
backgroundfile = 'data/Flight_Data/log_22Jul09/baf_img.jpg';
ybaf_img = imread(backgroundfile);
ybaf_info = imfinfo(backgroundfile);

base = BaseStation.ApproxPos;

im_scale = 9;
im_shift_x = -5500;
im_shift_y = -2100;

im_xmin = im_shift_x;
im_xmax = ybaf_info.Width * im_scale + im_shift_x;
im_ymin = im_shift_y;
im_ymax = ybaf_info.Height * im_scale + im_shift_y;

figure(); hold on; grid on;
%image([im_xmin im_xmax],[im_ymax im_ymin],ybaf_img);
plot((span(:,6)*pi/180-base(2))*RPh*cos(InitialPosition(1)),(span(:,5)*pi/180-base(1))*RMh,'r','LineWidth',2);
%plot((results_EKF_Tight.Pos_LLH(:,2) - base(2))*RM,(results_EKF_Tight.Pos_LLH(:,1)-base(1))*RP*cos(-0.47),'g');
%plot((results_UKF_Tight.UKF_x_hat_save(:,2)-base(2))*RM,(results_UKF_Tight.UKF_x_hat_save(:,1)-base(1))*RP*cos(-0.47),'r');
plot(0,0,'kd','MarkerFaceColor','k','MarkerSize',5);
xlabel('East Position (m)')
ylabel('North Position (m)');
axis([-2500 3500 -1500 2000]);
axis equal;


% what is this number?
% it must be because in the span logger I forgot to subtract the smistart
% variable?
span_toff = 74.118; 

figure(); hold on; grid on;
plot(span(:,2)*1e-9-span_toff,(span(:,5)*pi/180 - base(1)) * RM,'b');
plot(results_EKF_Tight.imudatatime(:,1),(results_EKF_Tight.Pos_LLH(1:300000,1) - base(1)) * RM,'g');

figure(); hold on; grid on;
plot(span(:,2)*1e-9-span_toff,span(:,11),'b');
plot(results_EKF_Tight.imudatatime(:,1),results_EKF_Tight.phi_q(1:300000,1)*180/pi - span_offset(1)*180/pi,'g');
xlabel('Test Time');
ylabel('Roll Angle (deg)');
legend('SPAN','Experiment');


figure(); hold on; grid on;
plot((span(:,2)*1e-9-span_toff) / 60,span(:,12));
plot(results_EKF_Tight.imudatatime(:,1) / 60,results_EKF_Tight.theta_q(1:300000)*180/pi,'r');



GARD_PlotStanford(sqrt(PosLLH_error_TEKF(4000:75600,1).^2 + PosLLH_error_TEKF(4000:75600,2).^2)*RM, ...
    EKF_HPL(4000:75600)'*RM,40,1,80,'Horizontal Stanford Plot - EKF',1);

figure(); grid on; hold on;
plot(span_time(4000:75600)/60,sqrt(PosLLH_error_TEKF(4000:75600,1).^2 + PosLLH_error_TEKF(4000:75600,2).^2)*RM);
plot(span_time(4000:75600)/60,EKF_HPL(4000:75600)*RM,'r');



for i=1:length(span2)
   span2(i,5:7) = span2(i,5:7) ./ [pi/180, pi/180, 1] ;
end

%LSQ_LLH_Error = LSQ_Solution_LLH(47:1546,:) - span2(1:1500,5:7);

LSQ_LLH_Error = LSQ_Solution_LLH_noDGPS(47:1546,:) - LSQ_Solution_LLH(47:1546,:);

figure(); grid on; hold on;
plot(span2(200:1500,2)*1e-9/60,sqrt(LSQ_LLH_Error(200:1500,1).^2 + LSQ_LLH_Error(200:1500,2).^2)*RM,'b');
plot(span2(200:1500,2)*1e-9/60,LSQ_RAIM_HPL(200:1500),'r');
xlabel('Test Time');
ylabel('Error (m)');
title('Least Squares Snapshot Solution');
legend('Horizontal Error','HPL');


figure(); grid on; hold on;
plot(span2(200:1500,2)*1e-9/60,sqrt(LSQ_LLH_Error(200:1500,1).^2 + LSQ_LLH_Error(200:1500,2).^2)*RM,'b');
plot(span2(200:1500,2)*1e-9/60,LSQ_RAIM_HPL(200:1500),'b--','LineWidth',2);
plot(span_time(4000:75600)/60,sqrt(PosLLH_error_TEKF(4000:75600,1).^2 + PosLLH_error_TEKF(4000:75600,2).^2)*RM,'r');
plot(span_time(4000:75600)/60,EKF_HPL(4000:75600)*RM,'r--','LineWidth',2);
legend('Snapshot Error','Snapshot HPL','EKF Error','EKF HPL');
xlabel('Test Time');
ylabel('Error (m)');

figure(); grid on; hold on;
plot(span_time(4000:75600),sqrt(PosLLH_error_TEKF(4000:75600,1).^2 + PosLLH_error_TEKF(4000:75600,2).^2)*RM,'r');
plot(span_time(4000:75600),EKF_HPL(4000:75600)*RM,'r--','LineWidth',2);
legend('EKF Error','EKF HPL');
xlabel('Test Time');
ylabel('Error (m)');

figure(); grid on; hold on;
plot(span2(200:1500,2)*1e-9/60,abs(LSQ_LLH_Error(200:1500,3)),'b');
plot(span2(200:1500,2)*1e-9/60,LSQ_RAIM_VPL(200:1500),'b--','LineWidth',2);
plot(span_time(4000:75600)/60,abs(PosLLH_error_TEKF(4000:75600,3)),'r');
plot(span_time(4000:75600)/60,EKF_VPL(4000:75600),'r--','LineWidth',2);
legend('Snapshot Error','Snapshot VPL','EKF Error','EKF VPL');
xlabel('Test Time');
ylabel('Error (m)');



figure(); grid on; hold on;
plot((span2(200:1500,2)*1e-9 - span_toff)/60,LSQ_RAIM_HPL(200:1500),'k-.','LineWidth',2);
plot(span_time(4000:75600)/60,EKF_HPL(4000:75600)*RM,'k--','LineWidth',1);
plot(span_time(4000:75600)/60,UKF_HPL(4000:75600)*RM,'k-','LineWidth',1);
line([0 30],[50 50],'LineWidth',2,'LineStyle','-','Color','k');
legend('Snapshot HPL','EKF HPL','UKF HPL','APV HAL');
xlabel('Test Time');
ylabel('Level (m)');
axis([0 30 0 60]);

figure(); grid on; hold on;
plot(span_time(4000:75600)/60,EKF_HPL(4000:75600)*RM,'k--','LineWidth',1);
plot(span_time(4000:75600)/60,UKF_HPL(4000:75600)*RM,'k-','LineWidth',1);
legend('EKF HPL','UKF HPL');
xlabel('Test Time');
ylabel('Level (m)');
axis([0 30 12 15]);


figure(); grid on; hold on;
plot((span2(200:1500,2)*1e-9 - span_toff)/60,LSQ_RAIM_VPL(200:1500),'k-.','LineWidth',2);
plot(span_time(4000:75600)/60,EKF_VPL(4000:75600),'k--','LineWidth',1);
plot(span_time(4000:75600)/60,UKF_VPL(4000:75600),'k-','LineWidth',1);
line([0 30],[50 50],'LineWidth',2,'LineStyle','-','Color','k');
line([0 30],[20 20],'LineWidth',1.5,'LineStyle','--','Color','k');
legend('Snapshot VPL','EKF VPL','UKF VPL','APV-1 VAL','APV-2 VAL');
xlabel('Test Time');
ylabel('Level (m)');
axis([0 30 0 100]);

figure(); grid on; hold on;
plot(span_time(4000:75600)/60,EKF_VPL(4000:75600),'k--','LineWidth',1);
plot(span_time(4000:75600)/60,UKF_VPL(4000:75600),'k-','LineWidth',1);
line([0 30],[20 20],'LineWidth',2,'LineStyle','--','Color','k');
legend('EKF VPL','UKF VPL','APV-2 VAL');
xlabel('Test Time');
ylabel('Level (m)');
axis([0 30 15 25]);


figure(); grid on; hold on;
plot(span_time(4000:75600)/60,span(4000:75600,7)*3.28-139);
xlabel('Test Time (mins)');
ylabel('Height (ft-AMSL)');
axis([0 30 0 1100]);


% figure of UKF HPL vs Error
figure(); grid on; hold on;
plot(PosLLH_error_T(76:1558)/60,sqrt(PosLLH_error_TEKF_CPDGPS(76:1558,1).^2 + PosLLH_error_TEKF_CPDGPS(76:1558,2).^2),'k');
plot(PosLLH_error_T(76:1558)/60,results_UKF_Tight.UKF_HPL(2:1484)*RMh,'k--');
%line([0 30],[50 50],'LineWidth',2,'LineStyle','-','Color','k');
line([0 30],[50 50],'LineWidth',1.5,'LineStyle','--','Color','k');
legend('UKF HPL','UKF Error','APV HAL');
xlabel('Test Time');
ylabel('Level (m)');
axis([0 30 0 55]);


figure(); grid on; hold on;
plot(PosLLH_error_T(126:1558)/60,abs(PosLLH_error_TEKF_CPDGPS(126:1558,3)),'k');
plot(PosLLH_error_T(126:1558)/60,results_UKF_Tight.UKF_VPL(52:1484),'k--');
line([0 30],[50 50],'LineWidth',2,'LineStyle','-','Color','k');
line([0 30],[20 20],'LineWidth',1.5,'LineStyle','--','Color','k');

legend('UKF VPL','UKF Error','APV VAL','APV-2 VAL');
xlabel('Test Time');
ylabel('Level (m)');
axis([0 30 0 55]);

% save data
% disp('Saving data...');
% save(datapath);
