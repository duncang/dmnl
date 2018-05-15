% script to organise the gps and mnav data into the format expected by the
% GPS INS Algorithm
% Written by Duncan Greer 18 Jun 2007
% $Id: OrganiseDCSData.m 907 2007-11-19 00:06:55Z greerd $
% 
% Instructions
% ============
% 1. Load data using DLMREAD and then run parsemnavtext.m
% 2. Run This Script
% 3. Run results = GARD_GPSINS_Loose_EKF(imudata,gpsdata2,barodata,magdata,IMU_RATE,1,MAG_RATE,BARO_RATE)
% 4. Save results.

% setup data rates
IMU_RATE = 100;
BARO_RATE = 100;
MAG_RATE = 100;
GPS_RATE = 20;

% micronav data splits into imu, mag and baro data
imudata(:,1) = mnav.imu(:,1) + mnav.imu(:,8)*1e-9;
imudata(:,5:7) = mnav.imu(:,10:12)*9.80; % acceleartion % convert from Gs to m/s/s


imudata(:,8:10) = mnav.imu(:,13:15)*pi/180/pi; % gyros - convert to radians




for i=1:length(mnav.imu(:,1));
    [imudata(i,11),imudata(i,12)] = GARD_GravityVector(imudata(i,5:7));
end

% baro data
barodata(:,1) = imudata(:,1);
barodata(:,2) = mnav.imu(:,22);


magdata(:,1) = imudata(:,1);
magdata(:,2:4) = mnav.imu(:,16:18);

% gps data
gpsdata(:,1) = NovatelData.GPSTime_Sec + ppsdata(1,1)+ppsdata(1,8)*1e-9 - NovatelData.GPSTime_Sec(1);
gpsdata(:,6:8) = NovatelData.PosECEF;
gpsdata(:,9:11) = NovatelData.PosECEF_Sigma;
gpsdata(:,14:16) = NovatelData.VelECEF;
gpsdata(:,17:19) = NovatelData.VelECEF_Sigma;

% convert gps data to 1hz
gpsdata2 = downsample(gpsdata,GPS_RATE);




