function mnav_data = parsemnavtext(mnav)
% function mnav_data = parsemnavtext(mnav)
% Written by Duncan Greer 18 Jun 2007 
%
% Inputs
% mnav - table of values read from micronav data using dlmread
%
% Outputs
% mnav_data - data struct containing mnav packet data
%
% $Id: parsemnavtext.m 852 2007-06-18 03:26:13Z greerd $
%
% MNAV IMU Packet
% ===============
% 1  - UnixTime
% 2  - yyyy
% 3  - mm
% 4  - dd
% 5  - hh
% 6  - mm
% 7  - ss
% 8  - nanosecond
% 9  - Packet ID
% 10 - X Acceleration
% 11 - Y Acceleration
% 12 - Z Acceleration
% 13 - X Gyro
% 14 - Y Gyro
% 15 - Z Gyro
% 16 - X Mag
% 17 - Y Mag
% 18 - Z Mag
% 19 - X Axis Temperature
% 20 - Y Axis Temperature
% 21 - Z Axis Temperature
% 22 - Static (Baro) Pressure
% 23 - Total (Pitot) Pressure
%
% MNAV GPS PACKET
% ===============
% 1  - UnixTime
% 2  - yyyy
% 3  - mm
% 4  - dd
% 5  - hh
% 6  - mm
% 7  - ss
% 8  - nanosecond
% 9  - Packet ID
% 10 - North Velocity (cm/sec)
% 11 - East Velocity (cm/sec)
% 12 - Down Velocity (cm/sec)
% 13 - Longitude
% 14 - Latitude
% 15 - Altitude
% 16 - ITOW
%



imupacket = 0;
gpspacket = 0;

lengthmnav = size(mnav,1);
for i=1:lengthmnav
     
    if i == 1
       mnav_imu = zeros(lengthmnav,size(mnav,2));
       mnav_gps = zeros(lengthmnav/25,size(mnav,2));
    end
   if mnav(i,9) == 1 ||mnav(i,9) == 3
       imupacket = imupacket + 1;
       mnav_imu(imupacket,:) = mnav(i,:);
   end
   
   if mnav(i,9) == 5
       gpspacket = gpspacket+1;
       mnav_gps(gpspacket,:) = mnav(i,:);
   end
   

    disp(sprintf('Processed %d packets: %d imu, %d gps',i,imupacket,gpspacket));
    
end


mnav_data.imu = mnav_imu(1:imupacket,:);
mnav_data.imu_count = imupacket;
mnav_data.gps = mnav_gps(1:gpspacket,:);
mnav_data.gps_count = gpspacket;
