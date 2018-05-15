
% load data
if ~exist('DMU') 
    disp('Loading DMU Data...');
    DMU.filename = 'data\DMUData.txt';
    data = load(DMU.filename);
    DMU.sec = data(:,1);
    DMU.nsec = data(:,2);
    DMU.roll = data(:,3);
    DMU.pitch = data(:,4);
    DMU.yaw = data(:,5);
    DMU.wx = data(:,6);
    DMU.wy = data(:,7);
    DMU.wz = data(:,8);
    DMU.ax = data(:,9);
    DMU.ay = data(:,10);
    DMU.az = data(:,11);
    DMU.temp = data(:,12);
    DMU.t= data(:,13);

    clear data;
end

if ~exist('GPS_PXYA')
    disp('Loading GPS Data...');
    load('data/GPSData2.mat'); 
end



dmu_time = DMU.sec + DMU.nsec * 1e-9;
dmu_time = dmu_time - dmu_time(1);

dmu_dt = 0.025;  %% 40 hz - seeplot(diff(dmu_time)) and take an average

NumberEpochsGPS = length(GPS_PXYA);
for(index=1:NumberEpochsGPS)
    
    gpsxyz_time(index) = GPS_PXYA(index).SystemHour * 3600 + ...
                         GPS_PXYA(index).SystemMinute * 60 + ...
                         GPS_PXYA(index).SystemSecond + ...
                         GPS_PXYA(index).SystemNanosecond*1e-9;
end

gps_dt = 0.1;

gpsxyz_time = gpsxyz_time - gpsxyz_time(1);



%% convert GPS data to LLH
disp('Converting GPS to LLH...');
GPS_LLH = zeros(NumberEpochsGPS,3);
for(index=1:NumberEpochsGPS)
   [GPS_LLH(index,1) GPS_LLH(index,2) GPS_LLH(index,3)] = ECEF2LLH([GPS_PXYA(index).XCoordinate;GPS_PXYA(index).YCoordinate;GPS_PXYA(index).ZCoordinate]);
end



%% convert GPS velocities
disp('Calculating GPS Velocities');
NumberEpochsGPSV = length(GPS_VLHA);
for(index=1:NumberEpochsGPSV)
    
    gpsvel_time(index) = GPS_VLHA(index).SystemHour * 3600 + ...
                         GPS_VLHA(index).SystemMinute * 60 + ...
                         GPS_VLHA(index).SystemSecond + ...
                         GPS_VLHA(index).SystemNanosecond*1e-9;
end



gpsvel_time = gpsvel_time - gpsvel_time(1);

GPS_VXYZ = zeros(NumberEpochsGPSV,3);
for(index=1:NumberEpochsGPSV)
   GPS_VXYZ(index,1) = GPS_VLHA(index).HorizontalSpeed * cos(GPS_VLHA(index).TrackOverGround*pi/180);
   GPS_VXYZ(index,2) = GPS_VLHA(index).HorizontalSpeed * sin(GPS_VLHA(index).TrackOverGround*pi/180);
   GPS_VXYZ(index,3) = -GPS_VLHA(index).VerticalSpeed;
   
end

%% calculate GPS accelerations - dogey but works ok
gpsacc_time = gpsvel_time;
GPS_AXYZ(:,1) = [0 diff(GPS_VXYZ(:,1))']/gps_dt;
GPS_AXYZ(:,2) = [0 diff(GPS_VXYZ(:,2))']/gps_dt;
GPS_AXYZ(:,3) = [0 diff(GPS_VXYZ(:,3))']/gps_dt;


%% calculate gps pseudo roll and pitch


g = [0;0;-9.79];
disp('Calculating GPS Pseudoroll');
pseudo_roll = zeros(1,NumberEpochsGPSV);
gps_flight_path_angle = zeros(1,NumberEpochsGPSV);

for(index = 1:NumberEpochsGPSV)
    % calculate flight path angle
   gps_flight_path_angle(index) = atan2(-GPS_VXYZ(index,3),sqrt(GPS_VXYZ(index,1)^2+GPS_VXYZ(index,2)^2)); 
    
   a_tilde(index,:) = GPS_AXYZ(index,:)' - ((GPS_AXYZ(index,:)*GPS_VXYZ(index,:)'/norm(GPS_VXYZ(index,:))^2) * GPS_VXYZ(index,:)');
   g_tilde(index,:) =  g - (g'*GPS_VXYZ(index,:)'/norm(GPS_VXYZ(index,:))^2)'*GPS_VXYZ(index,:)';
   
   pseudo_lift(index,:) = a_tilde(index,:) - g_tilde(index,:);
   p_tilde(index,:) = cross(g_tilde(index,:),GPS_VXYZ(index,:));
   
   pseudo_roll(index) = asin(pseudo_lift(index,:)*p_tilde(index,:)' / (norm(pseudo_lift(index,:)) * norm(p_tilde(index,:))));
   
end

%% filter the buggery out of pseudoroll
filtered_proll = filter(ones(1,30)/30,1,pseudo_roll);


%% plot pseudoattitude results
% figure();
% plot(gpsvel_time,pseudo_roll*180/pi)
% hold on;
% plot(gpsvel_time,filtered_proll*180/pi,'g');
% plot(dmu_time,DMU.roll*180/pi,'r');
% 
% figure();
% plot(gpsvel_time,gps_flight_path_angle*180/pi)
% hold on;
% plot(dmu_time,DMU.pitch*180/pi,'r');
% NumberEpochs = length(dmu_time);



phi_q = zeros(1,NumberEpochs);
theta_q = zeros(1,NumberEpochs);
psi_q = zeros(1,NumberEpochs);
gvec_phi = zeros(1,NumberEpochs);
gvec_theta = zeros(1,NumberEpochs);
phi_err = zeros(1,NumberEpochs);
phi_err2 = zeros(1,NumberEpochs);
theta_err = zeros(1,NumberEpochs);
for(Epoch_hi = 1:180000)

        % get sensor measurements
        omega_x = DMU.wx(Epoch_hi);
        omega_y = DMU.wy(Epoch_hi);
        omega_z = DMU.wz(Epoch_hi);
        A_xb = DMU.ax(Epoch_hi);
        A_yb = DMU.ay(Epoch_hi);
        A_zb = DMU.az(Epoch_hi);

        
        CurrentGPSEpoch = find(gpsvel_time > dmu_time(Epoch_hi),1);
        
        if(Epoch_hi == 1)


            phi_q(Epoch_hi) = 0;
            theta_q(Epoch_hi) = 0;
            psi_q(Epoch_hi) = 0;

            
            
            
        else
        

            phi_dot = (omega_y*sin(phi_q(Epoch_hi-1)) + omega_z * cos(phi_q(Epoch_hi-1)))*tan(theta_q(Epoch_hi-1)) + omega_x;
            theta_dot = omega_y*cos(phi_q(Epoch_hi-1)) - omega_z * sin(phi_q(Epoch_hi-1));
            psi_dot = (omega_y*sin(phi_q(Epoch_hi-1)) + omega_z * cos(phi_q(Epoch_hi-1))) * sec(theta_q(Epoch_hi-1));
            
            phi_q(Epoch_hi) = phi_q(Epoch_hi-1) + phi_dot * dmu_dt;
            theta_q(Epoch_hi) = theta_q(Epoch_hi-1) + theta_dot * dmu_dt;
            psi_q(Epoch_hi) = psi_q(Epoch_hi-1) + psi_dot * dmu_dt;
            
       
        end

        % find gravity vector for attitude estimate
        gvec_phi(Epoch_hi) =  atan2(A_yb,sqrt(A_xb^2 + A_zb^2));% roll
        gvec_theta(Epoch_hi) = atan2(-A_xb,A_zb); % pitch

        phi_err(Epoch_hi) = (phi_q(Epoch_hi) - gvec_phi(Epoch_hi));
        phi_err2(Epoch_hi) = phi_q(Epoch_hi) + filtered_proll(CurrentGPSEpoch);
        
        phi_q(Epoch_hi) = phi_q(Epoch_hi) - phi_err(Epoch_hi)*0.0001;
        
        theta_err(Epoch_hi) = (theta_q(Epoch_hi) - gvec_theta(Epoch_hi));
        theta_q(Epoch_hi) = theta_q(Epoch_hi) - theta_err(Epoch_hi)*0.01;
        

        if(mod(Epoch_hi,1000) == 0)
            disp(sprintf('Completed Epoch %d',Epoch_hi));
        end
end



figure();
subplot(3,1,1),plot(dmu_time,phi_q*180/pi);
hold on;
subplot(3,1,1),plot(dmu_time,DMU.roll*180/pi,'r');
grid on;
ylabel('Roll (deg)');
subplot(3,1,2),plot(dmu_time,theta_q*180/pi);
hold on;
subplot(3,1,2),plot(dmu_time,DMU.pitch*180/pi,'r');
grid on;
ylabel('Pitch (deg)');
subplot(3,1,3),plot(dmu_time,psi_q*180/pi);
hold on;
subplot(3,1,3),plot(dmu_time,DMU.yaw*180/pi,'r');
grid on;
ylabel('Yaw (deg)');

% plot position
figure();
subplot(3,1,1),plot(gpsxyz_time,GPS_LLH(:,1)*180/pi);
ylabel('Latitude (deg)');
grid on;
subplot(3,1,2),plot(gpsxyz_time,GPS_LLH(:,2)*180/pi);
ylabel('Longitude (deg)');
grid on;
subplot(3,1,3),plot(gpsxyz_time,GPS_LLH(:,3));
ylabel('Height MSL');
grid on;




% figure();
% plot(dmu_time,gvec_theta*180/pi);
% hold on;
% plot(dmu_time,DMU.pitch*180/pi,'r');

