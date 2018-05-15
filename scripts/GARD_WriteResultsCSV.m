%% write log data - assumes data to write is loaded to the workspace

% results_EKF_Loose = 
% 
%       GPSPos_E: [1886x3 double]
%      GPSPos_EE: [1886x3 double]
%     GPSVel_NED: [1886x3 double]
%      GPSVel_EE: [1886x3 double]
%        Pos_LLH: [377182x3 double]
%        Pos_NED: [377182x3 double]
%        Acc_NED: [377182x3 double]
%        Vel_NED: [377182x3 double]
%        Vel_UVW: [377182x3 double]
%            phi: [377182x1 double]
%          theta: [377182x1 double]
%            psi: [377182x1 double]
%          x_hat: [1886x15 double]
%          P_hat: [1886x15 double]
%              z: [1885x8 double]
%           C_BN: [3x3 double]
%        ObsRank: [1x1885 double]
%       GyroBias: [1885x3 double]
%      AccelBias: [1885x3 double]
%      PseudoAtt: [1886x3 double]


% for each trigger pulse, get the state vector

for i=1:length(trigger)

    imu_index = find(imu(:,2) > trigger(i,2),1);
    
    
    % calculate how old the inertial data is for this trigger
    imu_t_off = (imu(imu_index,2) - trigger(i,2)) * 1e-9;
    
    
    gps_index = find(gpsdata(:,1) > trigger(i,2)*1e-9,1); % the next GPS data available (1 sec epochs)
    
    if ~isempty(gps_index)

        gps_late(i) = gpsdata(gps_index,1) - trigger(i,2) * 1e-9;  % how much later the GPS data is (seconds)
    else
        % no more data - punch out
        
        break;
    end
    
    outmatrix(i,1) = trigger(i,2)*1e-9;  % the trigger time
    outmatrix(i,2) = imu(imu_index,2)*1e-9; % the corresponding IMU time
    outmatrix(i,3) = trigger(i,3) + trigger(i,4)*1e-6; % the unix time of the trigger pulse
    
    
    outmatrix(i,4) = gpsdata(gps_index,3) - gps_late(i); % the GPS second-of-week (approximate) of the trigger pulse
    
    
    outmatrix(i,5:7) = results_EKF_Loose.Pos_LLH(imu_index,:);
    outmatrix(i,8:10) = results_EKF_Loose.Vel_NED(imu_index,:);
    outmatrix(i,11) = results_EKF_Loose.phi(imu_index);
    outmatrix(i,12) = results_EKF_Loose.theta(imu_index);
    outmatrix(i,13) = results_EKF_Loose.psi(imu_index);

end;


outfile = 'processed.txt';

dlmwrite(outfile,outmatrix,'delimiter',',','precision',20);

