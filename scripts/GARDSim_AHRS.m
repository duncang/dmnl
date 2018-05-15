

vehicle.C_BN = eye(3);

vehicle.phi = zeros(length(imu),1);
vehicle.theta = zeros(length(imu),1);
vehicle.psi = zeros(length(imu),1);
vehicle.time = zeros(length(imu),1);
vehicle.rtTimestamp = zeros(length(imu),1);
GyroBias = zeros(ceil(length(imu)/200),3);

IMU_RATE = 200;
ins_dt = 1/IMU_RATE;

NumberMeasurements = 2;
NumberStates = 6;

Hk = zeros(NumberMeasurements,NumberStates);
PHIk = eye(NumberStates,NumberStates);
Pk = eye(NumberStates,NumberStates);
Rk = 10*eye(NumberMeasurements,NumberMeasurements);
Qk = zeros(NumberStates,NumberStates);
j = 0;

A_xb = zeros(length(imu),1);
A_yb = zeros(length(imu),1);
A_zb = zeros(length(imu),1);

xk = zeros(6,1);

filter_len=2000;


% initialise attitude
vehicle.C_BN = GARD_EulerToDCM(imu(1,12),imu(1,13),0);


for i=1:length(imu)
    
    if j > 0
    omega_x = imu(i,9)*IMU_RATE + GyroBias(j,1);
    omega_y = imu(i,10)*IMU_RATE + GyroBias(j,2);
    omega_z = imu(i,11)*IMU_RATE + GyroBias(j,3);
    else
    omega_x = imu(i,9)*IMU_RATE;
    omega_y = imu(i,10)*IMU_RATE;
    omega_z = imu(i,11)*IMU_RATE;
    end
    
    
    if i == 1
        
        
        
        
        
    A_xb(i) = imu(i,6);
    A_yb(i) = imu(i,7);
    A_zb(i) = imu(i,8);
    else
    A_xb(i) = (1-1/filter_len)*A_xb(i-1) + (1/filter_len)*imu(i,6);
    A_yb(i) = (1-1/filter_len)*A_yb(i-1) + (1/filter_len)*imu(i,7);
    A_zb(i) = (1-1/filter_len)*A_zb(i-1) + (1/filter_len)*imu(i,8);
    end
    % update time
    if i == 1
        vehicle.time(i) = 0;
    else
        vehicle.time(i) = vehicle.time(i-1) + ins_dt;
    end

    vehicle.rtTimestamp(i) = imu(i,2) * 1e-9;
    
    % update DCM
    vehicle.C_BN = GARD_DCMUpdate2(vehicle.C_BN,[omega_x;omega_y;omega_z]*ins_dt);

    % update PHI
    PHIk(1:3,4:6) = vehicle.C_BN;
    
    PHIk(4,4) = 1 - (1/30);
    PHIk(5,5) = 1 - (1/30);
    PHIk(6,6) = 1 - (1/30);
    
    Qk(4,4) = 2 * 0.05 / 30;
    Qk(5,5) = 2 * 0.05 / 30;
    Qk(6,6) = 2 * 0.05 / 30;
    
    
    Hk(1,1) = 1;
    Hk(2,2) = 1;
    
    

    % convert to euler angles for later analysis
    [vehicle.phi(i),vehicle.theta(i),vehicle.psi(i)] = GARD_DCMToEuler(vehicle.C_BN);

    
    if mod(i,200) == 0
        j = j+1;
        % calculate the gravity vector
        [vehicle.gvec_phi(j), vehicle.gvec_theta(j)] = GARD_GravityVector(A_xb(i),A_yb(i),A_zb(i));

        % form the measurement
        zk(1,1) = vehicle.gvec_phi(j) - vehicle.phi(i);
        zk(2,1) = vehicle.gvec_theta(j) - vehicle.theta(i);

        % predict
        xk_minus = PHIk * xk;
        Pk_minus = PHIk * Pk * PHIk' + Qk;
        Vk = Hk * Pk_minus * Hk' + Rk;
        Kk = Pk_minus * Hk' * pinv(Vk);
        
        % update
        xk = xk_minus + Kk * zk;
        Pk = Pk_minus - Kk * Hk * Pk_minus;

        xk_save(j,:) = xk;
        P_save(j,:) = diag(Pk);
        
        
        del_alpha = xk(1);
        del_beta = xk(2);
        del_gamma = xk(3);

        %if j > 1
        %GyroBias(j,1) = GyroBias(j-1,1) + xk(4);
        %GyroBias(j,2) = GyroBias(j-1,2) + xk(5);
        %GyroBias(j,3) = GyroBias(j-1,3) + xk(6);
        %else    
        GyroBias(j,1) = xk(4);
        GyroBias(j,2) = xk(5);
        GyroBias(j,3) = xk(6);
        %end
        xk(1:3) = 0;
        
        % correct atttitude - del alpha,beta,gamma are actually euler angle
        % deltas
        del_att_skew = [0         -del_gamma   del_beta; ...
                       del_gamma  0          -del_alpha; ...
                       -del_beta  del_alpha   0];

        vehicle.C_BN = vehicle.C_BN * (eye(3,3) + del_att_skew);

        vehicle.C_BN = GARD_OrthogonaliseDCM(vehicle.C_BN);
        
        % convert to euler angles for later analysis
        [vehicle.phi(i),vehicle.theta(i),vehicle.psi(i)] = GARD_DCMToEuler(vehicle.C_BN);

    end
    
    
end


figure(); grid on; hold on;
plot(span(:,2)*1e-9,span(:,11))
plot(vehicle.rtTimestamp,vehicle.phi*180/pi,'r')
xlabel('Time (sec)');
ylabel('Roll Angle (deg)');

figure(); grid on; hold on;
plot(span(:,2)*1e-9,span(:,12))
plot(vehicle.rtTimestamp,vehicle.theta*180/pi,'r')
xlabel('Time (sec)');
ylabel('Pitch Angle (deg)');

figure(); grid on; hold on;
plot(span(:,2)*1e-9,span(:,13))
plot(vehicle.rtTimestamp,vehicle.psi*180/pi,'r')
xlabel('Time (sec)');
ylabel('Heading Angle (deg)');


