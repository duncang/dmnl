%% wind estimator

clear;

%% generate truth data

% wind truth - constant
Vw_truth = 5; %m/s
Aw_truth = 110 * pi/180; % wind direction (to (i.e. opposite to forecast)) in radians

% vehicle truth - zig zag
Va_truth = 50; %m/s

% north, then east, then south-east
for t=1:40
    Aa_truth(t) = 0;
end
for t=41:70
    Aa_truth(t) = 90 * pi/180;
end
for t=71:100
    Aa_truth(t) = 135 * pi/180;
end

Pg_N_truth(1) = 0;
Pg_E_truth(1) = 0;

for t=1:100
    Va_N_truth(t) = Va_truth * cos (Aa_truth(t));
    Va_E_truth(t) = Va_truth * sin (Aa_truth(t));

    Vw_N_truth(t) = Vw_truth * cos (Aw_truth);
    Vw_E_truth(t) = Vw_truth * sin (Aw_truth);
    
    Vg_N_truth(t) = Va_N_truth(t) + Vw_N_truth(t);
    Vg_E_truth(t) = Va_E_truth(t) + Vw_E_truth(t);


    Pg_N_truth(t+1) = Pg_N_truth(t) + Vg_N_truth(t);
    Pg_E_truth(t+1) = Pg_E_truth(t) + Vg_E_truth(t);

end



%% simulate measurements
for t=1:100
   Vg_N_meas(t) = Vg_N_truth(t);% + randn(1);
   Vg_E_meas(t) = Vg_E_truth(t);% + randn(1);
end

%% setup state transition matrix

% PHI = [0,0,1,0,1,0 ;
%        0,0,0,1,0,1 ;
%        1,0,0,0,-1,0 ;
%        0,1,0,0,0,-1 ;
%        1,0,-1,0,0,0 ;
%        0,1,0,-1,0,0 ];

PHI = [0,0,1,0,1,0 ;
       0,0,0,1,0,1 ;
       0,0,1,0,0,0 ;
       0,0,0,1,0,0 ;
       0,0,0,0,1,0 ;
       0,0,0,0,0,1 ];


%% setup measurements matrix
H = [1,0,0,0,0,0 ;
     0,1,0,0,0,0 ];

 
Q = zeros(6,6);
Q(1,1) = 1;
Q(2,2) = 1;
Q(3,3) = 1;
Q(4,4) = 1;
Q(5,5) = 0.1;
Q(6,6) = 0.1;

R = zeros(2,2);
R(1,1) = 1;
R(2,2) = 1;

%% run filter
% setup state vector

x(:,1) = [50;0;50;0;0;0];
P = eye(6,6);
Px_(:,1) = diag(P);
for t=1:100
    
    % state propagation
    x(:,t+1) = PHI * x(:,t); 
    P = PHI*P*PHI' + Q;
    
    % setup measurements vector
    z(:,t+1) = H * x(:,t+1);
    z_meas = [Vg_N_meas(t);
              Vg_E_meas(t)];
    
    % generate innovation
    v(:,t+1) = z_meas - z(:,t+1);
    
    S = H*P*H' + R;
    
    K = P*H' * inv(S);
    
    x(:,t+1) = x(:,t+1) + K * v(:,t+1);
    P = (eye(6,6) - K*H) * P;
    
    P_x(:,t+1) = diag(P);


    windspeed(t) = sqrt(x(5,t)^2 + x(6,t)^2);
    winddirection(t) = atan2(x(6,t),x(5,t));
end

%% plot results

figure; hold on; grid on;
plot(x(5,:),'r');
plot(Vw_N_truth,'g');
xlabel('Time');
ylabel('Wind Velocity North');

figure; hold on; grid on;
plot(x(6,:),'r');
plot(Vw_E_truth,'g');
xlabel('Time');
ylabel('Wind Velocity East');



figure; hold on; grid on;
plot(x(1,:),'r');
plot(Vg_N_truth,'g');
xlabel('Time');
ylabel('Ground Velocity North');

figure; hold on; grid on;
plot(x(2,:),'r');
plot(Vg_E_truth,'g');
xlabel('Time');
ylabel('Ground Velocity East');


figure; hold on; grid on;
plot(x(3,:),'r');
plot(Va_N_truth,'g');
xlabel('Time');
ylabel('Air Velocity North');

figure; hold on; grid on;
plot(x(4,:),'r');
plot(Va_E_truth,'g');
xlabel('Time');
ylabel('Air Velocity East');

