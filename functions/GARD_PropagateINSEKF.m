function out = GARD_PropagateINSEKF(in,NumberStates,ins_dt)
%
%
% Definition of state vector
% in.x
% in.P
% 


% intiially copy in to out
out = in;

% apply sensor bias corrections
in.Acc = in.Acc - in.AccelBias;
in.Omega = in.Omega - in.GyroBias;

% propagate the INS and clock
g = GravityModel(in.x(1:3));
[out.x,output] = GARD_INSMechanisation(in.x(1:10),ins_dt,in.Acc,in.Omega,g);

lat_dot = output.llhdot(1);
long_dot = output.llhdot(2);
h_dot = output.llhdot(3);

q_INS = out.x(7:10);

C_BN = output.C_BN;

phi_q = atan2(C_BN(3,2),C_BN(3,3));
theta_q = asin(-C_BN(3,1));
psi_q = atan2(C_BN(2,1),C_BN(1,1));


out.C_BN = C_BN;
out.euler = [phi_q;theta_q;psi_q];

% clock update
out.UserClock(2) = in.UserClock(2);
out.UserClock(1) = in.UserClock(1) + in.UserClock(2)*ins_dt;

% propagate the state error

% generate PHI matrix
F_INS = zeros(NumberStates,NumberStates);
F_INS(1:9,1:9) = GARD_GenerateINSFMatrix(in.x(1:3),in.x(4:6),C_BN*in.Acc');

F_INS(16,17) = 1.0;

phi2 = eye(NumberStates,NumberStates);
phi2 = expm(F_INS*ins_dt);

% % augmented state matrix
phi2(7:9,13:15) = -C_BN*ins_dt; %% attitude to gyro
phi2(4:6,10:12) = C_BN*ins_dt; %%  velocity to acceleration

phi2(10,10) = exp( - in.x_accel_beta * ins_dt);
phi2(11,11) = exp( - in.y_accel_beta * ins_dt);
phi2(12,12) = exp( - in.z_accel_beta * ins_dt);
phi2(13,13) = exp( - in.x_gyro_beta * ins_dt);
phi2(14,14) = exp( - in.y_gyro_beta * ins_dt);
phi2(15,15) = exp( - in.z_gyro_beta * ins_dt);


phi2(16,17) = 1.0*ins_dt;



%Sf = 8*pi^2*2e-20*Speedoflight^2;  % frequency noise psd
%Sc = 2*2e-19*Speedoflight^2;  % bias noise psd
Sc = 1*0.14;
Sf = 1*0.0359;

Q = zeros(NumberStates,NumberStates);

Q(10,10) = 2*in.x_accel_beta*in.x_accel_Q^2;
Q(11,11) = 2*in.y_accel_beta*in.y_accel_Q^2;
Q(12,12) = 2*in.z_accel_beta*in.z_accel_Q^2;
Q(13,13) = 2*in.x_gyro_beta*in.x_gyro_Q^2;
Q(14,14) = 2*in.y_gyro_beta*in.y_gyro_Q^2;
Q(15,15) = 2*in.z_gyro_beta*in.z_gyro_Q^2;

Q(16,16) = (Sc)*ins_dt + Sf*(ins_dt^3)/3;
Q(17,17) = Sf*ins_dt;
Q(16,17) = Sf*(ins_dt^2)/2;
Q(17,16) = Sf*(ins_dt^2)/2;




G = zeros(NumberStates,NumberStates);
G(7:9,13:15) = -C_BN;  %% attitude to gyro
G(4:6,10:12) = C_BN; %%  velocity to acceleration

G(10:12,10:12) = eye(3,3);
G(13:15,13:15) = eye(3,3);

G(16,16) = 1;
G(16,17) = 1;
G(17,16) = 0;
G(17,17) = 1;

Q2d = phi2 * (G * Q * G') * phi2' * ins_dt;

Q2d(16,16) = 1*Sc*ins_dt+Sf*(ins_dt^3)/3;
Q2d(16,17) = 1*Sf*(ins_dt^2)/2;
Q2d(17,16) = 1*Sf*(ins_dt^2)/2;
Q2d(17,17) = 1*Sf*ins_dt;

%x_hat_minus = phi2 * x_hat_minus;

% propagate the covariance matrix
out.P = phi2 * in.P * phi2' + Q2d;





