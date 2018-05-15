function [x_out,output] = GARD_INSMechanisation(x_in,ins_dt,Acc_in,Omega_in,g)
% function [x_out,output] = GARD_INSMechanisation(x_in,ins_dt,Acc_in,Omega_in,g)
% INS Mechanisation Function Written by Duncan Greer 18 Dec 2007
% 
% Input State Vector
% ==================
% x1  = Latitude position (rad)
% x2  = Longitude position (rad)
% x3  = Height position (m)
% x4  = North velocity (m/s)
% x5  = East velocity (m/s)
% x6  = Down velicty (m/s)
% x7  = q0
% x8  = q1
% x9  = q2
% x10 = q3
% x11 = x acc bias (m/s/s)
% x12 = y acc bias (m/s/s)
% x13 = z acc bias (m/s/s)
% x14 = x gyro bias (m/s/s)
% x15 = y gyro bias (m/s/s)
% x16 = z gyro bias (m/s/s)
%
% Outputs 
% =======
% x_out = output state vector (10x1)
% outputs = data structure containing DCM and euler angles
%
% $Id: GARD_INSMechanisation.m 1850 2008-07-14 04:52:47Z greerd $
%

xs_0 = x_in;
RM = MeridianRadius(xs_0(1));
RP = PrimeRadius(xs_0(1));
RMh = RM + xs_0(3);
RPh = RP + xs_0(3);
OMEGAedot = 7.2921151467e-005;

A_xb = Acc_in(1);
A_yb = Acc_in(2);
A_zb = Acc_in(3);

omega_x = Omega_in(1);
omega_y = Omega_in(2);
omega_z = Omega_in(3);


% correct for bias estimate
% A_xb = A_xb - xs_0(11);
% A_yb = A_yb - xs_0(12);
% A_zb = A_zb - xs_0(13);
% omega_x = omega_x - xs_0(14);
% omega_y = omega_y - xs_0(15);
% omega_z = omega_z - xs_0(16);


%xs_0_k(7:10) = quatupdate(xs_0(7:10),omega_x*ins_dt,omega_y*ins_dt,omega_z*ins_dt);

% calculate DCM from updated quaternion
%C_BN = GARDSim_DCMfromQuat(xs_0_k(7:10));



C_BN = GARD_QuatToDCM(xs_0(7:10));
C_BN = GARD_DCMUpdate(C_BN,omega_x*ins_dt,omega_y*ins_dt,omega_z*ins_dt,ins_dt);

if ~isreal(C_BN)
%   warning('Imaginary C_BN');
   C_BN = real(C_BN);
end

% convert to euler angles for later analysis
phi_q = atan2(C_BN(3,2), C_BN(3,3));
theta_q = asin(-C_BN(3,1));
psi_q = atan2(C_BN(2,1), C_BN(1,1));

xs_0_k(7:10) = EulerToQuat([phi_q,theta_q,psi_q]);

% calculate rotation matrices and velocity vector used to calculate coriolis
Tecef2ned= T_ECEF2NED(xs_0(1),xs_0(2));  
Tned2ecef = Tecef2ned';

% rotate accelerometer measurements to nav frame using A_n = C_BN *
% A_b
A_n = C_BN * [A_xb;A_yb;A_zb];
A_xn = A_n(1);
A_yn = A_n(2);
A_zn = A_n(3) - g;


lat = xs_0(1);
long = xs_0(2);
hgt = xs_0(3);
Vn = xs_0(4);
Ve = xs_0(5);
Vd = xs_0(6);



% calculate coriolis correction
% CorMat(1,1) = 0;
% CorMat(1,2) = (2*OMEGAedot + long_dot)*sin(xs_0(1));
% CorMat(1,3) = -lat_dot;
% CorMat(2,1) = -(2*OMEGAedot + long_dot)*sin(xs_0(1));
% CorMat(2,2) = 0;
% CorMat(2,3) = -(2*OMEGAedot + long_dot)*cos(xs_0(1));
% CorMat(3,1) = lat_dot;
% CorMat(3,2) = (2*OMEGAedot + long_dot)*cos(xs_0(1));
% CorMat(3,3) = 0;
% 
% Coriolis = CorMat * xs_0(4:6);



Coriolis(1) = - 2*OMEGAedot*sin(lat)*Ve + ( Vn*Vd - Ve^2*tan(lat))/RMh;
Coriolis(2) = + 2*OMEGAedot*(Vn*sin(lat) + Vd*cos(lat)) + (Ve*Vd)/RPh + Ve*Vn*tan(lat)/RPh;
Coriolis(3) = - 2*OMEGAedot*Ve*cos(lat) - (Ve^2/RPh) - (Vn^2/RMh);  



A_xn = A_xn + Coriolis(1);
A_yn = A_yn + Coriolis(2);
A_zn = A_zn + Coriolis(3);

        

% RECTANGULAR INTEGRATION POSITION UPDATE       
% propogate velocity estimate
xs_0_k(4) = Vn + A_xn * ins_dt;
xs_0_k(5) = Ve + A_yn * ins_dt;
xs_0_k(6) = Vd + A_zn * ins_dt;


% determine latitude and longitude rates
lat_dot = xs_0_k(4) / (RMh);
long_dot = xs_0_k(5) / (RPh*cos(xs_0(1)));
hgt_dot = - xs_0_k(6);

% propogate position estimate   
xs_0_k(1) = xs_0(1) + lat_dot * ins_dt;
xs_0_k(2) = xs_0(2) + long_dot * ins_dt;
xs_0_k(3) = xs_0(3) + hgt_dot * ins_dt;  % note that x(6) is 
                                     % down velocity so 
                                     % we need to make 
                                     % negative before 
                                     % adding to height
                                     
                                     
x_out = xs_0_k';

output.C_BN = C_BN;
output.Euler(1) = phi_q;
output.Euler(2) = theta_q;
output.Euler(3) = psi_q;
output.llhdot(1) = lat_dot;
output.llhdot(2) = long_dot;
output.llhdot(3) = hgt_dot;
