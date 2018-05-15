%note you need to 


%forms the Jacobians for the INS equations

%declare state variables
syms delta_x phi theta psi vn ve vd lat lon hgt vn_dot ve_dot vd_dot phi_dot theta_dot psi_dot

 %sensor measurements
syms  omega_x omega_y omega_z A_b A_xb A_yb A_zb 

syms  Rmeridian Rnormal

syms  C_bn lon_dot lat_dot

syms  lat_ini lon_ini hgt_ini vn_ini ve_ini vd_ini  %these are the initial values

syms CorMat g

syms a b e Latitude Longitude Height Xs Ys Zs Xu Yu Zu N ele1 ele2 ele3 r_VecCalc p_b q_b r_b


syms  vne vee vde X Y Z

%previous state





%declare state parameter estimates for augmented state vector



% %declare other symbols representing constants (like gravity, mass etc)
% syms g_dash0 alpha beta M V_T
% 
% %
% syms u v w phi theta psi p q r Lat Lon Hgt


%differential equation for velocity and position





% WGS-84 ellipsoid parameters
% a = 6378137.0; % semi-major axis
% b = 6356752.3142; % semi-minor axis






%this is the state 
%delta_x = [phi, theta, psi, vne, vee, vde, X, Y, Z];  %velocities and position are in ECEF

%convert ECEF position and velocities into lat lon hgt and NED velocities 












% X = Position(1);
% Y = Position(2);
% Z = Position(3);

% % WGS-84 ellipsoid parameters
% a = 6378137.0; % semi-major axis
% b = 6356752.3142; % semi-minor axis


% P = sqrt(X^2 + Y^2);
% 
% Theta = atan((Z * a) / (P * b));
% 
% Esq = 1 - b^2 / a^2;
% EPsq = a^2 / b^2 - 1;
% 
% lat = atan(((Z + EPsq * b * sin(Theta)^3)) / (P- Esq * a * cos(Theta)^3));
% lon = atan(Y,X);
% n = a^2 / sqrt(a^2 * cos(lat)^2 + b^2 * sin(lat)^2);
% hgt = P / cos(lat) - n;
% 
% 
% 
  TMatrix = T_ECEF2NED(lon,lat);
% 
%  
  Vned = TMatrix*[vne, vee, vde]';
 
 
 


%this is the state 
delta_x = [phi, theta, psi, vn, ve, vd, lat, lon, hgt];



P_b = [p_b, q_b r_b]';  %roll pitch yaw rates around body axes 

P_n = C_bn*P_b   %convert to roll pitch and yaw around navigation axes


p_n = P_n(1);
q_n = P_n(2);
r_n = P_n(3); 


%phi dot
f(1) = omega_x +tan(theta)*(omega_y*sin(phi) + omega_z*cos(phi));

%theta dot
f(2) = omega_y*cos(phi) - omega_z*sin(phi);

%psi dot
f(3) = (omega_y*sin(phi) + omega_z*cos(phi))/cos(theta);


%  f(1) = phi_dot;
% % 
% % %theta dot
%  f(2) = theta_dot;
% % 
% % %psi dot
%  f(3) = psi_dot;


 
 vn = Vned(1);
 ve = Vned(2);
 vd = Vned(3);
 
 
V_n = [vn, ve, vd]';

A_b = [A_xb, A_yb, A_zb]'; %accelerations in body

% C_bn = GARDSim_DCMfromEuler(phi,theta,psi);



  % perform coriolis correction
       OMEGA_e = 7.292115e-5; % eearth rotation rate (inertial frame)
        
        
        CorMat(1,1) = 0;
        CorMat(1,2) = (2*OMEGA_e + lon_dot*sin(lat));
        CorMat(1,3) = -lat_dot;
        CorMat(2,1) = -(2*OMEGA_e + lon_dot*sin(lat));
        CorMat(2,2) = 0;
        CorMat(2,3) = -(2*OMEGA_e + lon_dot*cos(lat));
        CorMat(3,1) = lat_dot;
        CorMat(3,2) = (2*OMEGA_e + lon_dot*cos(lat));
        CorMat(3,3) = 0;
        
        Coriolis = CorMat * V_n;
 
        
        
        

       % g is local gravity vector based on WGS-84 earth gravitation model
       
v_ned_dot = C_bn*A_b + Coriolis + [0, 0, -g]'; 


% %vn dot
% f(4) = v_ned(1);
% 
% %ve dot
% f(5) = v_ned(2);
% 
% %vd dot
% f(6) =  v_ned(3);

 f(4) = v_ned_dot(1);
% 
% %ve dot
f(5) = v_ned_dot(2);
% 
% %vd dot
 f(6) =  v_ned_dot(3);




%convert to lat lon hgt dot


%find the meridian radius of curvature- this needs to verified.

% a = 6378137.0;   % semi-major axis (metres)
% f = 1/298.2572; % flattening
% e2 = f * (2-f); % eccentricity squared
% e = sqrt(e2);   % first eccentricity

% Rmeridian = a * (1 - e2) / (sqrt(1 - e2 * sin(Lat)^2))^3;
% %find the normal radius of curvature
% Rnormal = a / sqrt(1 - e2 * sin(Lat)^2);




%Lat_dot
f(7) = vn/(Rmeridian + hgt);

%Lon_dot
f(8) = ve/((Rnormal + hgt)*cos(lat));

%hgt_dot
f(9) = -vd;  %technically this should be -h_dot (according to aerosim) but wiht it like this it works for some reason. 




%Fk will be the phi matrix in the KF
Fk = jacobian(f,delta_x)

%Phik = expm(Fk);













