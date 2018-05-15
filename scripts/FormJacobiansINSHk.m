%note you need to 


%forms the Jacobians for the INS equations

%declare state variables
% syms delta_x phi theta psi vn ve vd lat lon hgt vn_dot ve_dot vd_dot phi_dot theta_dot psi_dot
% 
%  %sensor measurements
% syms  omega_x omega_y omega_z A_b A_xb A_yb A_zb 
% 
% syms  Rmeridian Rnormal
% 
% syms  C_bn lon_dot lat_dot
% 
% syms  lat_ini lon_ini hgt_ini vn_ini ve_ini vd_ini  %these are the initial values
% 
% syms CorMat g

syms a b e lat lon hgt phi theta psi Xs Ys Zs Xu Yu Zu N ele1 ele2 ele3 r_VecCalc vn ve vd Xs_vel Ys_vel Zs_vel R0

%previous state


y = cos(a+b);


jacobian(y,b)


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
delta_x = [phi, theta, psi, vn, ve, vd, lat, lon, hgt];




% WGS-84 ellipsoid parameters
% a = 6378137.0; % semi-major axis
% b = 6356752.3142; % semi-minor axis

% calculate the first eccentricity
%e = sqrt(a^2 - b^2) / a;

%p = 
%N = p / cos(Latitude);

% N = a^2 / sqrt(a^2 * cos(lat)^2 + b^2 * sin(lat)^2);
% 
% 
% 
% % 
% Xu = (N + hgt) * cos(lat) * cos (lon);
% Yu = (N + hgt) * cos(lat) * sin (lon);
% Zu = ((1 - e^2)*N + hgt) * sin(lat);




% 
Re = R0*(1+e*sin(lat)^2);

Xu = (Re + hgt) * cos(lat) * cos (lon);
Yu = (Re + hgt) * cos(lat) * sin (lon);
Zu = (Re*(1-2*e) + hgt) * sin(lat);


%Xu, Yu, Zu is user position in ECEF


% 
%   for k = 1:NumberMeasurements
%         for m = 1:3
%             ele(m) =  SatPos(i,k,m) - UserPos(m);
%         end
%         r_VecCalc(k) =  norm(ele);
        
        
        
        ele1 = Xs - Xu;
        ele2 = Ys - Yu;
        ele3 = Zs - Zu;


        r_VecCalc = sqrt(ele1^2 + ele2^2 +ele3^2);
        
        

f(1) = r_VecCalc;


%form Jacobian for position in delta_x


Hk = jacobian(f,delta_x)

%this should form part of the Hk matrix for one satellite, just repeat it
%for the others.


%form jacobian for velocity in NED

TMatrix = T_ECEF2NED(lon,lat);

%NED to ECEF velocities

VelECEF = TMatrix'*[vn ve vd]';
Xu_vel = VelECEF(1);
Yu_vel = VelECEF(2);
Zu_vel = VelECEF(3);

% z(k,1) = [Pseudoranges(i,k) - r_VecCalc(i,k)];
        %predicted relative velocity of sv and receiver
        r_VecCalcVel = (Xs_vel - Xu_vel)*(Xs-Xu) + (Ys_vel - Yu_vel)*(Ys-Yu) + (Zs_vel - Zu_vel)*(Zs-Zu);

        %   r_VecCalcVel(k) = (SVVel(k,1) - TempBase_B_Vel(1))*(SVPos(k,1)-TempBase_B_Pos(1)) + (SVVel(k,2) - TempBase_B_Vel(2))*(SVPos(k,2)-TempBase_B_Pos(2)) + (SVVel(k,3) - TempBase_B_Vel(3))*(SVPos(k,3)-TempBase_B_Pos(3));

      
        Relative_Velocity = r_VecCalcVel/r_VecCalc;



g(1) = Relative_Velocity;
Hkned = jacobian(g,delta_x)




% 
% 
% 
% 
% 
% % omega_b = [p_b, q_b r_b]';  %roll pitch yaw rates around body axes 
% % 
% % P_n = C_bn*P_b   %convert to roll pitch and yaw around navigation axes
% % 
% % 
% % p_n = P_n(1);
% % q_n = P_n(2);
% % r_n = P_n(3); 
% 
% 
% %phi dot
% % f(1) = omega_x +tan(theta)*(omega_y*sin(phi) + omega_z*cos(phi));
% % 
% % %theta dot
% % f(2) = omega_y*cos(phi) - omega_z*sin(phi);
% % 
% % %psi dot
% % f(3) = (omega_y*sin(phi) + omega_z*cos(phi))/cos(theta);
% 
% 
%  f(1) = phi_dot;
% % 
% % %theta dot
%  f(2) = theta_dot;
% % 
% % %psi dot
%  f(3) = psi_dot;
% 
% 
% V_n = [vn, ve, vd]';
% 
% A_b = [A_xb, A_yb, A_zb]'; %accelerations in body
% 
% % C_bn = GARDSim_DCMfromEuler(phi,theta,psi);
% 
% 
% 
%   % perform coriolis correction
%        OMEGA_e = 7.292115e-5; % eearth rotation rate (inertial frame)
%         
%         
%         CorMat(1,1) = 0;
%         CorMat(1,2) = (2*2 + lon_dot*sin(lat));
%         CorMat(1,3) = -lat_dot;
%         CorMat(2,1) = -(2*OMEGA_e + lon_dot*sin(lat));
%         CorMat(2,2) = 0;
%         CorMat(2,3) = -(2*OMEGA_e + lon_dot*cos(lat));
%         CorMat(3,1) = lat_dot;
%         CorMat(3,2) = (2*OMEGA_e + lon_dot*cos(lat));
%         CorMat(3,3) = 0;
%         
%         Coriolis = CorMat * V_n;
%  
%         
%         
%         
% 
%        % g is local gravity vector based on WGS-84 earth gravitation model
%        
% v_ned_dot = C_bn*A_b + Coriolis + [0, 0, -g]'; 
% 
% 
% % %vn dot
% % f(4) = v_ned(1);
% % 
% % %ve dot
% % f(5) = v_ned(2);
% % 
% % %vd dot
% % f(6) =  v_ned(3);
% 
%  f(4) = v_ned_dot(1);
% % 
% % %ve dot
% f(5) = v_ned_dot(2);
% % 
% % %vd dot
%  f(6) =  v_ned_dot(3);
% 
% 
% 
% 
% %convert to lat lon hgt dot
% 
% 
% %find the meridian radius of curvature- this needs to verified.
% 
% % a = 6378137.0;   % semi-major axis (metres)
% % f = 1/298.2572; % flattening
% % e2 = f * (2-f); % eccentricity squared
% % e = sqrt(e2);   % first eccentricity
% 
% % Rmeridian = a * (1 - e2) / (sqrt(1 - e2 * sin(Lat)^2))^3;
% % %find the normal radius of curvature
% % Rnormal = a / sqrt(1 - e2 * sin(Lat)^2);
% 
% 
% 
% 
% %Lat_dot
% f(7) = vn/(Rmeridian + hgt);
% 
% %Lon_dot
% f(8) = ve/((Rnormal + hgt)*cos(lat));
% 
% %hgt_dot
% f(9) = vd;  %technically this should be -h_dot (according to aerosim) but wiht it like this it works for some reason. 
% 
% 
% 
% 
% %Fk will be the phi matrix in the KF
% Fk = jacobian(f,delta_x)
% 
% 
% 
% 
% 
% 
% 
% 







