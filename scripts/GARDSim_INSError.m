%% INS FULL STATE ERROR ANALYSIS
% From J.A. Farrell & M. Barth "The Global Positioning System & Inertial
% Navigation" 1999 McGraw-Hill - Chapter 6


% State Vector
NumberEpochs = 3600*36;  % 36 hours

x = zeros(9,NumberEpochs);

latitude = 45.0;
longitude = 0.0;
h = 0;

g = 9.81;

f_D = -g;
f_N = 0;
f_E = 0;

v_N = 0;
v_E = 0;
v_D = 0;



%% position
lat = latitude*pi/180;
long = longitude*pi/180;

% earth rate calculations
omega_ie = 7.2921151467e-5;
a = 6378137.0;
b = 6356752.3;
f = (a-b) / a;
e = sqrt(f*(2-f));
R_lamda = a*(1-e^2) / (1 - e^2 * sin(lat)^2)^1.5;
R_phi = a / sqrt(1 - e^2 * sin(lat)^2);

OMEGA = omega_ie;
R = sqrt(R_lamda^2 + R_phi^2);

%  
%  OMEGA_N = omega_ie * cos(lamda);
%  OMEGA_D = -omega_ie * sin(lamda);
%  
%  omega_N = OMEGA_N + v_E / (R_phi + h);
%  omega_E = -v_N / (R_lamda + h);
%  omega_D = OMEGA_D + -v_E * tan(lamda) / (R_phi + h);
%  
%  F_p_p = zeros(3,3);
%  F_p_p(1,3) = -v_N / (R_lamda + h)^2;
%  F_p_p(2,1) = v_E * sin(lamda) / ((R_phi + h) * cos(lamda)^2);
%  F_p_p(2,3) = -v_E / ((R_phi + h)^2 * cos(lamda));
%  
%  F_p_v = zeros(3,3);
%  F_p_v(1,1) = 1/(R_lamda+h);
%  F_p_v(2,2) = 1/((R_phi + h)*cos(lamda));
%  F_p_v(3,3) = -1;
%  
%  
%  F_p_rho = zeros(3,3);
%  
%  F_v_p = zeros(3,3);
%  F_v_p(1,1) = -2*OMEGA_N*v_E - v_E^2 / ((R_phi + h)*cos(lamda)^2);
%  F_v_p(1,3) = -v_N*v_D / (R_lamda + h)^2 + v_E^2 * tan(lamda) / (R_phi+h)^2;
%  F_v_p(2,1) = 2*(OMEGA_N * v_N + OMEGA_D * v_D) + v_E * v_N / (R_phi+h)*cos(lamda)^2;
%  F_v_p(2,3) = -v_E * (v_N * tan(lamda) + v_D)/(R_phi+h)^2;
%  F_v_p(3,1) = -v_E * 2 * OMEGA_D;
%  F_v_p(3,3) = v_E^2 / (R_phi+h)^2 + v_N^2/(R_lamda + h)^2 - 2*g/a;
%  
%  F_v_v = zeros(3,3);
%  F_v_v(1,1) = v_D / (R_lamda+h);
%  F_v_v(1,2) = -2*v_E*sin(lamda) / ((R_phi+h)*cos(lamda)) + 2*OMEGA_D;
%  F_v_v(1,3) = v_N / (R_lamda+h);
%  F_v_v(2,1) = v_E * tan(lamda) / (R_phi+h) - 2 * OMEGA_D;
%  F_v_v(2,2) = v_N * tan(lamda) / (R_phi+h) + v_D/(R_phi+h);
%  F_v_v(2,3) = v_E/(R_phi+h) + 2 * OMEGA_N;
%  F_v_v(3,1) = -2 * v_N / (R_lamda+h);
%  F_v_v(3,2) = -2 * v_E / (R_phi+h) - 2 * OMEGA_N;
%  
%  F_v_rho = zeros(3,3);
%  F_v_rho(1,2) = f_D;
%  F_v_rho(1,3) = -f_E;
%  F_v_rho(2,1) = -f_D;
%  F_v_rho(2,3) = f_N;
%  F_v_rho(3,1) = f_E;
%  F_v_rho(3,2) = -f_N;
%  
%  F_rho_p = zeros(3,3);
%  F_rho_p(1,1) = omega_ie * sin(lamda);
%  F_rho_p(1,3) = v_E / ((R_phi + h)^2);
%  F_rho_p(2,3) = -v_N / ((R_lamda + h)^2);
%  F_rho_p(3,1) = omega_ie * cos(lamda) + v_E/((R_phi + h)*cos(lamda)^2);
%  F_rho_p(3,3) = -v_E * tan(lamda) / ((R_phi + h)^2);
%  
%  F_rho_v = zeros(3,3);
%  F_rho_v(1,2) = -1 / (R_phi + h);
%  F_rho_v(2,1) = 1/(R_lamda + h);
%  F_rho_v(3,2) = tan(lamda)/(R_phi + h);
%  
%  F_rho_rho = zeros(3,3);
%  F_rho_rho(1,2) = omega_D;
%  F_rho_rho(1,3) = omega_E;
%  F_rho_rho(2,1) = -omega_D;
%  F_rho_rho(2,3) = omega_N;
%  F_rho_rho(3,1) = omega_E;
%  F_rho_rho(3,2) = -omega_N;
%  
%  
%  
%  
%  F_INS = [F_p_p F_p_v F_p_rho;
%           F_v_p F_v_v F_v_rho;
%           F_rho_p F_rho_v F_rho_rho];
%       


F_INS = zeros(9,9);

F_INS(1,2) = -(OMEGA*sin(lat) + v_E/R * tan(lat));
F_INS(1,3) = v_N / R;
F_INS(1,5) = 1/R;
F_INS(1,7) = -OMEGA*sin(lat);
F_INS(1,9) = -v_E/R^2;
F_INS(2,1) = (OMEGA*sin(lat) + v_E/R * tan(lat));
F_INS(2,3) = OMEGA*cos(lat) + v_E/R;
F_INS(2,4) = -1/R;
F_INS(2,9) = v_N/R^2;
F_INS(3,1) = -v_N/R;
F_INS(3,2) = -OMEGA*cos(lat) - v_E/R;
F_INS(3,6) = -tan(lat) / R;
F_INS(3,7) = -OMEGA*cos(lat) - v_E/(R*cos(lat)^2);
F_INS(3,9) = v_E*tan(lat)/R^2;
F_INS(4,2) = -f_D;
F_INS(4,3) = f_E;
F_INS(4,4) = v_D/R;
F_INS(4,5) = -2*(OMEGA*sin(lat) + v_E/R * tan(lat));
F_INS(4,6) = v_N / R;
F_INS(4,7) = -v_E*(2*OMEGA*cos(lat) + v_E/(R*cos(lat)^2));
F_INS(4,9) = 1/R^2 * (v_E^2 * tan(lat) - v_N * v_D);
F_INS(5,1) = f_D;
F_INS(5,3) = -f_N;
F_INS(5,4) = 2*OMEGA*sin(lat) + v_E/R * tan(lat);
F_INS(5,5) = 1/R * (v_N*tan(lat) + v_D);
F_INS(5,6) = 2*OMEGA*cos(lat) + v_E/R;
F_INS(5,7) = 2*OMEGA*(v_N*cos(lat)-v_D*sin(lat)) + v_N*v_E/(R*cos(lat)^2);
F_INS(5,9) = -v_E/R^2 * (v_N * tan(lat) + v_D);
F_INS(6,1) = -f_E;
F_INS(6,2) = f_N;
F_INS(6,4) = -2*v_N / R;
F_INS(6,5) = -2*(OMEGA*cos(lat)+v_E/R);
F_INS(6,7) = 2*OMEGA*v_E*sin(lat);
F_INS(6,9) = 1/R^2 * (v_N^2 + v_E^2);
F_INS(7,4) = 1/R;
F_INS(7,9) = -v_N/R^2;
F_INS(8,5) = 1/(R*cos(lat));
F_INS(8,7) = v_E * tan(lat) / (R*cos(lat));
F_INS(8,9) = -v_E / (R^2 * cos(lat));
F_INS(9,6) = -1;
         
dt = 1;

NumberStates = size(F_INS);

PHI = eye(NumberStates) + F_INS * dt + 0.5*F_INS*dt^2;

x(:,1) = zeros(9,1);
x(1:3,1) = [0.1e-3,0.1e-3,1e-3]';


for Epoch=1:NumberEpochs-1
    x(:,Epoch+1) = PHI*x(:,Epoch);
end
    

