%note you need to 

%declare state variables
syms u v w phi theta psi p q r Lat Lon Hgt

%declare state parameter estimates for augmented state vector




%declare other symbols representing constants (like gravity, mass etc)
syms g_dash0 alpha beta M V_T

%mass and moments of inertia

syms m Jx Jy Jz Jxz

%atmospheric parameters

syms rho a 

%declare parameters

syms CL0 CLa CLq CLM CLde           %for CL
syms CD0 CL CLD0 k CDM CDde CDda CDdr  %for CD
syms CYbeta CYp CYr CYda CYdr           %for CY


syms Cm0 Cma Cmq Cmalphadot CmM Cmde  %for Cm
syms Clbeta Clp Clr Clda Cldr       %for Cl

syms Cnbeta Cnp Cnr Cnda Cndr           %for Cn

% CL = CL0 + CLa*alpha + (CLq*q*MAC)/(2*V_T) + CLM*M + CLde*de;  
% CD = CD0 + ((CL-CLD0)^2)*k + CDM*M + (CDde*de) + (CDda*da) + (CDdr*dr);  
% CY = CYbeta*beta+ ((CYp*r + CYr*r)*b)/(2*V_T) + CYda*da + CYdr*dr;


% Cm0 + Cma + ((Cmq+ Cmalphadot*alpha_dot)*MAC)/(2*V) + CmM*M + Cmde*de;
% Clbeta + ((Clp + Clr))) + Clda*da + Cldr*dr;
% Cnbeta + ((Cnp + Cnr)) + Cnda*da + Cndr*dr;

 
%declare control inputs
syms de da dr


syms Rmeridian Rnormal

syms c_bar b S MAC

syms alpha_dot



Q = 0.5*rho*V_T^2;


%  alpha = atan(w/u);
%  beta = asin(v/V_T);




CL = CL0 + CLa*alpha + (CLq*q*MAC)/(2*V_T) + CLM*M + CLde*de;  
CD = CD0 + ((CL-CLD0)^2)*k + CDM*M + (CDde*de) + (CDda*da) + (CDdr*dr);  
CY = CYbeta*beta+ ((CYp*r + CYr*r)*b)/(2*V_T) + CYda*da + CYdr*dr;



Fxa = -Q*S*CD; %Fxa = -D
Fya = Q*S*CY; %Fya = Y
Fza = -Q*S*CL; %Fza = -L


TMatrixBody2Wind = T_Body2Wind(alpha, beta);
Fbody = TMatrixBody2Wind'*[Fxa;Fya;Fza];

Fxb = Fbody(1);
Fyb = Fbody(2);
Fzb = Fbody(3);



% Calpha = cos(alpha);
% Cbeta = cos(beta)
% Salpha = sin(alpha);
% Sbeta = sin(beta);


%Fxb = Calpha*Cbeta*Fxa - Calpha*Sbeta*Fya - Salpha*Fza;



%Fxb = 

%u dot




f(1) = r*v - q*w - g_dash0*sin(theta) + Fxb/m;


%v dot


f(2) = -r*u + p*w + g_dash0*sin(phi)*cos(theta) + Fyb/m;


%w dot

f(3) = q*u - p*v + g_dash0*cos(phi)*cos(theta) + Fzb/m;





%form L M and N

Cm = Cm0 + Cma*alpha + ((Cmq*q+ Cmalphadot*alpha_dot)*MAC)/(2*V_T) + CmM*M + Cmde*de;
Cl = Clbeta*beta + ((Clp*p + Clr*r)*b)/(2*V_T) + Clda*da + Cldr*dr;
Cn = Cnbeta*beta + ((Cnp*p + Cnr*r)*b)/(2*V_T) + Cnda*da + Cndr*dr;


L = Q*S*b*Cl;

M = Q*S*c_bar*Cm;

N = Q*S*b*Cn;

MaeroMyModel = [L, M, N];








GAMMA_J = Jx*Jz - Jxz^2;
c1 = ((Jy-Jz)*Jz - Jxz^2)/GAMMA_J;
c2 = ((Jx-Jy+Jz)*Jxz)/GAMMA_J;
c3 = Jz/GAMMA_J;
c4 = Jxz/GAMMA_J;
c5 = (Jz-Jx)/Jy;
c6 = Jxz/Jy;
c7 = 1/Jy;
c8 = (Jx*(Jx-Jy) + Jxz^2)/GAMMA_J;
c9 = Jx/GAMMA_J;



%p q and r dot


f(4) = (c1*r + c2*p)*q + c3*L +c4*N;
f(5) = c5*p*r - c6*(p^2 - r^2) + c7*M;
f(6) = (c8*p - c2*r)*q + c4*L + c9*N;




%phi theta psi dot
f(7) = p +tan(theta)*(q*sin(phi) + r*cos(phi));
f(8) = q*cos(phi) - r*sin(phi);
f(9) = (q*sin(phi) + r*cos(phi))/cos(theta);



Cphi = cos(phi); 
Ctheta = cos(theta);
Cpsi = cos(psi);

Sphi = sin(phi);
Stheta = sin(theta);
Spsi = sin(psi);


pn_dot = u*Ctheta*Cpsi + v*(-Cphi*Spsi + Sphi*Stheta*Cpsi) + w*(Sphi*Spsi+Cphi*Stheta*Cpsi);
pe_dot = u*Ctheta*Spsi + v*(Cphi*Cpsi + Sphi*Stheta*Spsi) + w*(-Sphi*Cpsi + Cphi*Stheta*Spsi);
h_dot = u*Stheta - v*Sphi*Ctheta - w*Cphi*Ctheta; 




% a = 6378137.0;   % semi-major axis (metres)
% f = 1/298.2572; % flattening
% e2 = f * (2-f); % eccentricity squared
% e = sqrt(e2);   % first eccentricity
% 
% 
% Rmeridian = a * (1 - e2) / (sqrt(1 - e2 * sin(Lat)^2))^3;
% 
% %find the normal radius of curvature
% Rnormal = a / sqrt(1 - e2 * sin(Lat)^2);

%Lat, Lon and Hgt dot


f(10) = pn_dot/(Rmeridian + Hgt);
f(11) = pe_dot/((Rnormal + Hgt)*cos(Lat));
f(12) = h_dot;




%f = [ f(1), f(2),f(3),f(4),f(5),f(6),f(7),f(8),f(9),f(10),f(11),f(12)]';


%compute jacobian wrt the state

state = [u v w phi theta psi p q r Lat Lon Hgt];


%note the first parameter f in the jacobian function must be a column
%vector and the second 'state' must be a row vector. 

%Fk will be the phi matrix in the KF
Fk = jacobian(f,state)


%compute jacobian wrt the control inputs

control = [de da dr];  %note this 'control has to be a row vector not a column vector for use with symbolic toolbox


%note the first parameter f in the jacobian function must be a column
%vector and the second 'state' must be a row vector. 

Gk = jacobian(f,control)




%evaluate the Hk matrix



state = [u v w phi theta psi p q r Lat Lon Hgt];

%lat
h(1) = int(f(10));  %int takes the integral of f(10)

%lon
h(2) = int(f(11));

%hgt
h(3) = int(f(12));

%lat dot
h(4) = f(10);
h(5) = f(11);
h(6) = f(12);

%phi
h(7) = int(f(7));

%theta
h(8) = int(f(8));

%psi
h(9) = int(f(9));


%phi dot
h(10) = f(7);

%theta dot
h(11) = f(8);

%psi dot
h(12) = f(9);

%p
h(13) = int(f(4));

%q
h(14) =  int(f(5));

%r
h(15) =  int(f(6));



%form Hk matrix, which is the jacobian of the nonlinear h 
 Hk = jacobian(h,state)



% %evaluate whether its acceptable to ignore the fact that alpha = atan(w/u)
% %in the u v and w jacobians, because there are too many terms. stick some
% %numbers in for alpha and beta and see if you get a similar result as
% %without them
% 
% 
% 
% 
% %the Fk and Gk matrices are :
% 
% 
% %to find the Gm matrix use Jac = jacobian(f,control) where control contains
% %the control variables de da dr


