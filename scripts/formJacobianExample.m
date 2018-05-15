%note you need to 

%declare state variables
syms u v w phi 

%declare other symbols representing constants (like gravity, mass etc)
syms theta psi p q r lat lon hgt g alpha beta Fxa 



%now these are the equations of motion

%gp = u^2 + v^2;
%u dot
%f(1) = r*v - q*w - g*sin(theta) + (1/m)*(1/((w/u)^2 + 1)^0.5 )*((1 - (v^2/(u^2 + v^2 + w^2)))^0.5);






%Fxa = -0.5*rho*V_T^2*S*CD;
Fxa = u^2;
f(1) = r*v - q*w - g*sin(theta) + (1/m)*cos(alpha)*cos(beta)*Fxa



%f(1) = r*v - q*w - g*sin(theta) + (1/m)*cos(alpha)*cos(beta)*Fxa





Jac = jacobian(f,[u v])







