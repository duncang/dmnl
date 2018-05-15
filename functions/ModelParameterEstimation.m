function [x_DOT, AeroCoeffs, Forces, Moments] = ModelParameterEstimation(M0,Gravity,Atmos,x0, delta_x_control,AeroCoeffsTruth, EngCoeffsTruth, PropCoeffsTruth,ForcesTruth, MomentsTruth, MassTruth,InertiaTruth,CGposTruth,MachTruth,Wind_NED); 



%This function calculates the parameters in the model in real-time using Total Least
%Squares




% CL = CL0 + CLa*alpha + (CLq*q*MAC)/(2*V) + CLM*M + CLde*de;  
% CD = CD0 + ((CL-CLD0)^2)*k + CDM*M + (CDde*de) + (CDda*da) + (CDdr*dr);  %why is this abs values? doesnt matter for navion anyway since the coeffs. are zero.
% %CD = CD0 + ((CL-CLD0)^2)*k + CDM*M + (CDde*de) + (CDda*da) + (CDdr*dr);
% CY = CYbeta*beta+ ((CYp*r + CYr*r)*b)/(2*V) + CYda*da + CYdr*dr;
% % Pitch, Roll and Yaw Moment Coefficients
% Cm = Cm0 + Cma*alpha + ((Cmq*q+ Cmalphadot*alpha_dot)*MAC)/(2*V) + CmM*M + Cmde*de;
% Cl = Clbeta*beta + ((Clp*p + Clr*r)*b)/(2*V) + Clda*da + Cldr*dr;
% Cn = Cnbeta*beta + ((Cnp*p + Cnr*r)*b)/(2*V) + Cnda*da + Cndr*dr;



% %drag 
% DRAG/(Q*S) = CD;
% %lift
% LIFT = Q*S*CL;
% 
% %sideforce
% SIDEFORCE = Q*S*CY;


DRAG_Measured = 
LIFT_Measured = 
SIDEFORCE_Measured = 


%drag 
DRAG_Measured/(Q*S) = CD0 + ((CL-CLD0)^2)*k + CDM*M + (CDde*de) + (CDda*da) + (CDdr*dr);
%lift
LIFT_Measured/(Q*S) = CL0 + CLa*alpha + (CLq*q*MAC)/(2*V) + CLM*M + CLde*de;  

%sideforce
SIDEFORCE_Measured/(Q*S) = CYbeta*beta+ ((CYp*r + CYr*r)*b)/(2*V) + CYda*da + CYdr*dr;


% L = Q*S*b*Cl;
% 
% M = Q*S*c_bar*Cm;
% 
% N = Q*S*b*Cn;


L_Measured/(Q*S*b) = Clbeta*beta + ((Clp*p + Clr*r)*b)/(2*V) + Clda*da + Cldr*dr;

M_Measured/(Q*S*c_bar) = Cm0 + Cma*alpha + ((Cmq*q+ Cmalphadot*alpha_dot)*MAC)/(2*V) + CmM*M + Cmde*de;

N_Measured/(Q*S*b) = Cnbeta*beta + ((Cnp*p + Cnr*r)*b)/(2*V) + Cnda*da + Cndr*dr;




