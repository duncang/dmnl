function F = GARD_GenerateINSFMatrix(Pos_LLH,Vel_NED,Acc_NED)
% function F = GARD_GenerateINSFMatrix(Pos_LLH,Vel_NED,Acc_NED)
% Written by Duncan Greer (c) 9 November 2007
%
% $Id: GARD_GenerateINSFMatrix.m 4047 2010-10-19 06:19:05Z greerd $
% 
% This is a direct representation of the error model presented in Appendix
% A of Eck "Error Dynamics of MOdel Based INS/GPS Navigation for an
% Autonmously Flying Helicopter" AIAA 2000.
%

% load constants
GPSConstants;
LOCAL_GRAVITY = -GravityModel(Pos_LLH); 

RM = MeridianRadius(Pos_LLH(1));  % north-south
RP = PrimeRadius(Pos_LLH(1));     % east-west

R = sqrt(RM^2 + RP^2);

Rh = R + Pos_LLH(3);
RMh = RM + Pos_LLH(3);
RPh = RP + Pos_LLH(3);

lat = Pos_LLH(1);
lon = Pos_LLH(2);

%% F matrix from Titterton
F = zeros(9,9);

%% position error due to position error
F(1,3) = -Vel_NED(1)/RMh^2;
F(2,1) = Vel_NED(2)*tan(lat)/(RPh*cos(lat));
F(2,3) = -Vel_NED(2)/(RPh*RPh*cos(lat));

%% position error due to velocity error
F(1,4) = 1/RMh;
F(2,5) = 1/(RPh*cos(lat));
F(3,6) = -1;


%% velocity error due to position error
F(4,1) = -Vel_NED(2)*2*OMEGAedot*cos(lat) - Vel_NED(2)*Vel_NED(2)/(RPh*cos(lat)^2);
F(4,3) = Vel_NED(2)^2 * tan(lat)/RPh^2 - Vel_NED(1)*Vel_NED(3)/RMh^2;
F(5,1) = 2*OMEGAedot*(Vel_NED(1)*cos(lat) - Vel_NED(3)*sin(lat)) + (Vel_NED(1)*Vel_NED(2)/(RPh*cos(lat)^2));
F(5,3) = -Vel_NED(2)*(Vel_NED(1)*tan(lat) + Vel_NED(3))/RPh^2;
F(6,1) = 2*OMEGAedot*Vel_NED(2)*sin(lat);
F(6,3) = (Vel_NED(1)^2)/RMh^2 + (Vel_NED(2)^2)/RPh^2 - 2.0*LOCAL_GRAVITY/(sqrt(RM*RP));

%% velocity error due to velocity error
F(4,4) = Vel_NED(3)/RMh;
F(4,5) = -2*(OMEGAedot*sin(lat) + Vel_NED(2)*tan(lat)/RPh);
F(4,6) = Vel_NED(1)/RMh;
F(5,4) = 2*OMEGAedot*sin(lat) + Vel_NED(2)*tan(lat)/RPh;
F(5,5) = (Vel_NED(1)*tan(lat)+Vel_NED(3))/RPh;
F(5,6) = 2*OMEGAedot * cos(lat) + Vel_NED(2)/RPh;
F(6,4) = -2*Vel_NED(1)/RMh;
F(6,5) = -2*(OMEGAedot * cos(lat) + Vel_NED(2)/RPh);


%% velocity error due to tilt error
F(4,8) = -(Acc_NED(3));  %-fd
F(4,9) = (Acc_NED(2)); % fe
F(5,7) = (Acc_NED(3));   %fd
F(5,9) = -(Acc_NED(1)); % -fn
F(6,7) = -(Acc_NED(2)); % -fe
F(6,8) = (Acc_NED(1)); % fn

 
%% tilt error due to position error
F(7,1) = -OMEGAedot * sin(lat);
F(7,3) = -Vel_NED(2)/RPh^2;
F(8,3) = Vel_NED(1)/RMh^2;
F(9,1) = -OMEGAedot * cos(lat) - Vel_NED(2)/(RPh*cos(lat)^2);
F(9,3) = Vel_NED(2)*tan(lat)/RPh^2;

%% tilt error due to velocity error
F(7,5) = 1/RPh;
F(8,4) = -1/RMh;           
F(9,5) = -tan(lat)/RPh;



%% tilt error due to inertial rotation of the reference frame
F(7,8) = -OMEGAedot*sin(lat) - Vel_NED(2)*tan(lat)/RPh;
F(7,9) = Vel_NED(1)/RMh;
F(8,7) = OMEGAedot*sin(lat) + Vel_NED(2)*tan(lat)/RPh;
F(8,9) = OMEGAedot*cos(lat) + Vel_NED(2)/RPh;
F(9,7) = -Vel_NED(1)/RMh;
F(9,8) = -OMEGAedot*cos(lat) - Vel_NED(2)/RPh;


