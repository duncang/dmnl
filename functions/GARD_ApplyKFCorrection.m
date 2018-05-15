function [x_out,acc_bias_out,gyro_bias_out,clock_out] = GARD_ApplyKFCorrection(x_in,acc_bias_in,gyro_bias_in,clock_in,x_hat)

%
% function x_out = GARD_ApplyKFCorrection(x_in,x_hat)
% 
% This function applies the state correction, x_hat, to the
% full state vector x_in.
%
% x1  - Latitude (rad)
% x2  - Longitude (rad)
% x3  - Height (m)
% x4  - Vel North (m/s)
% x5  - Vel East (m/s)
% x6  - Vel Down (m/s)
% x7  - q0
% x8  - q1
% x9  - q2
% x10 - q3
% x11 - 
% 
% correction vector, x_hat
% x1:3 - LLH Error
% x4:6 - Vel Error
% x7:9 - Tilt Error
% x10:12 - Accel Bias
% x13:15 - Gyro Bias
% x16:17 - Clock Error
% 
% $Id$
%

%% apply correction
x_out(4:6) = x_in(4:6) + x_hat(4:6);
x_out(1:3) = x_in(1:3) + x_hat(1:3);

clock_out(1) =  clock_in(1) + x_hat(16);
clock_out(2) =  clock_in(2) + x_hat(17);


acc_bias_out = acc_bias_in - x_hat(10:12)';   
gyro_bias_out = gyro_bias_in - x_hat(13:15)';  

% convert tilt error to DCM Update
del_alpha = x_hat(7);
del_beta = x_hat(8);
del_gamma = x_hat(9);

% correct atttitude
del_att_skew = [0         -del_gamma   del_beta; ...
               del_gamma  0          -del_alpha; ...
               -del_beta  del_alpha   0];


C_BN = GARD_QuatToDCM(x_in(7:10));
           
C_BN = (eye(3,3) - del_att_skew) * C_BN; % DCM correction

C_BN = GARD_OrthogonaliseDCM(C_BN);


x_out(7:10) = GARD_DCMToQuat(C_BN);

