function [x_hat_out, P_out, v_out, s2_out] = GARD_EvaluateKF(dt, x_hat_in, P_in, phi, H, z, Q, R)
% function [x_hat_out, P_out, v_out, s2_out] = GARD_EvaluateKF(dt,
% x_hat_in, P_in, phi, H, z, Q, R)
%
% dt = arguments(1);
% x_hat_in = arguments(2);
% P_in = arguments(3);
% phi = arguments(4);
% H = arguments(5);
% z = arguments(6);
% Q = arguments(7);
% R = arguments(8);
%
% [x_hat_out, P_out] = GARD_EvaluateKF(dt, x_hat_in, P_in, phi, H, z, Q, R)
% evaluates the 8 State GPS position and velocity kalman filter solution
% written by Peter Roberts
%
% adapted by duncan greer 27 September 2005
% Last Update: $Id: GARD_EvaluateKF.m 1850 2008-07-14 04:52:47Z greerd $
%
% Inputs
%   dt          - time step
%   x_hat_in    - previous (posteriori) state estimate from previous epoch
%   P_in        - previous (posteriori) estimate covariance
%   phi         - state transition matrix
%   H           - observation matrix
%   z           - measurements matrix
%   R           - sensor noise covariance matrix
%   Q           - Process noise covariance matrix
%   
%
% Outputs
%   x_hat_out - posteriori state estimate
%   P_out     - posteriori estimate covariance 
%   v_out     - innovation vector
%   s2_out     - Normalised Sum Square Residuals

% setup equations



% Prediction step
x_hat_minus = phi * x_hat_in;
P_minus = phi * P_in * phi' + Q;

% correction step

% calculate the kalman gain
V = H * P_minus * H' + R;
K = P_minus * H' * inv(V);

% update equations
%v = z - H * x_hat_minus;
v = z;
x_hat_out = x_hat_minus + K * (v);
P_out = P_minus - K * H * P_minus;


s2_out = (v' * V * v);

v_out = v;
