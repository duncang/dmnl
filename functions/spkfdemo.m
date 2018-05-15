%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UNSCENTED TRANSFORM DEMO
%%% Written by Duncan Greer (c) 2007)
%%% Ref: Beyond The Kalman Filter, Ristic et al.  pp30-31
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clear workspace
clear;

% setup true rv
x_bar = [0;0];
P_x = [1 0.5; 0.5 1];

% calculate sigma points
sqrtPx = sqrtm(P_x);
xs0 = x_bar;
xs11 = x_bar + sqrtPx(1,:)';
xs12 = x_bar - sqrtPx(1,:)';
xs21 = x_bar + sqrtPx(2,:)';
xs22 = x_bar - sqrtPx(2,:)';

W0 = -1;
W1 = 0.5;
W2 = 0.5;

% propogate through system dynamics - in this case 1:1
bs0 = xs0;
bs11 = xs11;
bs12 = xs12;
bs21 = xs21;
bs22 = xs22;

% re-form rv from sigma points
b_bar = [0;0];
b_bar = b_bar + W0*bs0 + W1*bs11 + W1 * bs12 + W2*bs21 + W2 * bs22;
P_b = zeros(2,2);
P_b = P_b + W0 * ((bs0 - b_bar) * (bs0 - b_bar)');
P_b = P_b + W1 * ( (bs11 - b_bar) * (bs11 - b_bar)' );
P_b = P_b + W1 * ( (bs12 - b_bar) * (bs12 - b_bar)' );
P_b = P_b + W2 * ( (bs21 - b_bar) * (bs21 - b_bar)' );
P_b = P_b + W2 * ( (bs22 - b_bar) * (bs22 - b_bar)' );


% plot results
error_ellipse(P_x)
hold on;
grid on;
plot(x_bar(1),x_bar(2),'b+')
plot(xs11(1),xs11(2),'bd')
plot(xs12(1),xs12(2),'bd')
plot(xs21(1),xs21(2),'bd')
plot(xs22(1),xs22(2),'bd')
plot(b_bar(1),b_bar(2),'ko');
plot(bs11(1),bs11(2),'k+');
plot(bs12(1),bs12(2),'k+');
plot(bs21(1),bs21(2),'k+');
plot(bs22(1),bs22(2),'k+');
error_ellipse(P_b);
xlabel('x[1]');
ylabel('x[2]');
title('Selection of Sigma Points for 2-dimensional Random Variable');



