function [omega1, omega2, omega3] = vex(Q)

%inverse of the skew_symmetric function, returns the vector from the 3x3
%matrix Q
% Written by Troy B. 5.2.09
%



omega1 = Q(3,2);
omega2 = Q(1,3);
omega3 = Q(2,1);

