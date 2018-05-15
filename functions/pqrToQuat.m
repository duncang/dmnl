function [quat_dot] = pqrToQuat(pqr, quat)

%
%
% $Id: pqrToQuat.m 1883 2008-07-15 05:53:55Z n2523710 $
%

%from eq 1.4-28 page 43 of Aircraft Control and Simulation (Stevens)

p = pqr(1);
q = pqr(2);
r = pqr(3);


q0 = quat(1);
q1 = quat(2);
q2 = quat(3);
q3 = quat(4);


Q_dot = -0.5*[0, p, q, r; -p, 0, -r, q; -q, r, 0, -p; -r, -q, p, 0;]*[q0;q1;q2;q3];     





quat_dot(1) = Q_dot(1);
quat_dot(2) = Q_dot(2);
quat_dot(3) = Q_dot(3);
quat_dot(4) = Q_dot(4);


