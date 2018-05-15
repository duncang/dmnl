function [Q] = skewquat(q)
% function [Q] = skewquat(q)
% Returns the skew form of a quaternion
%
% Written by Duncan Greer 26 Jan 2007
%
% $Id: skewquat.m 1863 2008-07-14 07:02:29Z greerd $
%
%


Q = zeros(4,4);

Q(1,1) = q(1);
Q(1,2) = -q(2);
Q(1,3) = -q(3);
Q(1,4) = -q(4);

Q(2,1) = q(2);
Q(2,2) = q(1);
Q(2,3) = -q(4);
Q(2,4) = q(3);

Q(3,1) = q(3);
Q(3,2) = q(4);
Q(3,3) = q(1);
Q(3,4) = -q(2);

Q(4,1) = q(4);
Q(4,2) = -q(3);
Q(4,3) = q(2);
Q(4,4) = q(1);