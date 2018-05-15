function [Q] = skewbody(p,q,r)
% function [Q] = skewbody(q)
% Returns the skew form of body rates p,q,r
%
% Written by Duncan Greer 26 Jan 2007
%
% $Id: skewbody.m 1863 2008-07-14 07:02:29Z greerd $
%
%


Q = zeros(4,4);

Q(1,1) =  0.0; Q(1,2) = p;   Q(1,3) = q;    Q(1,4) = r;
Q(2,1) =  -p;  Q(2,2) = 0.0; Q(2,3) = -r;   Q(2,4) = q;
Q(3,1) =  -q;  Q(3,2) = r;   Q(3,3) =  0.0; Q(3,4) = -p;
Q(4,1) =  -r;  Q(4,2) = -q;  Q(4,3) =  p;   Q(4,4) = 0.0;