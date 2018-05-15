function q = GARD_DCMToQuat(C_BN)
% Usage: q = GARD_DCMToQuat(C_BN)
%
% From Titterton p45 and http://planning.cs.uiuc.edu/node153.html
% Updated method by Sheppard
% $Id$

M(1) = trace(C_BN);
M(2) = C_BN(1,1);
M(3) = C_BN(2,2);
M(4) = C_BN(3,3);

[maxval,index] = max(M);

p(index) = sqrt(1 + 2*M(index) - M(1));

switch index
    case 1
        p(2) = (C_BN(3,2) - C_BN(2,3))/p(1);
        p(3) = (C_BN(1,3) - C_BN(3,1))/p(1);
        p(4) = (C_BN(2,1) - C_BN(1,2))/p(1);
        
        
    case 2
        p(1) = (C_BN(3,2) - C_BN(2,3))/p(2);
        p(3) = (C_BN(1,3) - C_BN(3,1))/p(1);
        p(4) = (C_BN(2,1) - C_BN(1,2))/p(1);
    case 3
        p(1) = (C_BN(1,3) - C_BN(3,1))/p(3);
        p(2) = (C_BN(3,2) - C_BN(2,3))/p(1);
        p(4) = (C_BN(2,1) - C_BN(1,2))/p(1);
    case 4
        p(1) = (C_BN(2,1) - C_BN(1,2))/p(4);
        p(2) = (C_BN(3,2) - C_BN(2,3))/p(1);
        p(3) = (C_BN(1,3) - C_BN(3,1))/p(1);        
    otherwise
        warning('bad index');
end

q = p/2.0;

% q(1) = 0.5 * sqrt(1.0 + C_BN(1,1) + C_BN(2,2) + C_BN(3,3));
% 
% if q(1) == 0
%     blah = sqrt(C_BN(1,2)^2 * C_BN(1,3)^2 + C_BN(1,2)^2*C_BN(2,3)^2 + C_BN(1,3)^2*C_BN(2,3)^2);
%     
%     q(2) = C_BN(1,3)*C_BN(1,2) / blah;
%     q(3) = C_BN(1,2)*C_BN(2,3) / blah;
%     q(4) = C_BN(1,3)*C_BN(2,3) / blah;
%     
% else
%     q(2) = (1.0/(4*q(1))) * (C_BN(3,2) - C_BN(2,3));
%     q(3) = (1.0/(4*q(1))) * (C_BN(1,3) - C_BN(3,1));
%     q(4) = (1.0/(4*q(1))) * (C_BN(2,1) - C_BN(1,2));
% end
% 
% % check for norm
% if norm(q) < 0.9999 || norm(q) > 1.0001
%     disp('Warning - quaternion not unity');
% end