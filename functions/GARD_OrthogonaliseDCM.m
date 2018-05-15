function C_BN_out = GARD_OrthogonaliseDCM(C_BN_in)
% function C_BN_out = GARD_OrthogonaliseDCM(C_BN_in)
% Written by Duncan Greer 27 November 2007
% This function normalises the Direction Cosine Matrix
% using Bar-Ithacks Dual Method described in "Orthogonalization of
% Direction Cosine by Iterative Techniques" IEEE 1972.
%
% $Id: GARD_OrthogonaliseDCM.m 1850 2008-07-14 04:52:47Z greerd $
%
C_BN_n = C_BN_in;
iter = 0;
for i=1:10
    iter = iter+1;
    %% original method
    %C_BN_n = 1.5*C_BN_n - 0.5*C_BN_n*C_BN_n'*C_BN_n;
    
    
    %% dual method
    C_BN_n = 0.5 * inv(C_BN_n') + 0.5*C_BN_n;
        
    
    %% gradient projection method
    %C_BN_n = 0.5*C_BN_in + C_BN_n - 0.5*C_BN_n*C_BN_in'*C_BN_n;
    
    residual = eye(3,3) - C_BN_n * C_BN_n';
    
    if norm(residual) < 1e-15
        
        break;
    end
    
end

if(iter == 10)
    disp('[GARD_OrthogonaliseDCM: Warning - DCM did not converge');
end


C_BN_out = C_BN_n;