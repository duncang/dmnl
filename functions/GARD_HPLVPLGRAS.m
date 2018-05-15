function [HPL_H0 ,VPL_H0] = GARD_HPLVPLGRAS(N,UserPos, SVPos);

%Calculates HPL and VPL For GRAS, for ONE epoch only. 
% It is assumed that the correct SV positions and for the correct satellites are already known, ie above a certain mask angle.

%
% RAIM Availability Simulation - Written by Duncan Greer 17 Nov 2005
% Modified for use as a function by Troy Bruggemann 23 Nov 2005.


%SVPos, this is the satellites which are above a certain mask angle (eg 7.5 degrees)
    
       
    


% initialise storage vectors


    for k=1:N
        
        % calculate the geometry matrix - partial derivative unit vectors to each
        % SV
    
        [Az(k), El(k)] = AzEl(UserPos(1:3), [SVPos(k,1),SVPos(k,2),SVPos(k,3)]); %note, Az and El is in radians here.
        
            G(k,1) = cos(El(k)) * cos(Az(k)); 
            G(k,2) = cos(El(k)) * sin(Az(k));
            G(k,3) = sin(El(k));
            G(k,4) = 1;
       
    end
    
    % save number of SVs visible this epoch
    %N(Epoch) = i;

    % find the expected variance of satellite measurements
    var_i = 7.5 ^ 2;

    % form weighting matrix
    W_inv = eye(N) * var_i;
    W = inv(W_inv);

    % calculate S matrix
    S = inv(G' * W * G) * G' * W;
    
    % calculate s and d values (e,n,u,t)
    for k=1:N
        s_east_dot_var(k) = S(1,k)^2 * var_i;
        s_north_dot_var(k) = S(2,k)^2 * var_i;
        s_up_dot_var(k) = S(3,k)^2 * var_i;
        s_t_dot_var(k) = S(4,k)^2 * var_i;
        s_east_north_dot_var(k) = S(1,k) * S(2,k) * var_i;
    end
    
    d_east_2 = sum(s_east_dot_var);
    d_north_2 = sum(s_north_dot_var);
    d_up_2 = sum(s_up_dot_var);
    d_t_2 = sum(s_t_dot_var);
    d_EN_2 = sum(s_east_north_dot_var);
    
    d_major = sqrt(((d_east_2 + d_north_2) / 2) + sqrt(((d_east_2 - d_north_2) / 2) + (d_EN_2)));
    
    % calculate HPL and VPL values - see GRAS Mops Section 2.3.10.1
    K_ffmd = 5.8; % see GRAS Mops Section 2.3.11.5.2.1.2
    HPL_H0 = 10 * d_major;
    VPL_H0 = K_ffmd * sqrt(d_up_2);
    
    size(HPL_H0)
    
    HPL_H0 = abs(HPL_H0);  %um the resulting HPL is complex dunc says he just takes the absolute value..should it be complex? Check this...
    
    
    
    
    
    