
function [xk_hat, Pk] = GARD_KFGPSTB(xhat_prev,UserPos,UserVel,Pk_prev,N,PRMeasured,PRRateMeasured,SVPos, SVVel);

global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;
%Version 1.00
%Troy Bruggemann 30 June 2005

%This function does a least squares code solution for receiver velocities, single epoch only.

%INPUTS
% UserPos - [XPos,YPos,ZPos,dTPos] Estimated User Navigation State Vector(m)
% N - number of observations to use in solution (number of satellites)
% PRMeasured - Vector of Measured Pseudoranges [1..N]
% SVPos - Matrix of Satellite poisitions [1..N][Xs,Ys,Zs,dTs] (m)
%========================================================================
% OUTPUTS
%SolutionVec = [X,Y,Z,dt];
% X - estimated position ECEF (m)
% Y - estimated position ECEF (m)
% Z - estimated position ECEF (m)
% dt - estimated receiver clock bias (m)
%
%VarSolutionVec = [var_x,var_y,var_z,var_dt];
% var_x - estimated X position variance (m)
% var_y - estimated X position variance (m)
% var_z - estimated X position variance (m)
% var_dt - estimated X position variance (m)
% NumIterations - Number of iterations used for the solution
% ResidualVector - [1..N] Vector of pseudorange residuals (m)
% M - Least squares design matrix [1..M][1..4]

%Constants
% Speedoflight = 2.99792458e8; %m/s
% L1_Freq = 1575.42e6; %Hz
% L1_Wavelength = Speedoflight/L1_Freq; %Metres
%
% OMEGAedot = 7.2921151467e-5; % Earth rotation rate in radians per second



% setup error values
dt = 1;


% not sure what these values are - something to do with the process noise
Sp = 1.0;
Sf = 1.0;
Sg = 1.0;

% setup state transition - time invariant if dt is constant
Phik = eye(8,8);
Phik(1,5) = dt;
Phik(2,6) = dt;
Phik(3,7) = dt;
Phik(4,8) = dt;


% setup Q and R matrices (process and noise covariance)
Qk = zeros(8,8);

Sp_dt3 = Sp * (dt ^ 3) / 3;
Sp_dt2 = Sp * (dt ^ 2) / 2;
Sp_dt = Sp * dt;
Sf_dt = Sf * dt;
Sg_dt3 = Sg * (dt ^ 3) / 3;
Sg_dt2 = Sg * (dt ^ 2) / 2;
Sg_dt = Sg * dt;


% %form diagonal elements
% %note that this Q matrix is for state vector with states in the following order: [x, y , z , x_dot, y_dot, z_dot, dt, dt_dot]
% 


%from Peter's solution
% Qk(1,1) = Sp_dt3;
% 
% Qk(2,3) = Qk(1,1);
% Qk(3,5) = Qk(1,1);
% 
% Qk(1,2) = Sp_dt2;
% 
% Qk(2,4) = Sp_dt2;
%  
% Qk(3,6) = Sp_dt2;
%  
% Qk(5,1) = Sp_dt2;
% Qk(6,3) = Sp_dt2;
% Qk(7,5) = Sp_dt2;
% 
% 	
% Qk(5,2) = Sp_dt;
% Qk(6,4) = Sp_dt;   
% Qk(7,6) = Sp_dt;
%    
% Qk(4,7) = Sf_dt + Sg_dt3;
% Qk(4,8) =Sg_dt2;     
% Qk(8,7) = Sg_dt2;
% Qk(8,8) = Sg_dt;
% 












% % 
% Qk(1,1) = Sp_dt3; %the x's
% Qk(1,4) = Sp_dt3;
% Qk(4,1) = Qk(1,4);
% Qk(4,4) = Sp_dt;
% 
% 
% Qk(2,2) = Sp_dt3; %the y's
% Qk(2,5) = Sp_dt3;
% Qk(5,2) = Qk(2,5);
% Qk(5,5) = Sp_dt;
% 
% Qk(3,3) = Sp_dt3; %the z's
% Qk(3,6) = Sp_dt3;
% Qk(6,3) = Qk(3,6);
% Qk(6,6) = Sp_dt;
% 
% Qk(7,7) = Sf_dt + Sg_dt3; %the dt's
% Qk(7,8) = Sg_dt2;
% Qk(8,7) = Qk(7,8);
% Qk(8,8) = Sg_dt;



%duncans

Qk(1,1) = Sp_dt3;
Qk(2,2) = Sp_dt3;
Qk(3,3) = Sp_dt3;

Qk(1,5) = Sp_dt2;
Qk(2,6) = Sp_dt2;
Qk(3,7) = Sp_dt2;
Qk(4,8) = Sg_dt2;
Qk(5,1) = Sp_dt2;
Qk(6,2) = Sp_dt2;
Qk(7,3) = Sp_dt2;
Qk(8,4) = Sg_dt2;

Qk(4,4) = Sf_dt + Sg_dt3;
Qk(5,5) = Sp_dt;
Qk(6,6) = Sp_dt;
Qk(7,7) = Sp_dt;
Qk(8,8) = Sg_dt;







RangeNoiseVariance = 1.0;
    Rk = eye(N);
    Rk = Rk * RangeNoiseVariance;



%================================
%Kalman filter starts here
%================================
%this implementation is based upon An Introduction to the Kalman Filter
% Greg Welch
% 1
% and Gary Bishop
% 2


%Construct matrices and such



for k = 1:N

        for m = 1:3
            ele(m) =  SVPos(k,m) - UserPos(m);
        end

        r_VecCalc(k) =  norm(ele);

    end



    for k = 1:N

        % calculate hte earth rotation correction as per Kayton pg 228
        % eq 5.67

        delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) * UserPos(2) - SVPos(k,2) * UserPos(1));

    end

    for k = 1:N


        %Correct for the earth rotation here
        %%it has been verified that you SUBTRACT the earth rotation
        %%correction (as calculated above)from the calculated PR's,(this is equivalent to adding the correction to the
        %%measured PR's)
        %%In other words, a positive delta_pr_omegaedot means that this needs to be added to the
        %%measured PR, which is the same as subtracting from the calculated PR
        %%(ResCalc(k))

        ResCalc(k) =  r_VecCalc(k) +  UserPos(4) +  UserVel(4)*dt;% + SVPos(k,4);% - delta_pr_omegaedot(k);

  
        
         %predicted relative velocity of sv and receiver

        r_VecCalcVel(k) = (SVVel(k,1) - UserVel(1))*(SVPos(k,1)-UserPos(1)) + (SVVel(k,2) - UserVel(2))*(SVPos(k,2)-UserPos(2)) + (SVVel(k,3) - UserVel(3))*(SVPos(k,3)-UserPos(3));


        Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);

        
        
        
        
        
    end
    
    
        



    %Design Matrix
    %generate elements of M matrix

%     for k = 1:N
%         Hkk(k,1) =  -(SVPos(k,1) - UserPos(1))/r_VecCalc(k);
%         Hkk(k,2) =  -(SVPos(k,2) - UserPos(2))/r_VecCalc(k);
%         Hkk(k,3) =  -(SVPos(k,3) - UserPos(3))/r_VecCalc(k);
%         Hkk(k,4) = -(SVPos(k,1) - UserPos(1))/r_VecCalc(k);
%         Hkk(k,5) = -(SVPos(k,2) - UserPos(2))/r_VecCalc(k);
%         Hkk(k,6) = -(SVPos(k,3) - UserPos(3))/r_VecCalc(k);
%         Hkk(k,7) = 1;
%         Hkk(k,8) = 1;
%                
%         
%     end

    for k = 1:N
        Hk(k,1) =  -(SVPos(k,1) - UserPos(1))/r_VecCalc(k);
        Hk(k,2) =  -(SVPos(k,2) - UserPos(2))/r_VecCalc(k);
        Hk(k,3) =  -(SVPos(k,3) - UserPos(3))/r_VecCalc(k);
        Hk(k,4) = 1;
        Hk(k,5) = 0;
        Hk(k,6) = 0;
        Hk(k,7) = 0;
        Hk(k,8) = 0;
               
        
    end


%     for k = 1:N
%         Hk(k,1) =  -(SVPos(k,1) - UserPos(1))/r_VecCalc(k);
%         Hk(k,2) =  -(SVPos(k,2) - UserPos(2))/r_VecCalc(k);
%         Hk(k,3) =  -(SVPos(k,3) - UserPos(3))/r_VecCalc(k);
%         Hk(k,4) = (r_VecCalc(k)^2 - (SVPos(k,1) - UserPos(1))^2) / r_VecCalc(k)^3;
%         Hk(k,5) = (r_VecCalc(k)^2 - (SVPos(k,2) - UserPos(2))^2) / r_VecCalc(k)^3;
%         Hk(k,6) = (r_VecCalc(k)^2 - (SVPos(k,3) - UserPos(3))^2) / r_VecCalc(k)^3;
%         Hk(k,7) = 1;
%         Hk(k,8) = 1;       
%         
%     end


    
    
    
    %Hk = cat(1,Hkk,Hkk);
    
    
    
            
    

Hk



ResVec = PRMeasured' - ResCalc';

ResVecVel = PRRateMeasured' - Relative_Velocity';

zk = [ResVec]
%zk = [ResVec; ResVecVel]

%Qk = zeros(8,8);
%Rk = eye(N*2,N*2);


%xk is process state vector at time k
%Phik is the state transition matrix
% Qk is covariance matrix of the process noise
% wk is the process white noise vector
% 
% Hk is the relationship between xk and zk
% zk is the measurement at time k
% vk is the measurement white noise vector
% Rk is the measurement noise covariance matrix


%initial update



Xb = UserPos(1);
Yb = UserPos(2);
Zb = UserPos(3);
dtb = UserPos(4);

XbVel = UserVel(1);
YbVel = UserVel(2);
ZbVel = UserVel(3);
dtbVel = UserVel(4);


xk_hat_prev = xhat_prev';




%the current state is xk_minus(k) , previous is xk_minus1


    
    %Project error covariance ahead
    Pk_minus = Phik*Pk_prev*Phik' + Qk  %assume Wk not exists

     a = 50.5e-9; b = -50.5e-9;  %range of values for random number %note. b has to be the negative value one.
% 
% RandomNoise = a + (b-a) * rand(1);
    
    
    %compute kalman gain
    
    
    Kk = Pk_minus*Hk'*inv(Hk*Pk_minus*Hk' + Rk)  %assume Vk not exists



    %update estimate with measurement
   
    
    xk_hat = xk_hat_prev + Kk*(zk - Hk*xk_hat_prev)
    
    %update error covariance
    
    %Pk = (eye(8,8) - Kk*Hk)*Pk_minus;
      Pk = Pk_minus - Kk * Hk * Pk_minus;





xk_hat























