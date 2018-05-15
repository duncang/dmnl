
function [Qpdf, xk_hat_minus,DeltaX,wk] = GARD_PF1GPSTBGauss(Qpdf_prev, xk_hat_prev,wkprev, UserPos,UserVel,N,PRMeasured,PRRateMeasured,SVPos, SVVel);

global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;
%Version 1.00
%Troy Bruggemann 30 June 2005

%Standard Particle Filter. Can be compaared with the EKF and the least squares.

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


%assume Vk and Wk is zero, ie no white random noise included in the measurements.




RangeNoiseVariance = 20;
Rk = eye(N);
Rk = Rk * RangeNoiseVariance;


zk = PRMeasured';


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


        %Correct for the earth rotation here
        %%it has been verified that you SUBTRACT the earth rotation
        %%correction (as calculated above)from the calculated PR's,(this is equivalent to adding the correction to the
        %%measured PR's)
        %%In other words, a positive delta_pr_omegaedot means that this needs to be added to the
        %%measured PR, which is the same as subtracting from the calculated PR
        %%(ResCalc(k))

        ResCalc(k) =  r_VecCalc(k) +  UserPos(4) +  UserVel(4)*dt;% + SVPos(k,4);% - delta_pr_omegaedot(k);                    
        
        
        
    end
    
    
        



    %Design Matrix
    %generate elements of M matrix

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

   
    
   
    
            
    

% Hk
% 
% 
% 
% ResVec = PRMeasured' - ResCalc';
% 
% ResVecVel = PRRateMeasured' - Relative_Velocity';
% 
% zk = [ResVec]
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


%  xk_hat_prev = xk_prev';





%PF algorithm starts here


%Filter via SIS


NSamp = 100;  %Number of particles



%the current state is xk_minus(k) , previous is xk_minus1



%The particle filter algorithms start here:


%generate Normal probability function 

%Assuming that sum of many errors leads to gaussian, it is safe to assume that position errors are normally distributed.

% 
% Y = NORMPDF(X,MU,SIGMA);
% 
% 
% sigma = 1; %standard deviation
% mu = 0; %mean
% 
%   x = [-1:0.001:1];
% for i = 1:2000
%     
%   
% 
% y(i) = (1/(sqrt(2*pi*sigma^2)))*exp(-(x(i)-mu)/(2*sigma^2));
% 
% Ybah = normpdf(1,0,1)
% 
% end




%These are the delta x's, delta y's, ie the estimated - true state
% [x y z dt xdot ydot zdot dtdot]
% 
% x = [0 0 0 0 0 0 0 0]

%Sampling of the random  vectors with the distribution P0







%---------------------------------------------------------------------------
%Choose importance density as suboptimal choice as transitional prior
%---------------------------------------------------------------------------
%Form importance density q and sample according to q



%For additive zero-mean Gaussian process noise see p 47 of 'Beyond the Kalman Filter' by Ristic, Arumlampalam and Gordon.

%p_xk_G_xk_prev(i) = NORMPDF(X,xk_hat_minus,Qk)    %_G_ represents 'given' _C_ represents 'comma'.  


%Should Qk vary with time? how can u calculate it?

% xstddev = sqrt(Qk(1,1)); 
% 
% ystddev = sqrt(Qk(2,2)) ;
% 
% zstddev = sqrt(Qk(3,3)) ;
% 
% dtstddev = sqrt(Qk(4,4));
% 
% xdotstddev = sqrt(Qk(5,5)) ;
% 
% ydotstddev = sqrt(Qk(6,6)) ;
% 
% zdotstddev = sqrt(Qk(7,7)) ;
% 
% dtdotstddev = sqrt(Qk(8,8));
% 
% 




%PREDICTION STAGE 
% note that using the transitional prior as the suboptimal choice of q means the weights have to be calculated AFTER the projection of the particles
%generate particles for the current time (k)
for i = 1:NSamp

    %project the particles ahead, assume white noise is zero or negligable
    
 
    
    xk_hat_minustemp = Phik*xk_hat_prev(i,1:8)';   
     
     
     for jj = 1:8
         xk_hat_minus(i,jj) = xk_hat_minustemp(jj);
     end
    
end






% %Sample pdf to make the current importance density, based upon previous
% %Choose importance density as suboptimal choice as transitional prior, this is what the book says is the most popular one.
% 
% P0 = Qpdf_prev
% 
% 
% for i = 1:NSamp      
% 
% 
% mean_x = xk_hat_minus(i,1);
% mean_y = xk_hat_minus(i,2);
% mean_z = xk_hat_minus(i,3);
% mean_dt = xk_hat_minus(i,4);
% 
% mean_xvel = xk_hat_minus(i,5);
% mean_yvel = xk_hat_minus(i,6);
% mean_zvel = xk_hat_minus(i,7);
% mean_dtvel = xk_hat_minus(i,8);
% 
% %not sure if they should use the same random sample number for all x,y,z..etc or not. I've made it different for each at the moment.
% 
% q_xpos = normrnd(mean_x,xstddev) ;
% q_ypos = normrnd(mean_y,ystddev) ;
% q_zpos = normrnd(mean_z,zstddev) ;
% q_dtpos = normrnd(mean_dt,dtstddev) ;
% 
% q_xvel = normrnd(mean_xvel,xdotstddev) ;
% q_yvel = normrnd(mean_yvel,ydotstddev) ;
% q_zvel = normrnd(mean_zvel,zdotstddev) ;
% q_dtvel = normrnd(mean_dtvel,dtdotstddev) ;
% 
% 
% Q_xk_G_xk_prev_C_zk(i,1:8) = [q_xpos q_ypos q_zpos q_dtpos q_xvel q_yvel q_zvel q_dtvel];
% 
% end
P0 = Qpdf_prev;

%Sample pdf to make the current importance density, based upon previous
%Choose importance density as Gaussian Optimal Importance Function (p 45-46). 


for i = 1:NSamp
    
bk = Hk*xk_hat_minus(i,1:8)';            %don't know which 'i' to use here..bk is an N*1 i think

Sk = Hk*Qk*Hk' + Rk'  ;                  %Sk is an NxN, Hk is an N*8...remember inner matrix dimensions must be the same, and the resulting matrix has the dimension of the outer dimensions
Sigmak = Qk - Qk*Hk'*inv(Sk)*Hk*Qk;     %Sigma k is an 8x8

% Sigmak*Hk'*inv(Rk)*(zk - bk); %just to view the output
% 
% xk_hat_minus(i,1:8) %just to view the output

ak = xk_hat_minus(i,1:8)' + Sigmak*Hk'*inv(Rk)*(zk - bk);    %ak is an 8x1 , zk must be an Nx1 and bk must be an Nx1

%note that zk - bk is like the measured minus predicted pseudoranges




xstddev = sqrt(Sigmak(1,1)); 

ystddev = sqrt(Sigmak(2,2)) ;

zstddev = sqrt(Sigmak(3,3)) ;

dtstddev = sqrt(Sigmak(4,4));

xdotstddev = sqrt(Sigmak(5,5)) ;

ydotstddev = sqrt(Sigmak(6,6)) ;

zdotstddev = sqrt(Sigmak(7,7)) ;

dtdotstddev = sqrt(Sigmak(8,8));


mean_x = ak(1);
mean_y = ak(2);
mean_z = ak(3);
mean_dt = ak(4);

mean_xvel = ak(5);
mean_yvel = ak(6);
mean_zvel = ak(7);
mean_dtvel = ak(8);

%not sure if they should use the same random sample number for all x,y,z..etc or not. I've made it different for each at the moment.

q_xpos = normrnd(mean_x,xstddev) ;
q_ypos = normrnd(mean_y,ystddev) ;
q_zpos = normrnd(mean_z,zstddev) ;
q_dtpos = normrnd(mean_dt,dtstddev) ;

q_xvel = normrnd(mean_xvel,xdotstddev) ;
q_yvel = normrnd(mean_yvel,ydotstddev) ;
q_zvel = normrnd(mean_zvel,zdotstddev) ;
q_dtvel = normrnd(mean_dtvel,dtdotstddev) ;


Q_xk_G_xk_prev_C_zk(i,1:8) = [q_xpos q_ypos q_zpos q_dtpos q_xvel q_yvel q_zvel q_dtvel];






% to find the pdf p_zk_G_xkminusone, for the weighting matrix



    for k = 1:N
        
        %note that sqrt(Sk(k,k)) is the standard deviation. I've ignored cross correlation components
       p_zk_G_xkminusone(i,k) =  normrnd(bk(k),sqrt(Sk(k,k))); %not sure if they should use the same random sample number for all x,y,z..etc or not. I've made it different for each at the moment.


   end





wkprev(i,1:N)

%I need to work out what size this weighting matrix has to be
    
wk_est = wkprev(i,1:N).*p_zk_G_xkminusone(i,1:N)  %equation 3.39    %need to do .* because want to multiply each element by the same element in the other, ie weight each one, not sum them as in normal matrix multiplication
    
%         
% 
% %Calculate total weight
% totalW = sum(wk_est)
% %normalise wk
% 
% 
%  wk(i,1:N) = wk_est(i,:)/totalW;
% 
% %wk(i) = 1%wk_est/totalW;  %There is only one weight per particle!, ie weight is a 1xNSamp matrix or NSampx1 , depending, %total sum of the weights (sum wk) should be equal to 1! 
% 
% 

wk(i,1:N) = wkprev(i,1:N)


end



% 
% 
% 
% %Calculate Neff
% 
% for i = 1:NSamp
% 
% Neff(i) = sum(wk(i,1:N)^2)
% 
% end
% 
% 
% %Nthreshold is?
% Nthreshold = 100;
% 
% 
% if Neff > Nthreshold
%     
%     Resample = 1;
%     
% end
% 
% 
% if Resample == 1
%     
%     
%     %do resampling procedure
% 
%     
%     
% end
% 
% 





Qpdf = Q_xk_G_xk_prev_C_zk;



DeltaX = Q_xk_G_xk_prev_C_zk(1,1:8);




%pull the best particles out for the xk

















    



























