
function [SolutionVec, VarSolutionVec, NumIterations, ResidualVector, M,LSQ_Fail, limit] = GARD_LSQVelUDS(UserPos,UserVel,N,PRMeasured,SVPos, SVVel);

global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;
%Version 1.00
%Troy Bruggemann 30 June 2005

%This function does a least squares code solution for receiver velocities with an underdetermined solution, single epoch only.
%This works by using an estimated receiver clock drift calculated by
%differentiating the receiver clock bias
%this means only 3 range rate measurements are required to solve for
%x_dot,y_dot and z_dot because the dt_dot term is estimated

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

 %start least squares
 
 %initialise variables
    LSQ_Fail = 0;  
    
    %Initial user guess
    
 
    
    Xb = UserPos(1);
    Yb = UserPos(2);
    Zb = UserPos(3);
    dtb = UserPos(4); %This is the estimated receiver clock drift from differentiating the receiver clock bias
    
    XbVel = UserVel(1);
    YbVel = UserVel(2);
    ZbVel = UserVel(3);
    %dtbVel = UserVel(4);
    
      
    %TempBase_B_Vel = [XbVel,YbVel,ZbVel, dtbVel];
    TempBase_B_Vel = [XbVel,YbVel,ZbVel];
    
    
    TempBase_B_Pos = [Xb,Yb,Zb, dtb];
    
    N = 3 %only use 3 measurements to calculate velocity.
    
    
%     for k = 1:N   
%     
%     % calculate hte earth rotation correction as per Kayton pg 228
%     % eq 5.67    
%     
%     dt = PRMeasured(k)/Speedoflight;    
%     %delta_pr_omegaedot = -(OMEGAedot / c) * (SVPos(k,1) * x(2) - SVPos(k,2) * x(1));    
%     SV_Pos_Corrx = -OMEGAedot*dt*SVPos(k,2);
%     SV_Pos_Corry =  OMEGAedot*dt*SVPos(k,1);
% 
%     SVPos(k,1) = SVPos(k,1) + SV_Pos_Corrx;
%     SVPos(k,2) = SVPos(k,2) + SV_Pos_Corry;
% 
%     
% 
%     end

    
%Start loop here, loop for LSQ
for j = 1:200   %Least square iterations counter

%Calculated slant ranges

for k = 1:N
    
    
    for m = 1:3
         ele(m) =  SVPos(k,m) - TempBase_B_Pos(m);
    end    
    
   r_VecCalc(k) =  norm(ele);   


end


%Calculated relative velocities

for k = 1:N
    
%     for m = 1:3
%          ele(m) =  SVVel(k,m) - TempBase_B_Vel(m);
%     end    
%     
%    r_VecCalcVel(k) =  norm(ele);   


%predicted relative velocity of sv and receiver

 r_VecCalcVel(k) = (SVVel(k,1) - TempBase_B_Vel(1))*(SVPos(k,1)-TempBase_B_Pos(1)) + (SVVel(k,2) - TempBase_B_Vel(2))*(SVPos(k,2)-TempBase_B_Pos(2)) + (SVVel(k,3) - TempBase_B_Vel(3))*(SVPos(k,3)-TempBase_B_Pos(3));

% velocities = (GPS_SatVel[i][0] - VelX[0])*(GPS_SatPos[i][0] - X[0]) 
% 					+ (GPS_SatVel[i][1] - VelX[1])*(GPS_SatPos[i][1] - X[1]) 
% 					+ (GPS_SatVel[i][2] - VelX[2])*(GPS_SatPos[i][2] - X[2]);
    

Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);



end





%Use SP3 values for satellite clock error:
% 1     3     4    13    16    19    20    23    24    27

% 
%  SVClkSP3 = Speedoflight*[385.268028e-6 80.159890e-6 427.043212e-6 -13.131507e-6 1.878373e-6 -12.495888e-6  -88.990365e-6 211.762141e-6  91.458908e-6  517.820742e-6];



       



%Correct calculated slant ranges with receiver clock bias and satellite
%clock bias and earth rotation correction
 


%  delta_t = PR / c;
% 
%   SV_Pos_Corr[0] = -w_erde * delta_t * GPS_SV_Pos[1];
%   SV_Pos_Corr[1] =  w_erde * delta_t * GPS_SV_Pos[0];
%   SV_Pos_Corr[2] =  0.0;

for k = 1:N   
    
    % calculate hte earth rotation correction as per Kayton pg 228
    % eq 5.67    
    
% dt = PRMeasured(k)/Speedoflight;    
% %delta_pr_omegaedot = -(OMEGAedot / c) * (SVPos(k,1) * x(2) - SVPos(k,2) * x(1));    
% SV_Pos_Corrx = -OMEGAedot*dt*SVPos(k,2);
% SV_Pos_Corry =  OMEGAedot*dt*SVPos(k,1);
% 
% SVPos(k,1) = SVPos(k,1) - SV_Pos_Corrx;
% SVPos(k,2) = SVPos(k,2) - SV_Pos_Corry;

%ResCalc(k) =  r_VecCalc(k)  + TempBase_B_Pos(4)- SVClkSP3(k);  

% TempBase_B_Pos(4) is the estimated receiver clock drift
ResCalcVel(k) =  Relative_Velocity(k) + TempBase_B_Pos(4) - SVVel(k,4) ;  %462707.565238509


 
%ResCalc(k) =  r_VecCalc(k) - SVPos(k,4);
%  
%    ResCalc(k) =  r_VecCalc(k) ;%+ SVPos(k,4) ;  
end



%Design Matrix
%generate elements of M matrix

for k = 1:N
    
            
        M(k,1) =  -(SVPos(k,1) - TempBase_B_Pos(1))/r_VecCalc(k);
        M(k,2) =  -(SVPos(k,2) - TempBase_B_Pos(2))/r_VecCalc(k);
        M(k,3) =  -(SVPos(k,3) - TempBase_B_Pos(3))/r_VecCalc(k);
        %M(k,4) = 1;   
        
        
                
end
    

%residual vector
ResVec = PRMeasured(1:3)' - ResCalcVel';



%LSQ solution

A = M'*M;
b = M'*ResVec;

Delta_X = inv(A)*b;
%Delta_X = A\b;

% limit = sqrt(Delta_X(1)^2 +  Delta_X(2)^2  + Delta_X(3)^2 +Delta_X(4)^2);

limit = sqrt(Delta_X(1)^2 +  Delta_X(2)^2  + Delta_X(3)^2 );
numlimitAnt1_2(j) = limit;
if limit < 0.00000001;
    numiterAnt1_2 = j;
    break;
else
    numiterAnt1_2 = 200;
end

if j == 200
    LSQ_Fail = 1;  %Solution Cannot converge
end

TempBase_B_Vel(1) = TempBase_B_Vel(1) + Delta_X(1);
TempBase_B_Vel(2) = TempBase_B_Vel(2) + Delta_X(2);
TempBase_B_Vel(3) = TempBase_B_Vel(3) + Delta_X(3);
%TempBase_B_Vel(4) = TempBase_B_Vel(4) + Delta_X(4);
                                        



                                        


end

%final solution
XbVel = TempBase_B_Vel(1);
YbVel = TempBase_B_Vel(2);
ZbVel = TempBase_B_Vel(3);
dtbVel = TempBase_B_Pos(4);  %this is not calculated by least squares, but we've estimated it from differentiation receiver clock bias
 


%DOPvalues
var_x = A(1,1);
var_y = A(2,2);
var_z = A(3,3);
% var_dt = A(4,4);

%Outputs
SolutionVec = [XbVel,YbVel,ZbVel,dtbVel];
VarSolutionVec = [var_x,var_y,var_z];
NumIterations = numiterAnt1_2;
ResidualVector = ResVec ; %Last iteration residual vector, output for RAIM using Least squares residual method





