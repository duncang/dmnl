function [SolutionVec, VarSolutionVec, NumIterations, ResidualVector, M, LSQ_Fail, limit, DOP] = GARD_LSQClockCoast(UserPos,N,PRMeasured,SVPos,EstimatedClockBiasCoast,ClockAugmentation,k_squared);

%Version 1.00
%Troy Bruggemann 30 June 2005
%This function does a least squares code solution, single epoch only.
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
%
%DOPS
%put the DOPS in a vector, the order is [GDOP PDOP HDOP VDOP]
%
% $Id: GARD_LSQClockCoast.m 1850 2008-07-14 04:52:47Z greerd $
%

global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;

%start least squares


%Nclock = N+1; %increase N by one to incorporate the clock measurement
%assume that PRMeasured already has the clock term in it

%initialise variables
LSQ_Fail = 0;

%Initial user guess



ScaleFactor = sqrt(k_squared); %see p 228 of Lee, i believe this is sigma_clock/sigmaPr after comparison with Sturza


if ClockAugmentation == 0  %augmenting accordign to Lee, just add another time in there.


Xb = UserPos(1);   %just corrupt my position guess a bit, shouldnt matter though
Yb = UserPos(2);
Zb = UserPos(3);
dtb = UserPos(4);

TempBase_B_Pos = [Xb,Yb,Zb, dtb];

end;


if ClockAugmentation == 1 %augmenting by using 3 satellite navigation

Xb = UserPos(1);   %just corrupt my position guess a bit, shouldnt matter though
Yb = UserPos(2);
Zb = UserPos(3);
dtb = EstimatedClockBiasCoast; %assume i know this

TempBase_B_Pos = [Xb,Yb,Zb, dtb];

end;






%Start loop here, loop for LSQ
for j = 1:200   %Least square iterations counter
    
    
    
  if ClockAugmentation == 0
    
    %Calculated slant ranges

    for k = 1:N

        for m = 1:3
            ele(m) =  SVPos(k,m) - TempBase_B_Pos(m);
        end

        r_VecCalc(k) =  norm(ele);

    end



    for k = 1:N

        % calculate hte earth rotation correction as per Kayton pg 228
        % eq 5.67

        delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) * TempBase_B_Pos(2) - SVPos(k,2) * TempBase_B_Pos(1));

    end

    for k = 1:N


        %Correct for the earth rotation here
        %%it has been verified that you SUBTRACT the earth rotation
        %%correction (as calculated above)from the calculated PR's,(this is equivalent to adding the correction to the
        %%measured PR's)
        %%In other words, a positive delta_pr_omegaedot means that this needs to be added to the
        %%measured PR, which is the same as subtracting from the calculated PR
        %%(ResCalc(k))

        ResCalc(k) =  r_VecCalc(k) + TempBase_B_Pos(4) - SVPos(k,4) - delta_pr_omegaedot(k);

        %ResCalc(k) =  r_VecCalc(k) - SVPos(k,4);
        %
        %    ResCalc(k) =  r_VecCalc(k) ;%+ SVPos(k,4) ;
    end
    
    %ResCalc(Nclock) = TempBase_B_Pos(4);  %add receiver clock term, previous estimate when integrity is assured.        

    
   % ResCalc(Nclock) = TempBase_B_Pos(4);

    %Design Matrix
    %generate elements of M matrix    
   

    for k = 1:N
        M(k,1) =  -(SVPos(k,1) - TempBase_B_Pos(1))/r_VecCalc(k);
        M(k,2) =  -(SVPos(k,2) - TempBase_B_Pos(2))/r_VecCalc(k);
        M(k,3) =  -(SVPos(k,3) - TempBase_B_Pos(3))/r_VecCalc(k);
        M(k,4) = 1;
    end
    
%     %add the clock 
%     M(Nclock,1) = 0;
%       M(Nclock,2) = 0;
%         M(Nclock,3) = 0;
%           M(Nclock,4) = 1/ScaleFactor;
%           
    end%      if ClockAugmentation == 0
    
    
    
    
    
          
   %if atomic clock augmentation = true , only ahve to solve for  3 unknowns.     
    if ClockAugmentation == 1
        
       
    %Calculated slant ranges

    for k = 1:N

        for m = 1:3
            ele(m) =  SVPos(k,m) - TempBase_B_Pos(m);
        end

        r_VecCalc(k) =  norm(ele);

    end


    for k = 1:N

        % calculate hte earth rotation correction as per Kayton pg 228
        % eq 5.67

        delta_pr_omegaedot(k) = -(OMEGAedot / Speedoflight) * (SVPos(k,1) * TempBase_B_Pos(2) - SVPos(k,2) * TempBase_B_Pos(1));

    end

    for k = 1:N


        %Correct for the earth rotation here
        %%it has been verified that you SUBTRACT the earth rotation
        %%correction (as calculated above)from the calculated PR's,(this is equivalent to adding the correction to the
        %%measured PR's)
        %%In other words, a positive delta_pr_omegaedot means that this needs to be added to the
        %%measured PR, which is the same as subtracting from the calculated PR
        %%(ResCalc(k))

        ResCalc(k) =  r_VecCalc(k) + EstimatedClockBiasCoast - SVPos(k,4) - delta_pr_omegaedot(k);

        %ResCalc(k) =  r_VecCalc(k) - SVPos(k,4);
        %
        %    ResCalc(k) =  r_VecCalc(k) ;%+ SVPos(k,4) ;
    end
    
    %ResCalc(Nclock) = EstimatedClockBiasCoast;  %add receiver clock term        


    %Design Matrix
    %generate elements of M matrix         
          
        
        for k = 1:N
        M(k,1) =  -(SVPos(k,1) - TempBase_B_Pos(1))/r_VecCalc(k);
        M(k,2) =  -(SVPos(k,2) - TempBase_B_Pos(2))/r_VecCalc(k);
        M(k,3) =  -(SVPos(k,3) - TempBase_B_Pos(3))/r_VecCalc(k);
        %M(k,4) = 1;
         end
    
%     %add the clock 
%     M(Nclock,1) = 0;
%       M(Nclock,2) = 0;
%         M(Nclock,3) = 0;
%           %M(Nclock,4) = 1;
          
    end%      if ClockAugmentation == 1
        
        

      if ClockAugmentation == 0

    %residual vector
       ResVec = PRMeasured' - ResCalc';
       
       %scale Delta_T see p228 of Lee
       
%        ResVec(Nclock) = ResVec(Nclock)/ScaleFactor;
    
    
      end
      
        
       if ClockAugmentation == 1

    %residual vector
        ResVec = PRMeasured' - ResCalc';
    
    
        end

    %LSQ solution

   
    
    A = M'*M;
    b = M'*ResVec;

    Delta_X = inv(A)*b;
    %Delta_X = A\b;

    % limit = sqrt(Delta_X(1)^2 +  Delta_X(2)^2  + Delta_X(3)^2 +Delta_X(4)^2);

    limit = sqrt(Delta_X(1)^2 +  Delta_X(2)^2  + Delta_X(3)^2 );
    numlimitAnt1_2(j) = limit;
    if limit < 0.0000005;
        numiterAnt1_2 = j;
        break;
    else
        numiterAnt1_2 = 200;
    end

    if j == 200
        LSQ_Fail = 1;  %Solution Cannot converge
    end
    
    
    if ClockAugmentation == 0;    
        

    TempBase_B_Pos(1) = TempBase_B_Pos(1) + Delta_X(1);
    TempBase_B_Pos(2) = TempBase_B_Pos(2) + Delta_X(2);
    TempBase_B_Pos(3) = TempBase_B_Pos(3) + Delta_X(3);
    TempBase_B_Pos(4) = TempBase_B_Pos(4) + Delta_X(4);    
    

    end;
    
   if ClockAugmentation == 1;
       
    TempBase_B_Pos(1) = TempBase_B_Pos(1) + Delta_X(1);
    TempBase_B_Pos(2) = TempBase_B_Pos(2) + Delta_X(2);
    TempBase_B_Pos(3) = TempBase_B_Pos(3) + Delta_X(3);   
    TempBase_B_Pos(4) = EstimatedClockBiasCoast;
    
   end
       
    

end  %for j = ..



%final solution
Xb = TempBase_B_Pos(1);
Yb = TempBase_B_Pos(2);
Zb = TempBase_B_Pos(3);
dtb = TempBase_B_Pos(4);


if ClockAugmentation == 0

%Calculate DOPvalues



%Transform M matrix into LTP coordinates (NEU) (for HDOP/VDOP calculation) 
            
    

Position = [Xb, Yb, Zb];


[Latitude,Longitude,Height] = ECEF2LLH(Position);


Telev2 = T_ECEF2LTP(Longitude,Latitude);

%make Telev2 into a 4 by 4

Telev2(4,1:3) = [0 0 0];
Telev2(1:4,4) = [0 0 0 1];
           
    
    
H_LTP = M*Telev2';   %note this is transposed 
 
 

AA = (H_LTP'*H_LTP)^-1;



var_x = AA(1,1);
var_y = AA(2,2);
var_z = AA(3,3);
var_dt = AA(4,4);

GDOP = sqrt(var_x + var_y + var_z + var_dt);

PDOP = sqrt(var_x + var_y + var_z);

HDOP = sqrt(var_x + var_y);

VDOP = sqrt(var_z);

TDOP = sqrt(var_dt); %IS THIS RIGHT? CHECK THAT THIS IS HOW TDOP IS CALCULATED

end  %if ClockAugmentation == 0


 
if ClockAugmentation == 1  %if clock augmentation is turned on
    
   
    
    %Calculate DOPvalues

    

%Transform M matrix into LTP coordinates (NEU) (for HDOP/VDOP calculation) 
            
    

Position = [Xb, Yb, Zb];


[Latitude,Longitude,Height] = ECEF2LLH(Position);


Telev2 = T_ECEF2LTP(Longitude,Latitude);
   

%Transform A matrix into LTP coordinates (for HDOP/VDOP calculation)          
    
    
H_LTP = M*Telev2';   %equivalent to H in the paper.

 

G = (H_LTP'*H_LTP)^-1; %equivalent to G in the paper

%k_squared = 1 %Have to work this value out with the clock drift parameters of variance etc.


GG = G*(eye(3) + k_squared*H_LTP'*ones(N)*H_LTP*G);


var_x = GG(1,1);
var_y = GG(2,2);
var_z = GG(3,3);
var_dt = 1111111111111111111111; %have to give this a proper value yet.

    
HDOP = sqrt(var_x + var_y);

VDOP = sqrt(var_z);

GDOP = sqrt(var_x + var_y + var_z + var_dt);

PDOP = sqrt(var_x + var_y + var_z);


TDOP = sqrt(var_dt);



 end %if ClockAugmentation == 1  


%Outputs
SolutionVec = [Xb,Yb,Zb,dtb];
VarSolutionVec = [var_x,var_y,var_z,var_dt];
NumIterations = numiterAnt1_2;
ResidualVector = ResVec ; %Last iteration residual vector, output for RAIM using Least squares residual method
DOP = [GDOP PDOP HDOP VDOP TDOP];