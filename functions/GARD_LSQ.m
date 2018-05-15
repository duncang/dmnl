function [SolutionVec, VarSolutionVec, NumIterations, ResidualVector, M, LSQ_Fail, limit, DOP] = GARD_LSQ(UserPos,N,PRMeasured,SVPos,W,DGPS_flag)
%GPSPosLSQ VarSolutionVec_Observed NumIterations_Observed_Vel ResVec_Observed M_Observed LSQ_Fail_Observed limit_Observed DOP_Observed(1:5)] = GARD_LSQ(PosTruth,N,PRMeasured_Simulated,SatPos);
%
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
%put the DOPS in a vector, the order is [GDOP PDOP HDOP VDOP TDOP]
%
% $Id: GARD_LSQ.m 3546 2010-06-17 13:26:36Z greerd $
%

global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;

%start least squares

%initialise variables
LSQ_Fail = 0;

% check if we have a weights matrix
if ~exist('W','var')
    W = eye(N);
end


% check if we have differential corrected PRs
if ~exist('DGPS_flag','var')
    DGPS_flag = 0;
end

%Initial user guess


Xb = UserPos(1);   
Yb = UserPos(2);
Zb = UserPos(3);
dtb = UserPos(4);


TempBase_B_Pos = [Xb,Yb,Zb, dtb];


% %add orbital error (to model ephemeris errors)
% for k = 1:N
%     for m = 1:3
%     SVPos(k,m) = SVPos(k,m)+ abs(10*rand(1));  
%     end
% end



%Start loop here, loop for LSQ
for j = 1:200   %Least square iterations counter

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

        if DGPS_flag == 0
            ResCalc(k) =  r_VecCalc(k) + TempBase_B_Pos(4) - SVPos(k,4) - delta_pr_omegaedot(k);
        else
            ResCalc(k) =  r_VecCalc(k) + TempBase_B_Pos(4) - delta_pr_omegaedot(k);
        end
        

%just while im not adding corrections to the pseudoranges (for debugging gps/ins):

        %ResCalc(k) =  r_VecCalc(k) ;

  
        
        
    end



    %Design Matrix
    %generate elements of M matrix

    for k = 1:N
        M(k,1) =  -(SVPos(k,1) - TempBase_B_Pos(1))/r_VecCalc(k);
        M(k,2) =  -(SVPos(k,2) - TempBase_B_Pos(2))/r_VecCalc(k);
        M(k,3) =  -(SVPos(k,3) - TempBase_B_Pos(3))/r_VecCalc(k);
        M(k,4) = 1;
    end
    
%     
%       for k = 1:N
%         M(k,1) =  -(SVPos(k,1) - TempBase_B_Pos(1))/ResCalc(k);
%         M(k,2) =  -(SVPos(k,2) - TempBase_B_Pos(2))/ResCalc(k);
%         M(k,3) =  -(SVPos(k,3) - TempBase_B_Pos(3))/ResCalc(k);
%         M(k,4) = 1;
%     end


    %residual vector
    ResVec = PRMeasured' - ResCalc';
      

    %LSQ solution

    A = M'*W*M;
    b = M'*W*ResVec;

    Delta_X = inv(A)*b;
    %Delta_X = A\b;

    % limit = sqrt(Delta_X(1)^2 +  Delta_X(2)^2  + Delta_X(3)^2 +Delta_X(4)^2);

    limit = sqrt(Delta_X(1)^2 +  Delta_X(2)^2  + Delta_X(3)^2 );
    numlimitAnt1_2(j) = limit;
    %if limit < 0.0000005;
    if limit < 1e-8;
        numiterAnt1_2 = j;
        break;
    else
        numiterAnt1_2 = 200;
    end

    if j == 200
        LSQ_Fail = 1;  %Solution Cannot converge
    end

    TempBase_B_Pos(1) = TempBase_B_Pos(1) + Delta_X(1);
    TempBase_B_Pos(2) = TempBase_B_Pos(2) + Delta_X(2);
    TempBase_B_Pos(3) = TempBase_B_Pos(3) + Delta_X(3);
    TempBase_B_Pos(4) = TempBase_B_Pos(4) + Delta_X(4);


end

%final solution
Xb = TempBase_B_Pos(1);
Yb = TempBase_B_Pos(2);
Zb = TempBase_B_Pos(3);
dtb = TempBase_B_Pos(4);

% %Calculate DOPvalues
% %Use the true position to calculate these values
% X = UserPos(1);
% Y = UserPos(2);
% Z = UserPos(3);

%     %Calculate estimated lat , lon , hgt of user. 
%     el = atan2(Zb,sqrt(Xb^2 + Yb^2)); %lat, lon, hgt in ECEF
%     az = atan2(Yb,Xb);
%     hgt = sqrt(Xb^2 + Yb^2 +Zb^2);
% 
%     %direction cosines (page 202 of Kathon this is the transformation
%     %matrix for from ECEF to LTP. (East , North, Up)
% 
%     sel = sin(el);
%     cel = cos(el);
%     saz = sin(az);
%     caz = cos(az);
%     Telev(1,1) = -sel*caz;
%     Telev(1,2) = -sel*saz;
%     Telev(1,3) = cel;
%     Telev(1,4) = 0;
%     Telev(2,1) = -saz;
%     Telev(2,2) = caz;
%     Telev(2,3) = 0.0;
%     Telev(2,4) = 0;
%     Telev(3,1) = cel*caz;
%     Telev(3,2) = cel*saz;
%     Telev(3,3) = sel;
%     Telev(3,4) = 0;
%     Telev(4,1) = 0;
%     Telev(4,2) = 0;
%     Telev(4,3) = 0;
%     Telev(4,4) = 1;


%Transform M matrix into LTP coordinates (for HDOP/VDOP calculation)
        
Position = [Xb, Yb, Zb];


[Latitude,Longitude,Height] = ECEF2LLH(Position);



Telev2 = T_ECEF2LTP(Longitude,Latitude);

%make Telev2 into a 4 by 4

Telev2(4,1:3) = [0 0 0];
Telev2(1:4,4) = [0 0 0 1];
           
    
    
H_LTP = M*Telev2';   

%I have verified by comparing H_LTP with Hltp calculated from teh azimuth
%and elevation angles that you should multiply M by the transpose of
%Telev2. (see page 201 in Kayton).

%however the derivatino done by aaron dando and the book vallado shows that
%the transformation given on p.202 of kayton is not the transpose as is
%indicated. So either Kayton is wrong or the VDOPs are supposed to be less
%than the HDOPs in this particular case (which i think is unlikely).


%see GPS:Theory and practice by Hofmann-Wellenhof for a good explanation of
%DOP and the correct coordinate transformation.


% Hltp =  [ -0.985835558498596        -0.133366156149999        -0.101694247593142                         1;
%          0.038197909035887         -0.13434294997918        -0.990198410186654                         1;
%           0.74659555920773        -0.625180917858734        -0.227472835557583                         1;
%          0.329551552260845         0.785885557944164        -0.523239585866033                         1;
%          0.655049537492162         0.714890491946064        -0.244625607728257                         1;
%         -0.389261901831207         0.413323943023346        -0.823188003985958                         1;
%          0.528917533669639         0.292822631011521        -0.796555804287632                         1;
%         -0.697212676384312         0.650473340186093        -0.301295399228351                         1;
%       -0.00520300659254806        -0.817302778193395        -0.576184950757789                         1;
%          0.202850707251317        -0.822362184238041        -0.531575045033982                         1;
%         -0.853387610017255         0.245204086797665         -0.46000493789389                         1;]




%H_LTP = Telev2.*M;  

AA = (H_LTP'*H_LTP)^-1;

var_x = AA(1,1);
var_y = AA(2,2);
var_z = AA(3,3);
var_dt = AA(4,4);

GDOP = sqrt(var_x + var_y + var_z + var_dt);

PDOP = sqrt(var_x + var_y + var_z);

HDOP = sqrt(var_x + var_y);

VDOP = sqrt(var_z);

TDOP = sqrt(var_dt);




%Outputs
SolutionVec = [Xb,Yb,Zb,dtb];
VarSolutionVec = [var_x,var_y,var_z,var_dt];
NumIterations = numiterAnt1_2;
ResidualVector = ResVec ; %Last iteration residual vector, output for RAIM using Least squares residual method
DOP = [GDOP PDOP HDOP VDOP TDOP];