
function [SolutionVec, VarSolutionVec, NumIterations, ResidualVector, M, LSQ_Fail, limit] = GARD_LSQVel(UserPos,UserVel,N,PRMeasured,SVPos, SVVel);

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

%start least squares

%initialise variables
LSQ_Fail = 0;

%Initial user guess



Xb = UserPos(1);
Yb = UserPos(2);
Zb = UserPos(3);
dtb = UserPos(4);

XbVel = UserVel(1);
YbVel = UserVel(2);
ZbVel = UserVel(3);
dtbVel = UserVel(4);


TempBase_B_Vel = [XbVel,YbVel,ZbVel, dtbVel];

TempBase_B_Pos = [Xb,Yb,Zb, dtb];



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

        %predicted relative velocity of sv and receiver

        r_VecCalcVel(k) = (SVVel(k,1) - TempBase_B_Vel(1))*(SVPos(k,1)-TempBase_B_Pos(1)) + (SVVel(k,2) - TempBase_B_Vel(2))*(SVPos(k,2)-TempBase_B_Pos(2)) + (SVVel(k,3) - TempBase_B_Vel(3))*(SVPos(k,3)-TempBase_B_Pos(3));


        Relative_Velocity(k) = r_VecCalcVel(k)/r_VecCalc(k);



    end




    for k = 1:N


        ResCalcVel(k) =  Relative_Velocity(k) + TempBase_B_Vel(4) - SVVel(k,4) ;


    end



    %Design Matrix
    %generate elements of M matrix

    for k = 1:N


        M(k,1) =  -(SVPos(k,1) - TempBase_B_Pos(1))/r_VecCalc(k);
        M(k,2) =  -(SVPos(k,2) - TempBase_B_Pos(2))/r_VecCalc(k);
        M(k,3) =  -(SVPos(k,3) - TempBase_B_Pos(3))/r_VecCalc(k);
        M(k,4) = 1;


    end


    %residual vector
    ResVec = PRMeasured' - ResCalcVel';



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
    TempBase_B_Vel(4) = TempBase_B_Vel(4) + Delta_X(4);





end

%final solution
XbVel = TempBase_B_Vel(1);
YbVel = TempBase_B_Vel(2);
ZbVel = TempBase_B_Vel(3);
dtbVel = TempBase_B_Vel(4);



%DOPvalues
var_x = A(1,1);
var_y = A(2,2);
var_z = A(3,3);
var_dt = A(4,4);

%Outputs
SolutionVec = [XbVel,YbVel,ZbVel,dtbVel];
VarSolutionVec = [var_x,var_y,var_z,var_dt];
NumIterations = numiterAnt1_2;
ResidualVector = ResVec ; %Last iteration residual vector, output for RAIM using Least squares residual method





