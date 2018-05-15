

for i = startepochHighRate:endepochHighRate


    %=====================================================================
    %Calculated Wind from Truth data
    %=====================================================================

    %Convert ECEF velocities into NED

    %if i = startepochHighRate

    Position = [Xpos_truth(i-1), Ypos_truth(i-1), Zpos_truth(i-1)];

    %transform
    [Latitude,Longitude,Height] = ECEF2LLH(Position);

    T_ECEF2NEDConvert = T_ECEF2NED(Longitude,Latitude);
    VelNED(i,1:3) = T_ECEF2NEDConvert*[Xvel_truth(i-1),Yvel_truth(i-1),Zvel_truth(i-1)]';


    %calculate Vtrue airspeed

    %calculate true airspeed
    T_sealevel = 288.15; %K
    Tstatic =  AtmosTruth(2,i);
    V_TAS(i) = AtmosTruth(4,i)*Mach_truth(i)*sqrt(Tstatic/T_sealevel);

    V_TAS(i) = Airspeed_truth(i);

    Vnorthtemp(i) =  V_TAS(i)*cos(Pitch_truth(i)-Alpha_truth(i))*cos(Yaw_truth(i)+Beta_truth(i));
    Veasttemp(i) =  V_TAS(i)*cos(Pitch_truth(i)-Alpha_truth(i))*sin(Yaw_truth(i)+Beta_truth(i));

    %not sure if this is the correct formula but it works
    Vdowntemp(i) =  V_TAS(i)*sin(Alpha_truth(i)-Pitch_truth(i));  %Vup = V_TAS(i)*sin(Pitch_truth(i)-Alpha_truth(i));
    

    VWind_NorthCalc(i) =  VelNED(i,1) - Vnorthtemp(i);
    VWind_EastCalc(i) =  VelNED(i,2) - Veasttemp(i);
    VWind_DownCalc(i) =  VelNED(i,3) - Vdowntemp(i);


    %Transform the True wind into NED

    TMatrixBodytoNED = T_Body2NED(Roll_truth(i),Pitch_truth(i), Yaw_truth(i));
    VWind_Truth = TMatrixBodytoNED*[u_wind_truth(i), v_wind_truth(i), w_wind_truth(i)]';

    VWind_North_Truth(i) = VWind_Truth(1);
    VWind_East_Truth(i) = VWind_Truth(2);
    VWind_Down_Truth(i) =  VWind_Truth(3);


    %calculate the error in the estimate and INS

    VWind_North_Error(i) = VWind_North_Truth(i) - VWind_NorthCalc(i);
    VWind_East_Error(i) = VWind_East_Truth(i) - VWind_EastCalc(i);
    VWind_Down_Error(i) = VWind_Down_Truth(i) - VWind_DownCalc(i);



end





%Estimate Wind at 1 Hz, from aircraft sensors

for i = startepoch:endepoch
    
    %Convert ECEF velocities into NED

    %if i = startepochHighRate

    Position = [Xpos_truth(i-1), Ypos_truth(i-1), Zpos_truth(i-1)];

    %transform
    [Latitude,Longitude,Height] = ECEF2LLH(Position);

    T_ECEF2NEDConvert = T_ECEF2NED(Longitude,Latitude);
    VelNED(i,1:3) = T_ECEF2NEDConvert*[Xvel_truth(i-1),Yvel_truth(i-1),Zvel_truth(i-1)]';


    %calculate Vtrue airspeed

    %calculate true airspeed
    T_sealevel = 288.15; %K
    Tstatic =  AtmosTruth(2,i);
    V_TAS(i) = AtmosTruth(4,i)*Mach_truth(i)*sqrt(Tstatic/T_sealevel);

    V_TAS(i) = Airspeed_truth(i);

    Vnorthtemp(i) =  V_TAS(i)*cos(Pitch_truth(i)-Alpha_truth(i))*cos(Yaw_truth(i)+Beta_truth(i));
    Veasttemp(i) =  V_TAS(i)*cos(Pitch_truth(i)-Alpha_truth(i))*sin(Yaw_truth(i)+Beta_truth(i));

    %not sure if this is the correct formula but it works
    Vdowntemp(i) =  V_TAS(i)*sin(Alpha_truth(i)-Pitch_truth(i));  %Vup = V_TAS(i)*sin(Pitch_truth(i)-Alpha_truth(i));
    

    VWind_NorthCalc(i) =  VelNED(i,1) - Vnorthtemp(i);
    VWind_EastCalc(i) =  VelNED(i,2) - Veasttemp(i);
    VWind_DownCalc(i) =  VelNED(i,3) - Vdowntemp(i);



end




% [CurveFitNorth, BestFitNorth] = LeastSquaresBestFit(VWind_NorthCalcGPSINS, 1, 1);
%
% [CurveFitEast, BestFitNorth] = LeastSquaresBestFit(VWind_EastCalcGPSINS, 1, 1);
%
% [CurveFitDown, BestFitNorth] = LeastSquaresBestFit(VWind_DownCalcGPSINS, 1, 1);

%average out the wind estimate using least squares fitting to get rid of
%the high frequency noise on the wind estimate.


