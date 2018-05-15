function [GPSPosLSQ, GPSVelLSQ, SatPos,SatVel,SV_AboveElevationMask,PRMeasured_Simulated,PRRateMeasured_Simulated,N, Sat_PRN_Vec,Elevation,Azimuth,ElevationBody,AzimuthBody,ResVec,ResVec_Observed_Vel,limit_Observed,limit_Observed_Vel,M,NoSats,DOP_Observed,VarSolutionVec_Observed,AA_out,PR_SM_prev] = SimulatedGPSSub(PosTruth, VelTruth,PositionPrevious, VelocityPrevious, GPSWeek, GPSSec,ElevationMask,NavData,PHI,THETA,PSI,GPSPr_noise,GPSPr_Rate_noise,PRMeasPredict_Previous,i,MaxSlopeSatelliteT);

%$Id: SimulatedGPSSub.m 2531 2009-04-28 04:28:16Z bruggema $
%Troy Bruggemann  11 August 2005
GPSConstants; %All generic constants put in this script, as global variables;
%turn JIT accelerator on (speeds up code processing)
SensorNoiseParameters;

feature accel on %only works in matlab 7


%VelTruth is in ECEF x, y, z

%set these to size 32 , otherwise if the highest satellite is only 21 (22
%drops out) then the size of this array is only 21 and it conflicts with
%previous sizes of 22, so set it fixed to 32
SatPos(1:32,1:4) = zeros(32,4);
SatVel(1:32,1:4) = zeros(32,4);
%this is the GPS processor which calculates simulated pseudoranges and
%satellite positions for the  GPSandModelEKF.m


%------------------------------------------------------------------------
%Atmospheric Models
%-----------------------------------------------------------------------

%IONOSPHERIC PARAMETERS for IONO Model
ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA
BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA


%--------------------------
%Optimised GPS Constellation - DO229C WAAS MOPS
%---------------------------

n = 1;

Sat_PRN_Vec = zeros(1,32);  %initialise the vector with zeros. Columns of this vector is the PRN number of the satellite

for SV = 1:32

    [SV_X_Data(SV) SV_Y_Data(SV) SV_Z_Data(SV) SV_T_Data(SV) ValidData_Satellite(SV)] = GPSOrbitPropagatorOptimal(GPSWeek, GPSSec, SV, NavData);

    [SV_X_Data1(SV) SV_Y_Data1(SV) SV_Z_Data1(SV) SV_T_Data1(SV) SV_Xvel_Data(SV) SV_Yvel_Data(SV) SV_Zvel_Data(SV) SV_Tvel_Data(SV) SV_Xacc_Data(SV) SV_Yacc_Data(SV) SV_Zacc_Data(SV) SV_Tacc_Data(SV) ValidData_Satellite(SV)] = GPSOrbitPropagatorOptimalVelocities(GPSWeek, GPSSec, SV, NavData);

    
    
    
    %[SV_X_Data1(SV) SV_Y_Data1(SV) SV_Z_Data1(SV) SV_T_Data1(SV) SV_Xvel_Data(SV) SV_Yvel_Data(SV) SV_Zvel_Data(SV) SV_Tvel_Data(SV) SV_Xacc_Data(SV) SV_Yacc_Data(SV) SV_Zacc_Data(SV) SV_Tacc_Data(SV) ValidDataSatVels(SV)] = GPSOrbitPropagatorOptimalVelocitiesFixed(GPSWeek, GPSSec, SV, NavData);



    %             if ValidData_Satellite(SV) == 1 && ValidDataSatVels(SV) == 1
    %
    %
    %
    %
    %             end

    %calculate satellite elevations

    [Azimuth(SV), Elevation(SV)] = AzEl(PosTruth(1:3), [SV_X_Data(SV), SV_Y_Data(SV), SV_Z_Data(SV)]);
    Azimuth(SV) = rad2deg(Azimuth(SV));
    Elevation(SV) = rad2deg(Elevation(SV));


    [AzimuthBody(SV), ElevationBody(SV)] = AzElAntenna(PosTruth(1:3), [SV_X_Data(SV), SV_Y_Data(SV), SV_Z_Data(SV)],PHI,THETA,PSI);

    AzimuthBody(SV) = rad2deg(AzimuthBody(SV));
    ElevationBody(SV)  = rad2deg(ElevationBody(SV));
    
    
    
    %Calculate elevation mask which the earth creates
    
    
    [latEl,lonEl,hgtEl] = ECEF2LLH(PosTruth);
    
    [RnMa, ReMa] = WGS84_calcRnRe(latEl);     
    RnhMa = RnMa + hgtEl;
    RehMa = ReMa + hgtEl;
    
    RoMa = sqrt(RnMa*ReMa);
   
    
    ElevationMaskEarth = -(pi/2 - asin(RoMa/(RoMa+hgtEl))); %+ 5*pi/180;   %add 5 degrees to account for approximation of earth surface. %the negative is there because its a negative angle wrt LTP
    
        
    
         

    

    %if (Elevation(SV) > ElevationMask)  && (ValidData_Satellite(SV) == 1 && ValidDataSatVels(SV) == 1)
   % if (ElevationBody(SV) > ElevationMask)  && (ValidData_Satellite(SV) == 1)% && ValidDataSatVels(SV) == 1)
        
    if (Elevation(SV) > ElevationMaskEarth) && (ElevationBody(SV) > ElevationMask)  && (ValidData_Satellite(SV) == 1)% && ValidDataSatVels(SV) == 1)

        SV_AboveElevationMask(n) = SV ;  %PRN numbers of the satellites above the elevation mask

        SatPos(SV,1) = SV_X_Data(SV);   %adding error to simulate orbit determination error.
        
        SatPos(SV,2) = SV_Y_Data(SV);
        SatPos(SV,3) = SV_Z_Data(SV);
        SatPos(SV,4) = c*SV_T_Data(SV); %convert to metres

        SatVel(SV,1) = SV_Xvel_Data(SV);
        SatVel(SV,2) = SV_Yvel_Data(SV);
        SatVel(SV,3) = SV_Zvel_Data(SV);
        SatVel(SV,4) = c*SV_Tvel_Data(SV); %convert to metres

        Sat_PRN_Vec(SV) = 1;



        %set ionodelay to 0 for now
        %IonoDelay = 0;%10*randn(1);  %this gives the gps pr a random walking/markov thing        
                   
     
         
      % IonoDelay =  abs(GPSPr_noise*rand(1)*3);  %need this rand(1) here so that the ionodelay is differnet for each PR
         
      
       GPSPr_noiseError = GPSPr_noise(SV,i);
          
         GPSPr_Rate_noiseError = GPSPr_Rate_noise(SV,i);
         
         
          
          
        %  GPSPr_noiseError = 0;         
        % GPSPr_Rate_noiseError = 0;
          
          
        % IonoDelay = 0;
              [PRMeasured_Simulated(SV),PRRateMeasured_Simulated(SV),GeometricRange(SV)] = GARDSim_PRPredictSimulated([SV_X_Data(SV) SV_Y_Data(SV) SV_Z_Data(SV) SV_T_Data(SV)] ,[SV_Xvel_Data(SV) SV_Yvel_Data(SV) SV_Zvel_Data(SV) SV_Tvel_Data(SV)] ,PosTruth,VelTruth,GPSPr_noiseError,GPSPr_Rate_noiseError,0);
       
        %[PRMeasured_Simulated(SV),PRRateMeasured_Simulated(SV),GeometricRange(SV)] = GARDSim_PRPredictSimulated([SV_X_Data(SV) SV_Y_Data(SV) SV_Z_Data(SV) SV_T_Data(SV)] ,[SV_Xvel_Data(SV) SV_Yvel_Data(SV) SV_Zvel_Data(SV) SV_Tvel_Data(SV)] ,PosTruth,VelTruth,GPSPr_noiseError,GPSPr_Rate_noiseError,PRMeasPredict_Previous(SV));
        n = n+1;



% 
% %         =================================================
% %         Adding errors onto PR
% %         =================================================
%           if SV == 1 && (GPSSec >= 432020 && GPSSec <= 432022) %20 seconds after the start time.
%  
%              PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 30;  %this is on SV 1
% % 
% % 
%           end
% %  
%            if SV == SVwithError1 && (GPSSec >= TimeLow && GPSSec <= TimeHigh) %20 seconds after the start time.
% 
%             PRMeasured_Simulated(SV) = round( PRMeasured_Simulated(SV) + BiasError);  %this is on SV 1
% 
%               %PRMeasured_Simulated(SV) = round(PRMeasured_Simulated(SV) + BiasError);  %round this so i can tell which one i added the error too.
%               
%            end








% 

 % if (i >= 50 && i <= 52 && SV == SVwithError1 )
    
   % if (i == 30  && SV == SVwithError1 )
%if (i >=50 && i<=52 && (SV == SVwithError1 || SV == SVwithError2 ))
% %        
%  % if (i == 50 && (SV == SVwithError1 || SV == SVwithError2 ))
% 
          %PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 50;  %this is on SV 1
% 
       %   PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 5*(i-30);
%         
%         %random constant error:
%                 
%         %random walk error
%         %booboo = 1;

%  %PRMeasured_Simulated(SV) =  PRMeasured_Simulated(SV) + BiasError;
% 
% % %         
% %         





%ERRORS ADDED FOR THESIS RESULTS
%       %0.5 m/s is the smallest fault i can show is detected. 
%       
 %if (i >= 30 && SV == 1 )  %SV1 is most difficult to detect sat (max slope) i think
 

 %just pick whatever satellite is max slope when i = 3. 
 

 
     
 if (i >= 30 && SV == MaxSlopeSatelliteT )  %SV1 is most difficult to detect sat (max slope) i think
       
   
     
    %  PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 0.5*(i-30);
     
      %PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 2.5*(i-30);
      
     % PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 5.0*(i-30);
     
     
     
    % PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 1*(i-30);
     
     %PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 2.5*(i-30);
    
   
   
    %PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 0.5*(i-30);
      
  %  PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 50;
    
    
 %PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 1*(i-30);

 
 %step fault
 % PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 300;

   end
    
    









% 
% % 
% 
%   if (i >= 50 && i <= 52 && SV == SVwithError1 )
%     
%    % if (i == 30  && SV == SVwithError1 )
% %if (i >=50 && i<=52 && (SV == SVwithError1 || SV == SVwithError2 ))
% % %        
% %  % if (i == 50 && (SV == SVwithError1 || SV == SVwithError2 ))
% % 
%           %PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 50;  %this is on SV 1
% % 
%           PRMeasured_Simulated(SV) = PRMeasured_Simulated(SV) + 5*(i-30);
% %         
% %         %random constant error:
% %                 
% %         %random walk error
% %         %booboo = 1;
% 
% %  %PRMeasured_Simulated(SV) =  PRMeasured_Simulated(SV) + BiasError;
% % 
% % % %         
% % %         
%    end
%     
%     









        
    end

 end


%number of satellites above elevation mask
%N = length(SV_AboveElevationMask);


if exist('SV_AboveElevationMask') ~= 0
    %number of satellites above elevation mask
N = length(SV_AboveElevationMask);
else
    N = 0;
    
end

    

%N = 4; %always keep it at 9 satellites for now


%         %generate PR's for each satellite
%
%         for SV = 1:32
%
%
%             if Sat_PRN_Vec(SV) == 1
%
%                          IonoDelay = 0;
%                    [PRMeasured_SimulatedTemp(SV),GeometricRangeTemp(SV)] = GARDSim_PRPredictSimulated(SatPos(SV,:),PosTruth,IonoDelay);
%
%             end
%
%
%         end
%

%pick out the satellites we want   for the Least squares



k = 1;
for SV = 1:32

    if Sat_PRN_Vec(SV) == 1

        PRMeasured_SimulatedLSQ(k) =  PRMeasured_Simulated(SV);

        
          PRRateMeasured_SimulatedLSQ(k) = PRRateMeasured_Simulated(SV);
 
        SatPosLSQ(k,1) = SV_X_Data(SV); %+  SatOrbError(SV,1,i);    %5.7 is orbital error added.
        SatPosLSQ(k,2) = SV_Y_Data(SV);% +  SatOrbError(SV,2,i);
        SatPosLSQ(k,3) = SV_Z_Data(SV); %+  SatOrbError(SV,3,i);
        SatPosLSQ(k,4) = c*SV_T_Data(SV); %convert to metres
        
        
        
        SatVelLSQ(k,1) = SV_Xvel_Data(SV); %+  SatOrbError(SV,1,i);    %5.7 is orbital error added.
        SatVelLSQ(k,2) = SV_Yvel_Data(SV);% +  SatOrbError(SV,2,i);
        SatVelLSQ(k,3) = SV_Zvel_Data(SV); %+  SatOrbError(SV,3,i);
        SatVelLSQ(k,4) = c*SV_Tvel_Data(SV); %convert to metres

       
        
        
        
        
        %smooth the PR's, 
        
       % PRMeasured_Simulated(SV),PRRateMeasured_Simulated(SV)
        
        %p 40 of 229D
        
%         
%         if PRMeasPredict_Previous(SV) == 0  % if its zero, there has been no past measurement, eg start epoch
%             
%             
%             %PRMeasPredict_Previous(SV) = PRMeasured_Simulated(SV);     %if zero, initialise with the current PR value
%             
%          PRSmoothed(SV) = PRMeasured_Simulated(SV) ;
%         
%           PRMeasured_SimulatedLSQ_SM(k) =  PRSmoothed(SV);
%          
%            PR_SM_prev(SV) = PRSmoothed(SV);
%         else
%              
%                    
%         
%         Pproj = PRMeasPredict_Previous(SV) + (PRRateMeasured_Simulated(SV)-GPSPr_Rate_noise(SV,i));  % PRMeasPredict_Previous is the previous smoothed CP. I subtract the rate noise to test
%         %it , the smoothed least squares solution is worse than not using
%         %it...
%         
%         alpha_sm = 1/100;             
%         
%         PRSmoothed(SV) = alpha_sm*PRMeasured_Simulated(SV) + (1-alpha_sm)*Pproj;
%         
%          
%        PRMeasured_SimulatedLSQ_SM(k) =  PRSmoothed(SV);
%        
%        PR_SM_prev(SV) = PRSmoothed(SV);
%        
%         end
%        
%        
%        
        k = k+1;
       
               
    else   %else if no measurements, sent the SV index to 0
        
        
      PR_SM_prev(SV) = 0;  
        
        
    end


%I've commented this smoothing stuff out for now, since i havent accounted
%for N < 4 satellitse and i dont really use the smoothed PR's

PR_SM_prev = 0;
PRSmoothed= 0;
PRMeasured_SimulatedLSQ_SM = 0;
PRMeasPredict_Previous = 0;
    
    
if exist('PRMeasured_SimulatedLSQ') == 0
    %number of satellites above elevation mask
NoSats = 1; 
else
    NoSats = 0; 
    
end   
        
        
        
end




if NoSats == 1 %|| N < 4
    
    GPSPosLSQ = [0,0,0,0];
    GPSVelLSQ = [0,0,0,0];  
    SatPos = 0;
    SatVel= 0;
    
    SV_AboveElevationMask= 0;
    PRMeasured_Simulated= 0;
        PRRateMeasured_Simulated= 0;
        Sat_PRN_Vec= 0;
   
    
Elevation= 0;
Azimuth= 0;
AzimuthBody= 0;
ResVec= 0;
ElevationBody= 0;
limit_Observed = 0; 
 limit_Observed_Vel = 0;

 ResVec_Observed_Vel =0;
   LSQ_Fail_Observed =0;
    M_Observed =0;
    ResVec_Observed =0;
    NumIterations_Observed =0;
    VarSolutionVec_Observed =0;
    
    DOP_Observed = 0;
    AA_out= 0;
    M = 0;
    
    
    
elseif  N < 4 && NoSats ~= 1
    GPSPosLSQ = [0,0,0,0];
    GPSVelLSQ = [0,0,0,0];  
   % SatPos = 0;
   % SatVel= 0;
    
   % SV_AboveElevationMask= 0;
   % PRMeasured_Simulated= 0;
    %    PRRateMeasured_Simulated= 0;
      %  Sat_PRN_Vec= 0;
   
    
%Elevation= 0;
%Azimuth= 0;
%AzimuthBody= 0;

%ElevationBody= 0;

ResVec= 0;
limit_Observed = 0; 
 limit_Observed_Vel = 0;

 ResVec_Observed_Vel =0;
   LSQ_Fail_Observed =0;
    M_Observed =0;
    ResVec_Observed =0;
    NumIterations_Observed =0;
    VarSolutionVec_Observed =0;
    
    M = 0;
    
NoSats  = 0;


DOP_Observed  = 0;

VarSolutionVec_Observed = 0;

AA_out= 0;

PR_SM_prev= 0;

%     ,,,,,DOP_Observed,VarSolutionVec_Observed,AA_out,PR_SM_prev] = SimulatedGPSSub(PosTruth, VelTruth,PositionPrevious, VelocityPrevious, GPSWeek, GPSSec,ElevationMask,NavData,PHI,THETA,PSI,GPSPr_noise,GPSPr_Rate_noise,PRMeasPredict_Previous,i,MaxSlopeSatelliteT);
% 
%     
%     
%     
%     
    
    
    
    
    
    
else
    
    
    
    




[GPSPosLSQ VarSolutionVec_Observed NumIterations_Observed ResVec_Observed M_Observed LSQ_Fail_Observed limit_Observed DOP_Observed(1:5),AA_out] = GARD_LSQ_NSS(PositionPrevious,N,PRMeasured_SimulatedLSQ,SatPosLSQ);


%using smoothed PR's for measurements
%[GPSPosLSQ VarSolutionVec_Observed NumIterations_Observed ResVec_Observed M_Observed LSQ_Fail_Observed limit_Observed DOP_Observed(1:5),AA_out] = GARD_LSQ_NSS(PositionPrevious,N,PRMeasured_SimulatedLSQ_SM,SatPosLSQ);






ResVec = ResVec_Observed;
M = M_Observed;

%UserPos = ([PosTruth(1,i),PosTruth(1,i),PosTruth(1,i),c*dTposTCXO(i)]); %this will be a copy of the values set by the user in the flight planning section


%[SolutionVec, VarSolutionVec, NumIterations, ResidualVector, M, LSQ_Fail, limit] = GARD_LSQVel(UserPos,UserVel,N,PRMeasured,SVPos, SVVel);


[GPSVelLSQ VarSolutionVec_Observed_Vel NumIterations_Observed_Vel ResVec_Observed_Vel M_Observed_Vel LSQ_Fail_Observed_Vel limit_Observed_Vel] = GARD_LSQVel(GPSPosLSQ,VelocityPrevious,N,PRRateMeasured_SimulatedLSQ,SatPosLSQ,SatVelLSQ);

end

