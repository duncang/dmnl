% GARDSIM_FlightModel
%
% Written by Duncan Greer (c) 2005 GARD Project
% Created 27 September 2005
%
% this funciton/script? generates a set of pseudorange measurements based on the
% flight position data provided by teh Aerosim model. The simulated PRs are
% then fed into a position solution function (GARD_LSQ) and a GPS position
% solution is calculated.
%
% $Id: GARDSIm_FlightModel.m 1877 2008-07-15 05:00:01Z n2523710 $
%
GPSConstants; %All generic constants put in this script, as global variables;
MaskAngle = 7.5; %degrees
% load data to the workspace

DataPath = 'data/rnav_approach/';
%DataPath = 'data/DebugAeromodelpqrchanged11.1.08/';
%DataPath = 'data/DebugAeroModel7.1.08/';
%DataPath = 'data/Feb0108FlightNoWind/';
%DataPath = 'data/rnav_approach_long/';


FlightDataFile = strcat(DataPath,'pos_truth_ecef.mat');
load(FlightDataFile);

VelocityDataFile = strcat(DataPath,'vel_truth.mat');
load(VelocityDataFile);

DataRate = 100; %hz
% TODO:  receiver clock model to generate a time offset


% load the SP3 precise data for this day
FilenameSP3 = 'data/sp3/NG05JA05.SP3';
[NumberSVs_SP3, VehicleIDs_SP3, EpochNumber_SP3, GPSTime_SP3, SV_X_Data_SP3, SV_Y_Data_SP3, SV_Z_Data_SP3, SV_T_Data_SP3] = readsp3(FilenameSP3);

% interpolate sp3 orbit to 1 second epochs for 1 hours
for (SV=1:length(VehicleIDs_SP3))
    if(VehicleIDs_SP3(SV) == 1)
        % get the orbit for teh satellite if the data is available
        SV_X_Data_True(:,SV) = spline([1:15*60:12*3600],SV_X_Data_SP3(1:4*12,SV),[1:1:12*3600]);
        SV_Y_Data_True(:,SV) = spline([1:15*60:12*3600],SV_Y_Data_SP3(1:4*12,SV),[1:1:12*3600]);
        SV_Z_Data_True(:,SV) = spline([1:15*60:12*3600],SV_Z_Data_SP3(1:4*12,SV),[1:1:12*3600]);
        SV_T_Data_True(:,SV) = spline([1:15*60:12*3600],SV_T_Data_SP3(1:4*12,SV),[1:1:12*3600]);
    end
end

% load a nav file and find visible satellites for the desired epochs

FilenameNav = 'data/rinex/EESE0051_mod.05N';

NavigationData = freadnav(FilenameNav);
% iono model parameters for the above Nav file
%ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA           
%BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA  

% find the nav data time toe and generate the desired start and stop times
% that we can calculate the prs for.

GPS_Week_Start = NavigationData(1,32);
GPS_Second_Start = round(NavigationData(1,2));

SimulationLength = pos_truth_ecef(1,length(pos_truth_ecef)); % the last time stamp in seconds 

GPS_Second_End = GPS_Second_Start + SimulationLength;

% check for week rollover
if(GPS_Second_End > 302400)
    GPS_Second_End = GPS_Second_End - 604800;
    GPS_Week_End = GPS_Week_Start + 1;
else
    GPS_Week_End = GPS_Week_Start;
end

% Initial User Clock Bias
InitialUserClockBias = 0;

% for each time interval epoch of 1 second, calculate a PR to match the
% vehicle position given by the pos_truth_ecef data.

% we will need to predict all of hte sattellite positions, and apply a mask
% angle to determine which satellites would be visible to the flight.  this
% will be time varying.. eep!!]
% TODO: preallocate matrices for SPEED

NumberEpochs = floor(SimulationLength);
NumberSVs = 32;  % total number of SVs (by max PRN) in teh constellation

PR_Sim = zeros([NumberEpochs,NumberSVs]);
PRR_Sim = zeros([NumberEpochs,NumberSVs]);
CP_Sim = zeros([NumberEpochs,NumberSVs]);
NAmb_Sim = zeros([NumberEpochs,NumberSVs]);
GeometricRange_Sim = zeros([NumberEpochs,NumberSVs]);

% read in the ionosphere data

% get the iono tec map
IonexFilename = 'data/ionex/igsg0050_mod.05i';
if ~exist('IonoTECMap')
    disp(sprintf('[GARDSim_FlightModel] Reading Ionosphere data from %s',IonexFilename));
    [IonoTECMap, GPSTime_IonoObs, IonoLongitude, IonoLatitude] = ReadIonex(IonexFilename);
end
% use the 6th iono map for this simulation
IonoTECMap_Sim = IonoTECMap(:,:,6);
IonoLatitude_Sim = IonoLatitude(6,:);
IonoLongitude_Sim = IonoLongitude(6,:);

time_step = 1;

disp(sprintf('[GARDSim_FlightModel] Generating GPS measurements for %d Epochs',SimulationLength));
for Epoch = 1:SimulationLength
    SVTime = GPS_Second_Start + Epoch - 1;  % assuming 1 second epochs
    SVWeek = GPS_Week_Start;
    if (SVTime > 302400)
        SVTime = SVTime - 604800;
        SVWeek = SVWeek + 1;
    end
     
    GPSTime_Week(Epoch) = SVWeek;
    GPSTime_Sec(Epoch) = SVTime;
    
    % This uses the orbit propogator to predict SV locations.
    % TODO: should actually use SP3 precise ephemeris for the day instead. 
    for SV=1:32
        [SV_X_Data(Epoch,SV) SV_Y_Data(Epoch,SV) SV_Z_Data(Epoch,SV) SV_T_Data(Epoch,SV) ValidData(Epoch,SV)] = GPSOrbitPropagator(SVWeek, SVTime, SV, NavigationData);
        [SV_Xvel_Data(Epoch,SV), SV_Yvel_Data(Epoch,SV), SV_Zvel_Data(Epoch,SV), SV_Tvel_Data(Epoch,SV), SV_Xacc_Data(Epoch,SV), SV_Yacc_Data(Epoch,SV), SV_Zacc_Data(Epoch,SV), SV_Tacc_Data(Epoch,SV), ValidVelData(Epoch,SV)] = GPSOrbitPropagatorVelocities(SVWeek,SVTime, SV, NavigationData);
    end
    
    % find the user position,  and LOS vector to each satellite -data rate
    % is 100 Hz
    UserPos_X(Epoch) = pos_truth_ecef(2,(Epoch * DataRate) + 1);
    UserPos_Y(Epoch) = pos_truth_ecef(3,(Epoch * DataRate) + 1);
    UserPos_Z(Epoch) = pos_truth_ecef(4,(Epoch * DataRate) + 1);

    
    % find user lat, long, height
    [UserLatitude(Epoch), UserLongitude(Epoch),UserHeight(Epoch)] = ECEF2LLH([UserPos_X(Epoch), UserPos_Y(Epoch), UserPos_Z(Epoch)]);
    
    % add clock model
    if Epoch == 1
        [UserPos_T(Epoch),UserClockDrift(Epoch)] = GARD_ClockModel(InitialUserClockBias);
    else
        [UserPos_T(Epoch),UserClockDrift(Epoch)] = GARD_ClockModel(UserPos_T(Epoch-1)/c);
    end
    
   UserPos_T(Epoch) = UserPos_T(Epoch) * c;
    %UserClockDrift(Epoch) = UserClockDrift(Epoch) * c;
    
    
    % calculate velocity
    Tecef2ned = T_ECEF2NED(UserLatitude(Epoch),UserLongitude(Epoch));
    Tned2ecef = Tecef2ned';
    UserVel = Tned2ecef * vel_truth(2:4,(Epoch*DataRate)+1);
    UserVel_X(Epoch) = UserVel(1);
    UserVel_Y(Epoch) = UserVel(2);
    UserVel_Z(Epoch) = UserVel(3);
    UserVel_T(Epoch) = 0.0;
    
    for SV=1:32
        [SV_Azimuth(Epoch,SV), SV_Elevation(Epoch,SV)] = AzEl([UserPos_X(Epoch) UserPos_Y(Epoch) UserPos_Z(Epoch)], [SV_X_Data(Epoch,SV) SV_Y_Data(Epoch,SV) SV_Z_Data(Epoch,SV)]);

            
        % if the sv is above teh mask angle, calculate the predicted
        % pseudorange

        if(SV_Elevation(Epoch,SV) > deg2rad(MaskAngle))
            
            % calculate the ionospheric delay
            IonoDelay(Epoch,SV) = GARD_IonoDelay(SV_Elevation(Epoch,SV), SV_Azimuth(Epoch,SV), UserLatitude(Epoch), UserLongitude(Epoch), IonoTECMap_Sim, IonoLatitude_Sim, IonoLongitude_Sim);
            
            
            if Epoch == 1
                PR_previous = 0;
%                CP_previous = 0;
%                ICP_previous = 0;
            else
                PR_previous = PR_Sim(Epoch-1,SV);
%                CP_previous = CP_Sim(Epoch-1,SV);
%                ICP_previous = ICP_Sim(Epoch-1,SV);
            end
           % [PRMeasPred, PRRMeasPredict, ICPMeasPredict,
           % CPMeasPredict,NAmb, GeometricRange] = 
           %GARDSim_PRPredict(SatPos,UserPos, IonoDelay, time_step, PRMeas_previous, ICPMeas_previous,CPMeas_previous,GeometricRange_Previous,EstimatedClockBias);
           
           if(ValidData(Epoch,SV) == 1)
               
%                [PR_Sim(Epoch,SV), PRR_Sim(Epoch,SV), ICP_Sim(Epoch,SV), CP_Sim(Epoch,SV), NAmb_Sim(Epoch,SV), GeometricRange_Sim(Epoch,SV)] =  ...
%                           GARDSim_PRPredict([SV_X_Data_True(Epoch,SV) SV_Y_Data_True(Epoch,SV) SV_Z_Data_True(Epoch,SV) SV_T_Data_True(Epoch,SV)], ...
%                                             [UserPos_X(Epoch) UserPos_Y(Epoch) UserPos_Z(Epoch) UserPos_T(Epoch)], ...
%                                             IonoDelay(Epoch,SV), time_step, PR_previous, ICP_previous, CP_previous);
%                 [PR_Sim(Epoch,SV), PRR_Sim(Epoch,SV), ICP_Sim(Epoch,SV), CP_Sim(Epoch,SV), NAmb_Sim(Epoch,SV), GeometricRange_Sim(Epoch,SV)] =  ...
%                           GARDSim_PRPredict([SV_X_Data(Epoch,SV) SV_Y_Data(Epoch,SV) SV_Z_Data(Epoch,SV) SV_T_Data(Epoch,SV)], ...
%                                             [UserPos_X(Epoch) UserPos_Y(Epoch) UserPos_Z(Epoch) UserPos_T(Epoch)], ...
%                                             IonoDelay(Epoch,SV), time_step, PR_previous, ICP_previous, CP_previous);
                [PR_Sim(Epoch,SV), PRR_Sim(Epoch,SV), GeometricRange_Sim(Epoch,SV)] = ...
                        GARDSim_PRPredictSimulated([SV_X_Data(Epoch,SV) SV_Y_Data(Epoch,SV) SV_Z_Data(Epoch,SV) SV_T_Data(Epoch,SV)], ...
                                                   [SV_Xvel_Data(Epoch,SV) SV_Yvel_Data(Epoch,SV) SV_Zvel_Data(Epoch,SV) SV_Tvel_Data(Epoch,SV)], ...
                                                   [UserPos_X(Epoch) UserPos_Y(Epoch) UserPos_Z(Epoch) UserPos_T(Epoch)], ...
                                                   [UserVel_X(Epoch) UserVel_Y(Epoch) UserVel_Z(Epoch) UserVel_T(Epoch)], ...
                                                   0,0, ...
                                                   PR_previous);
                                        
                % check if we have new SV measurements - if so, the PRR_Sim
                % will be the size of the PR since the previous PR was zero
                if(PR_previous == 0)
                    PRR_Sim(Epoch,SV) = 0;
                end
           end
        end
    end
    % add on varous errors to the PR prediction such as iono/tropo errors,
    % receiver noise, multipath etc

    if(mod(Epoch,100) == 0)
        disp(sprintf('[GARDSim_FlightModel] Generated GPS Measurements for %d Epochs',Epoch));
    end

end


% TODO: output the pr measurements to a measurements vector for furhter
% processing

SV_Ephemeris = NavigationData;

outputfile = strcat(DataPath,'PR_Simulation.mat');
save(outputfile, ...
        'UserPos_X', 'UserPos_Y', 'UserPos_Z', 'UserPos_T', ...
        'UserVel_X', 'UserVel_Y', 'UserVel_Z', 'UserVel_T', ...
        'PR_Sim', 'PRR_Sim', ...
        'SV_X_Data', 'SV_Y_Data', 'SV_Z_Data', 'SV_T_Data', ...
        'SV_Xvel_Data', 'SV_Yvel_Data', 'SV_Zvel_Data', 'SV_Tvel_Data', ...
        'SV_Xacc_Data', 'SV_Yacc_Data', 'SV_Zacc_Data', 'SV_Tacc_Data', ...
        'ValidData', 'ValidVelData', 'GPSTime_Week','GPSTime_Sec','SV_Ephemeris');

% perform position solution
disp(sprintf('[GARDSim_FlightModel] Evaluating flight for %d Epochs',NumberEpochs));
for Epoch=1:NumberEpochs
    
    % use the positoin solution from the previous epoch as the predicted
    % position for the next solution
    
    % format the input vectors
    SVIndex = 0;
    for SV=1:NumberSVs
        if(PR_Sim(Epoch,SV) ~= 0)
            % add to PR vector
            SVIndex = SVIndex + 1;
            PR_Vec(SVIndex) = PR_Sim(Epoch,SV);
            SV_Vec(SVIndex,:) = [SV_X_Data(Epoch,SV) SV_Y_Data(Epoch,SV) SV_Z_Data(Epoch,SV) SV_T_Data(Epoch,SV)];
        end
    end
    
    
    
    if(Epoch == 1)
          UserPos = [UserPos_X(1) UserPos_Y(1) UserPos_Z(1) UserPos_T(1)];
    else
        UserPos = SolutionVec(Epoch-1,:);
    end
    
    % Perform Least-Squares position solution
    [SolutionVec(Epoch,:) VarSolutionVec(Epoch,:) NumIterations ResidualVector(Epoch,1:SVIndex) M LSQ_Fail(Epoch,:) limit DOP(Epoch,:)] = GARD_LSQ(UserPos,SVIndex,PR_Vec,SV_Vec);


    % perform RAIM Check
    [BadGeometry(Epoch), RAIM_ALARM(Epoch), SLOPE_Max(Epoch), r, rd, ARP(Epoch)] = GARD_RAIM(SVIndex,1/15000,7.5,300,ResidualVector(Epoch,:),M);
    
    % formulate solution error function
    SolutionErr_X(Epoch) = pos_truth_ecef(2,Epoch*DataRate+1) - SolutionVec(Epoch,1);
    SolutionErr_Y(Epoch) = pos_truth_ecef(3,Epoch*DataRate+1) - SolutionVec(Epoch,2);
    SolutionErr_Z(Epoch) = pos_truth_ecef(4,Epoch*DataRate+1) - SolutionVec(Epoch,3);
    
    % reset LSQ variables
    clear PR_Vec SV_Vec
    
    % 
    if(mod(Epoch,100) == 0)
        disp(sprintf('[GARDSim_FlightModel] Completed evaluation for %d Epochs',Epoch));
    end
end

%plot3(SolutionVec(:,1),SolutionVec(:,2),SolutionVec(:,3),'r*');
%grid on;

% convert to LLH
disp('[GARDSim_FlightModel] Converting results to LLH');
for Epoch=1:NumberEpochs
    [Lat(Epoch),Long(Epoch),Height(Epoch)] = ECEF2LLH(SolutionVec(Epoch,:));
    Lat(Epoch) = rad2deg(Lat(Epoch));
    Long(Epoch) = rad2deg(Long(Epoch));
    
    [Lat_truth(Epoch),Long_truth(Epoch),Height_truth(Epoch)] = ECEF2LLH(pos_truth_ecef(2:4,Epoch*DataRate+1));
    Lat_truth(Epoch) = rad2deg(Lat_truth(Epoch));
    Long_truth(Epoch) = rad2deg(Long_truth(Epoch));
    
    % rotate the x,y, z errors to ENU coords
    TMatrix_ECEF2ENU = T_ECEF2ENU(Long(Epoch),Lat(Epoch));
    
    Err_XYZ = [SolutionErr_X(Epoch);SolutionErr_Y(Epoch);SolutionErr_Z(Epoch)];
    Err_ENU = TMatrix_ECEF2ENU * Err_XYZ;
    SolutionErr_E(Epoch) = Err_ENU(1);
    SolutionErr_N(Epoch) = Err_ENU(2);
    SolutionErr_U(Epoch) = Err_ENU(3);
    
end

disp('[GARDSim_FlightModel] Done!');

%figure();


% plot(SolutionVec(:,1),'r');
% hold on;
% plot(pos_truth_ecef(1,:),pos_truth_ecef(2,:),'b');
% xlabel('X Position (m)');
% ylabel('Y Position (m)');
% 
% 

% load the australian coastline and airports
if(~exist('AustAirports'))
    load data/AustAirports.mat;
end

if ~exist('austcoast')
    load 'data/austcoast.dat';
end

plot(Long,Lat,'r*');
hold on;
plot(austcoast(:,1),austcoast(:,2),'k');
plot(AustAirports(:,6)*180/pi,AustAirports(:,5)*180/pi,'k.');
plot(Long_truth,Lat_truth,'b');
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
grid on;
hold off;
% 
% % error plot
figure();
title('Simulated Pseudorange Based Position Solution');
subplot(3,1,1), plot(SolutionErr_X); grid on; ylabel('X Residual Error (m)');
subplot(3,1,2), plot(SolutionErr_Y); grid on; ylabel('Y Residual Error (m)');
subplot(3,1,3), plot(SolutionErr_Z); grid on; ylabel('Z Residual Error (m)');
xlabel('Time (s)');

figure();
plot(UserPos_T,'b')
hold on
plot(SolutionVec(:,4),'r')
grid on
xlabel('Simulation Time (s)');
ylabel('Receiver Clock Bias (m)');
legend('Simulated Clock Bias','LSQ Estimated Clock Bias');


% plot the residual between the clock bias estimate adn the real clock bias
figure;
residual = SolutionVec(:,4) - UserPos_T';
plot(residual);
grid on
xlabel('Simulation Time (s)');
ylabel('Clock Bias Estimate Residual (m)');
% 
figure();
plot(ARP)
grid on;
xlabel('Simulation Time (s)')
ylabel('Approximate Radial Error Protected (metres)');
ylabel('ARP (metres)');
title('Approximate Radial Error Protected for Simulated Flight');

