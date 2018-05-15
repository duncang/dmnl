
ReadSimTruth = 0;
ReadBestXYZ = 0;
ReadObs = 1;
ReadNav = 1;

%DataPath = 'data/Simulator_Data/Apr_08/data/sg_rnav_approach/';
%FileNamePrefix = 'novatel_080410022139';

%DataPath = 'data/Simulator_Data/Apr_08/data/rnav_approach_long/';
%FileNamePrefix = 'novatel_080410234528';

DataPath = 'data/Simulator_Data/Apr_08/data/rnav_approach_long/';
FileNamePrefix = 'novatel_080411050311';



Year = '05';

NavFile = strcat(DataPath,FileNamePrefix,'.',Year,'N');
ObsFile = strcat(DataPath,FileNamePrefix,'.',Year,'O');





if ReadBestXYZ == 1
    % read bestxyz data from ASC file
end

if ReadObs == 1
    % load observation and nav data
    [GPSTime_Week, GPSTime_Sec,NumberRinexObsTypes,ValidDataRinexObs,ApproxPos,...
        Novatel_C1, Novatel_L1, Novatel_D1, Novatel_S1, Novatel_P2, Novatel_L2, ...
        Novatel_D2, Novatel_S2] = ReadRinexNovatel(ObsFile);


    NumberEpochsGPS = size(Novatel_C1,2);
    gps_dt = 1;

    NumberSVs = size(Novatel_C1,1);
    SVDontUse = zeros(1,NumberSVs);


end

if ReadNav == 1
    [SV_Ephemeris, ALPHA, BETA] = freadnav(NavFile);
end



% propogate satellite orbits and velocities
% disp('Calculating Satellite Orbits...');
% 
% % preallocate storage vectors for speed
% SV_X_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Y_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Z_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_T_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Xvel_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Yvel_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Zvel_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Tvel_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Xacc_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Yacc_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Zacc_Data = zeros(NumberEpochsGPS,NumberSVs);
% SV_Tacc_Data = zeros(NumberEpochsGPS,NumberSVs);
% 
% for i=1:NumberEpochsGPS
%     for j=1:NumberSVs
%         [SV_X_Data(i,j), SV_Y_Data(i,j), SV_Z_Data(i,j), SV_T_Data(i,j), ValidPosData(i,j)] = ...
%                GPSOrbitPropagator(GPSTime_Week(i), GPSTime_Sec(i), j, SV_Ephemeris, 7500);
%         [SV_Xvel_Data(i,j), SV_Yvel_Data(i,j), SV_Zvel_Data(i,j), SV_Tvel_Data(i,j), ...
%          SV_Xacc_Data(i,j), SV_Yacc_Data(i,j), SV_Zacc_Data(i,j), SV_Tacc_Data(i,j), ValidVelData(i,j)] = ...
%          GPSOrbitPropagatorVelocities(GPSTime_Week(i),GPSTime_Sec(i), j, SV_Ephemeris);
%     end
%     
%     if(mod(i,1000) == 0)
%         disp(sprintf('Completed Epoch %d',i));
%     end
% end


if ReadSimTruth == 1
    % read sim truth data

    SimTruthFile = 'motion_V1-mod.csv';
    sim_truth = dlmread(strcat(DataPath,SimTruthFile),',');
    sim_time = sim_truth(:,1)/1000.0;

    pos_truth_llh(1,:) = sim_time';
    pos_truth_llh(2,:) = sim_truth(:,14)';
    pos_truth_llh(3,:) = sim_truth(:,15)';
    pos_truth_llh(4,:) = sim_truth(:,16)';


    % convert simtruth time (col 1) to gps seconds
    for i=1:length(sim_truth)
        hrs_in = floor(sim_truth(i,1)/3600/1000);
        min_in = floor(sim_truth(i,1)/60/1000);
        sec_in = sim_truth(i,1)/1000 - hrs_in*3600 - min_in * 60;
        [gpsweek_out(i), gpssec_out(i)] = ftime([07 03 01 hrs_in min_in sec_in]);
        sim_truth(i,1) = gpssec_out(i);
    end

end
