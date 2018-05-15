function Data = GARD_LoadRinexDirectory(Directory)
% function Data = GARD_LoadRinexDirectory(Directory)
% This script loads all of the rinex data contained in the directory
% specified by 'Directory' and returns the data structure 'Data'.
% Written by Duncan Greer 12 April 2007
% $Id: GARD_LoadRinexDirectory.m 1850 2008-07-14 04:52:47Z greerd $
%



%% first get a listing of all observation and nav files in the directory
nav_files = dir(fullfile(Directory,'*.*N'));
number_nav_files = size(nav_files,1);
obs_files = dir(fullfile(Directory,'*.*O'));
number_obs_files = size(obs_files,1);


% provide a warning if hte number of obs files and nav files does not match
% up
if number_nav_files ~= number_obs_files
    disp(sprintf('Warning: %d Nav files and %d Obs files',number_nav_files,number_obs_files));
end


%% loop through nav files
for i = 1:number_nav_files
    Filename = fullfile(Directory,nav_files(i).name);
    disp(sprintf('Processing %s',Filename));
    [nav_data, ALPHA, BETA] = freadnav(Filename);
    Data.NavigationData(i).Ephemeris = nav_data;
    Data.NavigationData(i).Iono_ALPHA = ALPHA;
    Data.NavigationData(i).Iono_BETA = BETA;
    
    clear nav_data ALPHA BETA;
end

for i = 1:number_obs_files
    Filename = fullfile(Directory,obs_files(i).name);
    disp(sprintf('Processing %s',Filename));
    [GPSTime_Week, GPSTime_Sec,NumberRinexObsTypes,ApproxPos, DATA_STRUCT] = ReadRinexGRS(Filename);
    Data.ObsData(i).GPSTime_Week = GPSTime_Week;
    Data.ObsData(i).GPSTime_Sec = GPSTime_Sec;
    Data.ObsData(i).NumberRinexObsTypes = NumberRinexObsTypes;
    Data.ObsData(i).ApproxPos = ApproxPos;
    Data.ObsData(i).Observations = DATA_STRUCT;
    
    clear GPSTime_Week GPSTime_Sec NumberRinexObsTypes ApproxPos DATA_STRUCT;
end