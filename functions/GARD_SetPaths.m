
% Setup the paths for gardsim.   This script assumes that you are in the
% gard directory, and that the entire ARCAA repository is checked out

% Please add other subdirectories as required.
GardsimPath = pwd;
Gardsim_UsefulFunctions_Path = strcat(pwd, '/UsefulFunctions');
LogFileReaderPath = strcat(pwd, '\..\..\..\..\datacollection\LogFileReader');
MatlabNMEAPath = strcat(pwd, '/../MatlabNMEA');

% Set the paths
path(path, GardsimPath);
path(path, Gardsim_UsefulFunctions_Path);
path(path, LogFileReaderPath);
path(path, MatlabNMEAPath);

