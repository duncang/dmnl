% GARDSim_WriteMOTFile

%DataPath = 'data/rnav_approach/';
%DataPath = 'data/DebugAeromodelpqrchanged11.1.08/';
%DataPath = 'data/DebugAeroModel7.1.08/';
%DataPath = 'data/Feb0108FlightNoWind/';
DataPath = 'data/rnav_approach_long/';


load(strcat(DataPath,'pos_truth_ecef.mat'));
load(strcat(DataPath,'vel_truth.mat'));
load(strcat(DataPath,'att_truth.mat'));

load(strcat(DataPath,'sensors_clean.mat'));
sensors = sensors_clean;
clear 'sensors_clean';

OutputFileName = strcat(DataPath,'simremote.csv');

GARD_WriteSimGenMotionFile(OutputFileName,pos_truth_ecef,vel_truth,att_truth,sensors);

clear pos_truth_ecef vel_truth att_truth sensors;
