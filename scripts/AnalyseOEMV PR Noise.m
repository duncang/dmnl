% 
% 
% 
% 
% %Read in pseudoranges
% 
% 
% 
% Filename = 'GPSl0770.txt'
% 
% 
% 
% [GPSTime_Week, GPSTime_Sec,NumberRinexObsTypes,ApproxPos, DATA_STRUCT] = ReadRinexGRS(Filename)
% 
% 
% 
% 
% %calculate satellite position
% 
% 
% 
% 
% 
% %True user position
% 
% 
% 
% 
% %least squares best fit
% 
% 
% Signal = DATA_STRUCT.C1(8,:)
%  [CurveFit, BestFit] = LeastSquaresBestFit(Signal, 2, 1);
%  
%  
 SP3FILE = 'RAPID';
 %read SP3 file
 
 
 for i = 1:100
 
  [NumberSVs, VehicleIDs, EpochNumber, GPSTime, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data] = readsp3(SP3FILE);
  
  
 end