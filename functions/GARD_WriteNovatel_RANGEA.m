function GARD_WriteNovatel_RANGEA(rangedata,outputfile)



%% This is how the data is read
%   gps.RangeData(RangeRecord).rtTimeStamp = linedata(2);
%   gps.RangeData(RangeRecord).GPSWeek = linedata(3);
%   gps.RangeData(RangeRecord).GPSSec = linedata(4);
%   gps.RangeData(RangeRecord).lNumberObservations = linedata(5);
% 
%   for i = 1:gps.RangeData(RangeRecord).lNumberObservations
%      gps.RangeData(RangeRecord).Obs(i).usPRN = linedata(5+((i-1)*10)+1);
%      gps.RangeData(RangeRecord).Obs(i).usGlonassFrequency = linedata(5+((i-1)*10)+2);
%      gps.RangeData(RangeRecord).Obs(i).dPseudorange = linedata(5+((i-1)*10)+3);
%      gps.RangeData(RangeRecord).Obs(i).fPseudorangeSigma = linedata(5+((i-1)*10)+4);
%      gps.RangeData(RangeRecord).Obs(i).dCarrierPhase = linedata(5+((i-1)*10)+5);
%      gps.RangeData(RangeRecord).Obs(i).fCarrierPhaseSigma = linedata(5+((i-1)*10)+6);
%      gps.RangeData(RangeRecord).Obs(i).fDoppler = linedata(5+((i-1)*10)+7);
%      gps.RangeData(RangeRecord).Obs(i).fCNo = linedata(5+((i-1)*10)+8);
%      gps.RangeData(RangeRecord).Obs(i).fLockTime = linedata(5+((i-1)*10)+9);
%      gps.RangeData(RangeRecord).Obs(i).ulTrackingStatus = linedata(5+((i-1)*10)+10);
%   end


%% for each record, write line of ascii.  output should look like:

% HEADER:           #RANGEA,FILE,0,58.0,FINESTEERING,1541,269619.000,00400000,5103,4158;
% Number Meas.      28,
% measurements:     29,0,20726731.627,0.064,-108974053.852023,0.006,-1703.984,50.9,34.928,08109c04,
%                   29,0,20726728.457,0.121,-84914985.463813,0.008,-1327.785,47.7,20.460,01309c0b,
%                   22,0,22364074.985,0.073,-117578350.433353,0.007,2504.965,49.0,40.920,08109c24,
%                   22,0,22364069.148,0.211,-91619624.137904,0.008,1951.922,41.6,26.780,01309c2b,
%                   26,0,22447130.079,0.074,-118014811.285210,0.007,2579.305,48.9,40.640,08109c44,
%                   26,0,22447127.668,0.205,-91959731.211816,0.008,2009.848,41.9,26.880,01309c4b,
%                   21,0,20809439.836,0.059,-109408687.687882,0.006,-1.430,50.8,40.970,08109c64,
%                   21,0,20809434.869,0.137,-85253655.732017,0.007,-1.113,45.4,26.780,01309c6b,
%                   30,0,24647173.309,1.131,-129521766.701673,0.035,-2917.320,41.8,2.766,08009c84,
%                   18,0,20282554.984,0.060,-106639892.376733,0.006,763.918,50.6,41.490,08109ca4,
%                   18,0,20282551.075,0.131,-83096152.103046,0.007,595.262,45.8,26.780,01309cab,
%                   16,0,23354078.818,0.099,-122780851.913145,0.010,418.469,46.4,40.640,08109ce4,
%                   16,0,23354074.317,0.265,-95673526.537322,0.011,326.078,39.7,26.780,01309ceb,
%                   24,0,21765467.328,0.070,-114432649.581148,0.008,-1211.422,49.6,39.180,08109d24,
%                   24,0,21765465.878,0.170,-89168442.247743,0.009,-943.969,43.5,26.780,01309d2b,
%                   15,0,22452858.174,0.074,-118044907.382894,0.008,-1031.930,49.0,39.780,08109d84,
%                   15,0,22452855.344,0.148,-91983183.831472,0.009,-804.102,44.7,26.780,01309d8b,
%                   40,12,24025382.734,0.869,-128609748.886036,0.016,3556.289,42.1,5.890,08119604,
%                   40,12,24025386.634,1.902,-100029846.883690,0.010,2765.996,35.4,5.896,00b19e0b,
%                   45,13,23823453.521,1.374,-127573502.115455,0.024,-3863.336,38.5,5.900,08019644,
%                   39,8,20073531.943,0.332,-107304529.226592,0.007,1442.777,50.5,5.900,08119704,
%                   39,8,20073535.748,0.542,-83459094.218913,0.003,1122.160,46.3,5.900,10b19f0b,
%                   61,9,19384458.452,0.321,-103657409.226413,0.006,1296.902,50.6,5.900,08119724,
%                   61,9,19384460.231,0.766,-80622437.309534,0.002,1008.699,43.2,5.896,00b19f2b,
%                   60,10,20805184.789,0.446,-111293710.621868,0.009,-2608.469,47.8,5.900,08119744,
%                   60,10,20805187.932,0.774,-86561784.269979,0.003,-2028.812,43.1,5.896,10b19f4b,
%                   54,11,22267953.505,0.639,-119160310.312140,0.011,3757.102,44.6,5.900,08119764,
%                   54,11,22267953.639,0.899,-92680247.809782,0.004,2922.188,41.8,5.896,00b19f6b
% Checksum:         *070f6970
%
%
% Each measurment field contains:
%
%   PRN                 29,
%   Glonass Freq        0,
%   Pseudorange         20726731.627,
%   Pseudorange-S       0.064,
%   Carrier Phase (ADR) -108974053.852023,
%   Carrier Phase-S     0.006,
%   Doppler (Hz)        -1703.984,
%   C/No                50.9,
%   Locktime            34.928,
%   Channel Tracking Status 08109c04,
%

of = fopen(outputfile,'w');

number_records = length(rangedata);

for i = 1:number_records
    
    % write header
    % write sync char
    fprintf(of,'#'); 
    outstring = sprintf('RANGEA,COM1,0,0,FINESTEERING,%d,%10.3f,00000000,5103,4158;%d', ...
                rangedata(i).GPSWeek,   ...
                rangedata(i).GPSSec,    ...
                rangedata(i).lNumberObservations);

                
    % write data
    for j = 1:rangedata(i).lNumberObservations
        outstring = sprintf('%s,%d,%d,%.3f,%.3f,%.6f,%.3f,%.3f,%.1f,%.3f,%08x', outstring ,      ...
                rangedata(i).Obs(j).usPRN,        ...
                rangedata(i).Obs(j).usGlonassFrequency,        ...
                rangedata(i).Obs(j).dPseudorange,        ...
                rangedata(i).Obs(j).fPseudorangeSigma,        ...
                rangedata(i).Obs(j).dCarrierPhase,        ...
                rangedata(i).Obs(j).fCarrierPhaseSigma,        ...
                rangedata(i).Obs(j).fDoppler,        ...
                rangedata(i).Obs(j).fCNo,        ...
                rangedata(i).Obs(j).fLockTime,        ...
                rangedata(i).Obs(j).ulTrackingStatus);
    end
    
    % write checksum
    crc = novatelCRC(outstring);
    outstring = sprintf('%s*%08x\r\n',outstring,crc);
    
    
    fprintf(of,'%s',outstring);
    
    
end % for i=1:number_records

fclose(of);


