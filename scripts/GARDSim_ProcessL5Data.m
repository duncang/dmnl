
% rangedata = GARD_ReadRANGE('data/novatel_9july10/friday 9 july.ASC',2)
% [navdata, ION_ALPHA, ION_BETA] = freadnav('data/novatel_9july10/friday 9 july.10N')

%load data/novatel_9july10/navdata.mat;
%load data/novatel_9july10/rangedata.mat;

GPSConstants;

% look through data and find L5 records - we only have 1 prn to worry about, so no need for a 3d array yet

for i=1:size(rangedata,2)
    
    for j = 1:rangedata(i).NumberMeasurements 
        if rangedata(i).Measurements(j).signal_type == 14 && rangedata(i).Measurements(j).prn == 25
            L5Data(i,1) = rangedata(i).GPSSec;
            L5Data(i,2) = rangedata(i).Measurements(j).prn;
            L5Data(i,3) = rangedata(i).Measurements(j).pseudorange;
            L5Data(i,4) = rangedata(i).Measurements(j).carrier_phase;
            L5Data(i,5) = rangedata(i).Measurements(j).doppler;
            L5Data(i,6) = rangedata(i).Measurements(j).CNo;
            L5Data(i,7) = rangedata(i).Measurements(j).locktime;

        end

        if rangedata(i).Measurements(j).signal_type == 0 && rangedata(i).Measurements(j).prn == 25
            L1Data(i,1) = rangedata(i).GPSSec;
            L1Data(i,2) = rangedata(i).Measurements(j).prn;
            L1Data(i,3) = rangedata(i).Measurements(j).pseudorange;
            L1Data(i,4) = rangedata(i).Measurements(j).carrier_phase;
            L1Data(i,5) = rangedata(i).Measurements(j).doppler;
            L1Data(i,6) = rangedata(i).Measurements(j).CNo;
            L1Data(i,7) = rangedata(i).Measurements(j).locktime;

        end
        
        if rangedata(i).Measurements(j).signal_type == 9 && rangedata(i).Measurements(j).prn == 25
            L2Data(i,1) = rangedata(i).GPSSec;
            L2Data(i,2) = rangedata(i).Measurements(j).prn;
            L2Data(i,3) = rangedata(i).Measurements(j).pseudorange;
            L2Data(i,4) = rangedata(i).Measurements(j).carrier_phase;
            L2Data(i,5) = rangedata(i).Measurements(j).doppler;
            L2Data(i,6) = rangedata(i).Measurements(j).CNo;
            L2Data(i,7) = rangedata(i).Measurements(j).locktime;

        end
        
    end
    
end


% figure(); grid on; hold on;
% plot(L5Data(:,6))
% xlabel('Track Time (sec)');
% ylabel('CNo (dB-Hz)')



% lets look at the dual frquency correction
f1 = 1575.42e6;
f2 = 1227.60e6;
f5 = 1176.45e6;
alpha_L1L5 = 1 - f1^2 / f5^2;
alpha_L1L2 = 1 - f1^2 / f2^2;


start = 1;
stop = 10000;
SV = 25;
UserPos = [-5053036.7962,2562955.2710,-2919254.1297,0];
[UserPos_LLH(1),UserPos_LLH(2),UserPos_LLH(3)] = ECEF2LLH(UserPos(1:3));

for i=start:stop
    % check times match
    if L1Data(i,1) ~= L5Data(i,1)
        disp('L1 Time Mismatch');
    end
    if L1Data(i,1) ~= L2Data(i,1)
        disp('L2 Time Mismatch');
    end
    
    
    % calculate orbits
    [SVPos(i,1), SVPos(i,2), SVPos(i,3), SVPos(i,4), ValidPosData(SV), URA(i,SV)] = ...
               GPSOrbitPropagator(rangedata(i).GPSWeek, L1Data(i,1) - L1Data(i,3)/Speedoflight, SV, navdata, 7500);
    SVPos(i,4) = SVPos(i,4) * Speedoflight;

    [SV_Azimuth(i,SV), SV_Elevation(i,SV)] = AzEl(UserPos(1:3), SVPos(i,1:3));
    
    
    % calculate the iono delay correction - single frequency user
    % model from ICD 200

    ionodelay(i,SV) = ionomodel(L1Data(i,1), UserPos(1:3), SVPos(i,1:3), ION_ALPHA, ION_BETA);
    TropoDelay(i,SV) =  GARD_TropoDelay(SV_Elevation(i,SV),UserPos_LLH(3));
    EarthRotation(i,SV) = -(OMEGAedot / Speedoflight) * (SVPos(i,1) * UserPos(2) - SVPos(i,2) * UserPos(1));
    
    pr_diff_L1L5(i) = (L1Data(i,3) - L5Data(i,3));
    pr_diff_L1L2(i) = (L1Data(i,3) - L2Data(i,3));
    pr_correction_L1L5(i) = 1/alpha_L1L5 * pr_diff_L1L5(i);
    pr_correction_L1L2(i) = 1/alpha_L1L2 * pr_diff_L1L2(i);
end


