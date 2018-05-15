%  Validation of Ionospheric Model
% Written by Duncan Greer 27 July 2005


% methodology - compare predicted iono activity with estimated iono delay
% on the day

format long g;

% GPS constants
GPSConstants;

c = 2.99792458e8;  % speed of light m/s
x_true = [-5046773.309;2568446.690;-2925288.932;0];  % surveyed position of base station receiver

% get user geodetic latitude and longitude
[phi_u,lamda_u,h_u] = ECEF2LLH(x_true);


SP3DataFileName = 'data/SP3/NG05JA05.SP3';
[NumberSVs_SP3, VehicleIDs_SP3, NumberEpochs_SP3, GPSTime_SP3, SV_X_Data_SP3, SV_Y_Data_SP3, SV_Z_Data_SP3, SV_T_Data_SP3] = readsp3(SP3DataFileName);


NavDataFileName = 'data/rinex/EESE0051_mod.05N';
NavigationData = freadnav(NavDataFileName);

% iono model parameters for the above Nav file
ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA           
BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA  

ObservationFilename = 'data/rinex/EESE0051_mod.05O';
% load observation data if its not there already

if(~exist('GPStimeAMicro'))
    [GPStimeAMicro,ValidDataRinexObs, L1_PRNAMicro, L2_PRNAMicro, C1_PRNAMicro, P1_PRNAMicro, P2_PRNAMicro] = ReadRinexSblockRoof(ObservationFilename);
end

% PART A - for each satellite at each epoch, calculate the iono delay based
% on the ICD200 model

for Epoch = 1:NumberEpochs_SP3
    for SV = 1:length(VehicleIDs_SP3)
        % if SV data is available, calculate the iono
        if VehicleIDs_SP3(SV)
            IonoDelay(Epoch,SV) = ionomodel(GPSTime_SP3(Epoch), x_true(1:3)', [SV_X_Data_SP3(Epoch,SV) SV_Y_Data_SP3(Epoch,SV) SV_Z_Data_SP3(Epoch,SV)], ALPHA, BETA);
            [Azimuth(Epoch,SV) Elevation(Epoch,SV)] = AzEl( x_true(1:3)',[SV_X_Data_SP3(Epoch,SV) SV_Y_Data_SP3(Epoch,SV) SV_Z_Data_SP3(Epoch,SV)]);
        end
    end
end

% PART B - calculate the dual-frequency pseudorange and find the iono delay based on
% Measured L1 PR - Corrected PR

for Epoch = 1:4%NumberEpochs_SP3 % every 15 minutes
    for SV = 1%:NumberSVs_SP3
        % check if a measurement exists
        if ValidDataRinexObs(SV)
            PR1 = P1_PRNAMicro(SV,Epoch * 60 * 15);
            PR2 = P2_PRNAMicro(SV,Epoch * 60 * 15);
            PR_Corrected = DualFreqIono(PR1,PR2);
            
            IonoDelay_Calculated(Epoch,SV) = PR1 - PR_Corrected;
        end
    end
end



% PART C - calculate the iono delay based on IGS IONEX TEC Maps

% get the iono tec map
IonexFilename = 'data/ionex/igsg0050_mod.05i';
[IonoTECMap, GPSTime_IonoObs, IonoLongitude, IonoLatitude] = ReadIonex(IonexFilename);

% for each satellite at each epoch,...
IonoEpoch = 1;  % use the first TEC map

for Epoch = 1:4%NumberEpochs_SP3 % every 15 minutes - use the same TEC map for each obs for now - TODO: interpolate map
    for SV = 1%:NumberSVs_SP3
            % use azimuth and elevation calculated in part A
            
            % calculate sub-iono point
            [lamda_i(Epoch,SV), phi_i(Epoch,SV)] = getsubionopoint(lamda_u, phi_u, Elevation(Epoch,SV), Azimuth(Epoch,SV));
                        
            % lookup VTEC value from TEC Map
             long_index = find(IonoLongitude(IonoEpoch,:) > lamda_i(Epoch,SV)*180/pi, 1, 'first'); % return the first iono tec map bin that matches the iono longitude
             lat_index = find(IonoLatitude(IonoEpoch,:) > phi_i(Epoch,SV)*180/pi, 1, 'last'); % return the first inoo tec map bin that matches the iono latitude
             VTEC(Epoch,SV) = IonoTECMap(lat_index,long_index,IonoEpoch);
            
            % convert to Slant TEC
            z = pi/2 - Elevation(Epoch,SV);
            STEC(Epoch,SV) = vtec2stec(VTEC(Epoch,SV),z,6371.8,6371.8,450);
            
            %  calculate iono delay using 16.2cm per TECU
            
            IonoDelay_TECMap(Epoch,SV) = STEC(Epoch,SV) * 0.162;
        
    end
end

% print results
% for SV=1:length(VehicleIDs_SP3)
%    if VehicleIDs_SP3(SV)
%     figure;
%     subplot(2,1,1), plot(IonoDelay(1:NumberEpochs_SP3,SV));
%     title(sprintf('Ionospheric Delay SV: %d, 5 Jan 2005',SV));
%     xlabel('Time Epoch (15 mins)');
%     ylabel('Delay (m)')';
%     grid on;
%     axis([0 NumberEpochs_SP3 0 10]);
%     subplot(2,1,2), plot(Elevation(1:NumberEpochs_SP3,SV)* 180 / pi);
%     xlabel('Time Epoch (15 mins)');
%     ylabel('Elevation (deg)')';
%     grid on;
%     axis([0 NumberEpochs_SP3 0 90]);
%    end
% end

plot(IonoDelay(1:4,1),'b');
hold on;
plot(IonoDelay_Calculated(1:4,1),'r');
plot(IonoDelay_TECMap(1:4,1),'g');
grid on;
legend('Predicted Iono Delay (single f)','Calculated Iono Delay (dual f)','Calculated Iono Delay (TEC Map)');
xlabel('Epoch Number');
ylabel('Delay (m)');
