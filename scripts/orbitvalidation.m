%  Validation of Orbit Propogation
% Written by Duncan Greer 14 July 2004
% $Id: orbitvalidation.m 1884 2008-07-15 05:54:33Z n2523710 $


% methodology - compare propogated orbits with SP3 precise ephemeris and
% see if it is within spec.

% load GPS constants
GPSConstants;
%c = 2.99792458e8;  % speed of light m/s

SP3DataFileName = 'data/SP3/NG05JA05.SP3';
[NumberSVs_SP3, VehicleIDs_SP3, NumberEpochs_SP3, GPSTime_SP3, SV_X_Data_SP3, SV_Y_Data_SP3, SV_Z_Data_SP3, SV_T_Data_SP3] = readsp3(SP3DataFileName);


NavDataFileName = 'data/rinex/EESE0051_mod.05N';
NavigationData = freadnav(NavDataFileName);

% ObservationFilename = 'data/EESE0051_mod.05O';
% [GPStimeAMicro, L1_PRNAMicro, L2_PRNAMicro, C1_PRNAMicro, P1_PRNAMicro, P2_PRNAMicro] = ReadRinexSblockRoof(ObservationFilename);

% every 15 minutes calculate the SV orbit positions and compare with SP3

for EpochNumber = 1:NumberEpochs_SP3
    for SV = 1:31
        SVTxTime = GPSTime_SP3(EpochNumber,2);
        SVWeek = GPSTime_SP3(EpochNumber,1);

        [SV_X_Data_Eph(EpochNumber,SV) SV_Y_Data_Eph(EpochNumber,SV) SV_Z_Data_Eph(EpochNumber,SV) SV_T_Data_Eph(EpochNumber,SV) VehicleIDs_Eph(EpochNumber,SV)] = GPSOrbitPropagator(SVWeek, SVTxTime, SV, NavigationData,21600);
   end
end

% use all sats
for SV = 1:31
    % calculate errors
    X_error(:,SV) = SV_X_Data_Eph(:,SV) - SV_X_Data_SP3(:,SV);
    Y_error(:,SV) = SV_Y_Data_Eph(:,SV) - SV_Y_Data_SP3(:,SV);
    Z_error(:,SV) = SV_Z_Data_Eph(:,SV) - SV_Z_Data_SP3(:,SV);
    T_error(:,SV) = (SV_T_Data_Eph(:,SV) - SV_T_Data_SP3(:,SV)) * c; % convert to metres
end

% plot results
Time = [0:15:15*11];

for SV = 1:31;
    
    if VehicleIDs_Eph(1,SV) == 1
    
        figure;
        subplot(4,1,1),plot(Time,X_error(1:12,SV));
        grid on;
        xlabel('Time Since Ephemeris (mins)');
        ylabel('X Position Error (m)');
        title(sprintf('Predicted to Precise Orbit Determination Errors (SV%d)',SV));

        subplot(4,1,2),plot(Time,Y_error(1:12,SV));
        grid on;
        xlabel('Time Since Ephemeris (mins)');
        ylabel('Y Position Error (m)');

        subplot(4,1,3),plot(Time,Z_error(1:12,SV));
        grid on;
        xlabel('Time Since Ephemeris (mins)');
        ylabel('Z Position Error (m)');

        subplot(4,1,4),plot(Time,T_error(1:12,SV));
        grid on;
        xlabel('Time Since Ephemeris (mins)');
        ylabel('Time Offset Error (m)');
    end
    
end

% calculate the RMS positoin errors


% figure;
% plot(X_error);
% grid on;
% title('Predicted to Precise Orbit Determination Errors');
% xlabel('Time Epoch');
% ylabel('X Position Error (m)');
% 
% figure;
% plot(Y_error);
% grid on;
% title('Predicted to Precise Orbit Determination Errors');
% xlabel('Time Epoch');
% ylabel('Y Position Error (m)');
% 
% figure;
% plot(Z_error);
% grid on;
% title('Predicted to Precise Orbit Determination Errors');
% xlabel('Time Epoch');
% ylabel('Z Position Error (m)');
% 
% figure;
% plot(T_error);
% grid on;
% title('Predicted to Precise Orbit Determination Errors');
% xlabel('Time Epoch');
% ylabel('Time Offset Error (s)');

%plot3(SV_X_Data_SP3(:,1), SV_Y_Data_SP3(:,1), SV_Z_Data_SP3(:,1), 'r*');
%plot3(SV_X_Data_Eph(:),SV_Y_Data_Eph(:),SV_Z_Data_Eph(:), 'r*');

% x data
% figure;
% plot(SV_X_Data_SP3(:,1),'g+');
% hold on;
% plot(SV_X_Data_Eph(:),'r+');
% grid on;
% hold off;
% xlabel('Time Epoch');
% ylabel('X Position (m)');
% 
% % y data
% figure;
% plot(SV_Y_Data_SP3(:,1),'g+');
% hold on;
% plot(SV_Y_Data_Eph(:),'r+');
% grid on;
% hold off;
% xlabel('Time Epoch');
% ylabel('Y Position (m)');
% 
% % z data
% figure;
% plot(SV_Z_Data_SP3(:,1),'g+');
% hold on;
% plot(SV_Z_Data_Eph(:),'r+');
% grid on;
% hold off;
% xlabel('Time Epoch');
% ylabel('Z Position (m)');
% 
% % t data
% figure;
% plot(SV_T_Data_SP3(:,1),'g+');
% hold on;
% plot(SV_T_Data_Eph(:),'r+');
% grid on;
% hold off;
% xlabel('Time Epoch');
% ylabel('T Position (m)');

