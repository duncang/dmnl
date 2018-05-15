% RAIM Availability Simulation - Written by Duncan Greer 17 Nov 2005
%
% This script calculates the HPL_FD and VPL_FD values for a location over a
% 24 hour period
% $Id: GARDSim_HPLVPLSimulation.m 1879 2008-07-15 05:20:21Z n2523710 $


%Define constants as global variables

GPSConstants; %All generic constants put in this script, as global variables;



% set the user locaiton
%Xu = [-5046773.574,2568446.555,-2925288.451];  % S-Block Roof
% Brisbane - S27 23.0, E153 07.1 13' AMSL
%phi_u = (-27 - (23.0 / 60)) * pi / 180;
%lamda_u = (153 + (07.1 / 60)) * pi / 180;
%height_u = 13 / 3.2;  % note.. this should actually be adjusted to be the wgs-84 height, which is different to MSL
%Xu = LLH2ECEF(phi_u, lamda_u, height_u);

% Canberra - S35 18.4, E149 11.7 1886' AMSL
phi_u = (-35 - (18.4 / 60)) * pi / 180;
lamda_u = (149 + (11.7 / 60)) * pi / 180;
height_u = 1886 / 3.2;  % note.. this should actually be adjusted to be the wgs-84 height, which is different to MSL
Xu = LLH2ECEF(phi_u, lamda_u, height_u);

% use precalculated pbias and thresholds from Brown paper
% sigma = 33m, Pfa = 1e-6, Pmd = 1e-3
%thresh=[0 0 0 0 161.423 173.456 182.738 190.650 197.692];
%pbias=[0 0 0 0 263.399 272.950 279.917 285.636 290.606];

% load the a and lamda values from file
NPA20DOFValuesForParityMethod;


% set HALs in metres
HAL_APV = 20;
VAL_APV = 50;

HAL_NPA = 555.6; % m
HAL_TER = 1852; 

% load satellite positions from SP3
disp('Loading SP3 Data');
SP3DataFileName = 'data/sp3/igs13393.sp3'; % data for week 1339, day 3 - 7 sept 05 
if ~exist('NumberSVs_SP3')
    [NumberSVs_SP3, VehicleIDs_SP3, NumberEpochs_SP3, GPSTime_SP3, SV_X_Data_SP3, SV_Y_Data_SP3, SV_Z_Data_SP3, SV_T_Data_SP3] = readsp3(SP3DataFileName);
end

t_SP3 = [0:900:NumberEpochs_SP3*900 - 1];

% set the desired simulation time step  in seconds
Starttime = 0;      % start of data (second of day, UTC)
Endtime = NumberEpochs_SP3 * 900 - 1;    % 24 hours (seconds)
Timestep = 10;       % seconds


t_Sim =[Starttime:Timestep:Endtime];  
NumberEpochs_Sim = length(t_Sim);

% interpolate SP3 data to the time step
disp('Interpolating SP3 Data');
for SV=1:length(VehicleIDs_SP3)
    SV_X_Data_Sim(SV,:) = spline(t_SP3,SV_X_Data_SP3(:,SV), t_Sim);
    SV_Y_Data_Sim(SV,:) = spline(t_SP3,SV_Y_Data_SP3(:,SV), t_Sim);
    SV_Z_Data_Sim(SV,:) = spline(t_SP3,SV_Z_Data_SP3(:,SV), t_Sim);
    SV_T_Data_Sim(SV,:) = spline(t_SP3,SV_T_Data_SP3(:,SV), t_Sim);
end


% load nav data to determine which SVs were available on the day
% TODO - reload nav data file when this one times out - about every 2 hours
NavDataFileName='data/rinex/navfile_7sep05.05N';
[NavData] = freadnav(NavDataFileName);   % field 27 of NavData has the SVHealthy flag


% initialise storage vectors
HPL_H0 = zeros(1,NumberEpochs_Sim);
VPL_H0 = zeros(1,NumberEpochs_Sim);
N = zeros(1,NumberEpochs_Sim);

LineAPV = ones(1,NumberEpochs_Sim) * HAL_APV;
LineNPA = ones(1,NumberEpochs_Sim) * HAL_NPA;
LineTER = ones(1,NumberEpochs_Sim) * HAL_TER;
LineAPV_Vert = ones(1,NumberEpochs_Sim) * VAL_APV;
Avail_APV = zeros(1,NumberEpochs_Sim);
Avail_APV_Vert = zeros(1,NumberEpochs_Sim);
Avail_NPA = zeros(1,NumberEpochs_Sim);
Avail_TER = zeros(1,NumberEpochs_Sim);

disp(sprintf('Beginning Simulation with %d Epochs', NumberEpochs_Sim));

% for each time step
for Epoch = 1:NumberEpochs_Sim
    % for each satellite, calculate the azimuth and elevation.  if hte sv
    % is above the elevation mask (7.5 deg) add it to the current viewable
    % svs.  
    i = 0;  % index of visible satellits
    for SV=1:length(VehicleIDs_SP3)
        
        % index of this SV in the nav data
        SVnavindex = find(NavData(:,1) == SV);
        % if this SV is healthy
        if(NavData(SVnavindex,27) == 0)
            
            % calculate the geometry matrix - partial derivative unit vectors to each
            % SV

            [Az(SV,Epoch), El(SV,Epoch)] = AzEl(Xu, [SV_X_Data_Sim(SV,Epoch),SV_Y_Data_Sim(SV,Epoch),SV_Z_Data_Sim(SV,Epoch)]);
            if(El(SV,Epoch) > (7.5 * pi / 180))
                i = i+1;
                G(i,1) = cos(El(SV,Epoch)) * cos(Az(SV,Epoch));
                G(i,2) = cos(El(SV,Epoch)) * sin(Az(SV,Epoch));
                G(i,3) = sin(El(SV,Epoch));
                G(i,4) = 1;
            end
                  % limit number of measurements to 9
            if(i == 9)
                break;
            end
        else
            %disp('SV Unhealthy');
        end
    end
    
    % save number of SVs visible this epoch
    N(Epoch) = i;
    
    % find the expected variance of satellite measurements
    SigmaS = 2.0;
    var_i = SigmaS ^ 2;

    % find threshold and pbias values
         
    %Index to the a and lambda array
    a_ind = i-4;

     %-------------------------------
     %determine threshold value
     %------------------------------ 
    %normalised Td
    Td_norm = sqrt(a(a_ind));
    %Unnormalised Td
    Td = SigmaS*Td_norm; %metres

     %-------------------------
     %calculate pbias values
     %-------------------------
     %normalised pbias 

     pbias_norm = sqrt(lambda(a_ind));
     pbias = pbias_norm*SigmaS;

    
    % form weighting matrix
    W_inv = eye(N(Epoch)) * var_i;
    W = inv(W_inv);

    % calculate hte partiy transform matrix, P
    [Q,R] = qr(G);
    QT = Q';
    P = QT(5:i,:);
    
    A = inv(G' * G) * G';
    %B = G * A;
    S = eye(i) - G*A;
    
    SLOPE = zeros(1,i);
    
    for(SV = 1:i)
        SLOPE(SV) = sqrt(A(1,SV)^2 + A(2,SV)^2) / sqrt(S(SV,SV));
    end
    
    % sort the slopes
    [Y,I] = sort(SLOPE);
    SLOPEMAX(Epoch) = Y(i);
    biassat = I(i);
    
    rbias = pbias / norm(P(:,biassat));
    epsilon = zeros(i,1);
    epsilon(biassat) = rbias;
    poserr = A * epsilon;
    
    HPL(Epoch) = norm([poserr(1) poserr(2)]);
    VPL(Epoch) = norm(poserr(3));
    
    % if the HPL is above the HAL, mark this epoch as unavailable
    if(HPL(Epoch) > LineNPA(Epoch))
        Avail_NPA(Epoch) = 0;
    else
        Avail_NPA(Epoch) = 1;
    end
    
    if(HPL(Epoch) > LineTER(Epoch))
        Avail_TER(Epoch) = 0;
    else
        Avail_TER(Epoch) = 1;
    end
    
    if(HPL(Epoch) > LineAPV(Epoch))
        Avail_APV(Epoch) = 0;
    else
        Avail_APV(Epoch) = 1;
    end
    
    if(VPL(Epoch) > VAL_APV)
        Avail_APV_Vert(Epoch) = 0;
    else
        Avail_APV_Vert(Epoch) = 1;
    end
    
    
    %%% CODE FROM GRAS MOPS COMMENTED OUT %%%
%     % calculate S matrix
%     S_2 = inv(G' * W * G) * G' * W;
%     
%     % calculate s and d values (e,n,u,t)
%     for i=1:N(Epoch)
%         s_east_dot_var(i) = S_2(1,i)^2 * var_i;
%         s_north_dot_var(i) = S_2(2,i)^2 * var_i;
%         s_up_dot_var(i) = S_2(3,i)^2 * var_i;
%         s_t_dot_var(i) = S_2(4,i)^2 * var_i;
%         s_east_north_dot_var(i) = S_2(1,i) * S_2(2,i) * var_i;
%     end
%     
%     d_east_2 = sum(s_east_dot_var);
%     d_north_2 = sum(s_north_dot_var);
%     d_up_2 = sum(s_up_dot_var);
%     d_t_2 = sum(s_t_dot_var);
%     d_EN_2 = sum(s_east_north_dot_var);
%     
%     d_major = sqrt(((d_east_2 + d_north_2) / 2) + sqrt(((d_east_2 - d_north_2) / 2) + (d_EN_2)));
%     
%     % calculate HPL and VPL values - see GRAS Mops Section 2.3.10.1
%     K_ffmd = 5.8; % see GRAS Mops Section 2.3.11.5.2.1.2
%     HPL_H0(Epoch) = 10 * d_major;
%     VPL_H0(Epoch) = K_ffmd * sqrt(d_up_2);
   


    % clear matrix variables
    clear A B Q R Y I P W S W_inv G s_east_dot_var s_north_dot_var s_up_dot_var s_t_dot_var s_east_north_dot_var 
    
    % display status message
    if(mod(Epoch,100) == 0)
        disp(sprintf('Epoch: %d',Epoch));
    end
end

% find hte availability of NPA for this day
Avail_NPA_PC = (sum(Avail_NPA) / NumberEpochs_Sim) * 100
Avail_TER_PC = (sum(Avail_TER) / NumberEpochs_Sim) * 100
Avail_APV_PC = (sum(Avail_APV) / NumberEpochs_Sim) * 100
Avail_APV_Vert_PC =  (sum(Avail_APV_Vert) / NumberEpochs_Sim) * 100

% plot results

figure();
hold on;
plot(HPL,'b');
%plot(abs(HPL_H0),'r')
plot(LineAPV,'r-');
plot(LineNPA,'g-');
plot(LineTER,'y-');
grid on;
xlabel('Time of Day (s)');
ylabel('Metres');
legend('HPL FD','HAL APV','HAL NPA','HAL TER');
title(sprintf('HPL for Canberra 7 Sept 2005 with sigma=%2.1dm, Pfa = 2.22e-7, Pmd = 1e-3',SigmaS));
hold off;

figure();
hold on;
plot(VPL,'b');
plot(LineAPV_Vert,'r-');
grid on;
xlabel('Time of Day (s)');
ylabel('Metres');
legend('VPL FD','VAL APV');
title(sprintf('VPL for Canberra 7 Sept 2005 with sigma=%2.1dm, Pfa = 2.22e-7, Pmd = 1e-3',SigmaS));
hold off;




