
%7.2.06

feature accel on


%The data is in this order:


%
%                   Station Locations:
% Alice Springs
% Ceduna
% Darwin

% Hillarys
% Hobart
% Jabiru
% Karratha
% Melbourne
% Mount Stromlo
% Parkes
% Sydney
% Tidbinbilla
% Townsville
% Yaragadee
%
%

Filename_1 = 'data\argn\alic0330.07o' ;      % Alice Springs
Filename_2 = 'data\argn\cedu0330.07o';       % Ceduna
Filename_3 = 'data\argn\darw0330.07o';      %Darwin
Filename_4 = 'data\argn\hob20330.07o';     % Hobart
Filename_5 = 'data\argn\jab10330.07o';     % Jabiru
Filename_6 = 'data\argn\karr0330.07o';      % Karratha
Filename_7 = 'data\argn\mobs0330.07o';     % Melbourne
Filename_8 = 'data\argn\sydn0330.07o';  %Sydney
Filename_9 = 'data\argn\tow20330.07o';   %Townsville




[GPSTime_Week_1, GPSTime_Sec_1,NumberRinexObsTypes_1,ApproxPos_1, DATASTRUCT_1] = ReadRinexGRS(Filename_1);
[GPSTime_Week_2, GPSTime_Sec_2,NumberRinexObsTypes_2,ApproxPos_2, DATASTRUCT_2] = ReadRinexGRS(Filename_2);
[GPSTime_Week_3, GPSTime_Sec_3,NumberRinexObsTypes_3,ApproxPos_3, DATASTRUCT_3] = ReadRinexGRS(Filename_3);
[GPSTime_Week_4, GPSTime_Sec_4,NumberRinexObsTypes_4,ApproxPos_4, DATASTRUCT_4] = ReadRinexGRS(Filename_4);
[GPSTime_Week_5, GPSTime_Sec_5,NumberRinexObsTypes_5,ApproxPos_5, DATASTRUCT_5] = ReadRinexGRS(Filename_5);
[GPSTime_Week_6, GPSTime_Sec_6,NumberRinexObsTypes_6,ApproxPos_6, DATASTRUCT_6] = ReadRinexGRS(Filename_6);
[GPSTime_Week_7, GPSTime_Sec_7,NumberRinexObsTypes_7,ApproxPos_7, DATASTRUCT_7] = ReadRinexGRS(Filename_7);
[GPSTime_Week_8, GPSTime_Sec_8,NumberRinexObsTypes_8,ApproxPos_8, DATASTRUCT_8] = ReadRinexGRS(Filename_8);
[GPSTime_Week_9, GPSTime_Sec_9,NumberRinexObsTypes_9,ApproxPos_9, DATASTRUCT_9] = ReadRinexGRS(Filename_9);




%    Filename_10 = 'data\argn\mobs0330.07o';
%     Filename_11 = 'data\argn\tow20330.07o';
%     Filename_12 = 'data\argn\tow20330.07o';
%       Filename_13 = 'data\argn\alic0330.07o' ;
%     Filename_14 = 'data\argn\mobs0330.07o';
%     Filename_15 = 'data\argn\tow20330.07o';
%     Filename_16 = 'data\argn\tow20330.07o';


%
%
% coco0330.07o
% darw0330.07o
% dav10330.07o
% jab10330.07o
% mac10330.07o
% mobs0330.07o
% str10330.07o
% tid10330.07o
% yar20330.07o
% tow20330.07o
% maw10330.07o
% sydn0330.07o
% karr0330.07o
% hob20330.07o
% cas10330.07o




%note that this data is at 10 Hz!!!
Filename_GVS = 'data\Ground_Test_Data\2Feb2007\02020330.07O';

[GPSTime_Week_GVS, GPSTime_Sec_GVS,NumberRinexObsTypes_GVS,ApproxPos_GVS, DATASTRUCT_GVS] = ReadRinexGRS(Filename_GVS);





%read NAV file

%brdc0330.07n

%Calculate Satellite Positions based on GVS data:

Filename_Nav = 'data\Ground_Test_Data\2Feb2007\02020330.07N';

%read optimised constellation data from file
[NavData, ALPHA,BETA] = freadNav(Filename_Nav);

MaxGVS = length(GPSTime_Sec_GVS);
NumSatEpochs = MaxGVS-5;

%Pre-allocate
SV_X_Data = zeros(NumSatEpochs,32);
SV_Y_Data = zeros(NumSatEpochs,32);
SV_Z_Data = zeros(NumSatEpochs,32);
SV_T_Data = zeros(NumSatEpochs,32);
ValidData_Satellite = zeros(NumSatEpochs,32);

SV_Xvel_Data = zeros(NumSatEpochs,32);
SV_Yvel_Data = zeros(NumSatEpochs,32);
SV_Zvel_Data = zeros(NumSatEpochs,32);
SV_Tvel_Data = zeros(NumSatEpochs,32);
SV_Xacc_Data = zeros(NumSatEpochs,32);
SV_Yacc_Data = zeros(NumSatEpochs,32);
SV_Zacc_Data = zeros(NumSatEpochs,32);
SV_Tacc_Data = zeros(NumSatEpochs,32);
ValidDataSatVels = zeros(NumSatEpochs,32);



max_temp = size(DATASTRUCT_GVS1.C1);
max_number_of_sats = max_temp(1);
GPSConstants;

for i = 1:NumSatEpochs  %epoch counter in seconds , number of points to calculate satellite positions for.



    SVWeek(i) = GPSTime_Week_GVS(i);    
    

    for SV = 1:max_number_of_sats

        PR = DATASTRUCT_GVS1.C1(SV,i);

        SVTime_GPSSecs(i) = GPSTime_Sec_GVS1(i) - PR/c; %pass the time of transmission to the GPSOrbitPropagator


        [SV_X_Data(i,SV) SV_Y_Data(i,SV) SV_Z_Data(i,SV) SV_T_Data(i,SV) ValidData_Satellite(i,SV)] = GPSOrbitPropagator(SVWeek(i), SVTime_GPSSecs(i), SV, NavData);



    end
end








%====================================
% GRS FUNCTION
%====================================

%for each station calculate the ionospheric delay


%Do least squares solution , note that for the GVS's it is every 30
%seconds.


%read NAV file



Filename_Nav1 = 'data\Ground_Test_Data\2Feb2007\brdc0330.07n';

%read optimised constellation data from file
[NavData_1, ALPHA_1,BETA_1] = freadNav(Filename_Nav1);

Max_1 = length(GPSTime_Sec_1);
NumSatEpochs_1 = Max_1-5;

%Pre-allocate
SV_X_Data_1 = zeros(NumSatEpochs_1,32);
SV_Y_Data_1 = zeros(NumSatEpochs_1,32);
SV_Z_Data_1 = zeros(NumSatEpochs_1,32);
SV_T_Data_1 = zeros(NumSatEpochs_1,32);
ValidData_Satellite_1 = zeros(NumSatEpochs_1,32);




max_temp_1 = size(DATASTRUCT_1.C1);
max_number_of_sats_1 = max_temp_1(1);


for i = 1:NumSatEpochs_1  %epoch counter in seconds , number of points to calculate satellite positions for.



    SVWeek_1(i) = GPSTime_Week_1(i);
    
    

    for SV = 1:max_number_of_sats_1

        PR = DATASTRUCT_1.C1(SV,i);

        SVTime_GPSSecs_1(i) = GPSTime_Sec_1(i) - PR/c; %pass the time of transmission to the GPSOrbitPropagator


        [SV_X_Data_1(i,SV) SV_Y_Data_1(i,SV) SV_Z_Data_1(i,SV) SV_T_Data_1(i,SV) ValidData_Satellite_1(i,SV)] = GPSOrbitPropagator(SVWeek_1(i), SVTime_GPSSecs_1(i), SV, NavData_1);



    end
end



%ALICE
%initial starting
UserPos_1 = ApproxPos_1;
    UserPos_1(4) = 0; 
    
for i = 1:length(GPSTime_Sec_1)-5

    k = 1;
    for SV = 1:32

        if  (ValidData_Satellite_1(i,SV) == 1) && (DATASTRUCT_1.C1(SV,i) ~= 0)

            SVPos_1(k,1) = SV_X_Data_1(i,SV);
            SVPos_1(k,2) = SV_Y_Data_1(i,SV);
            SVPos_1(k,3) = SV_Z_Data_1(i,SV);
            SVPos_1(k,4) = c*SV_T_Data_1(i,SV);
            
            
            PRMeasured_Observed_1(k) = DATASTRUCT_1.C1(SV,i);
            
            k = k+1;
            
            
        end         
        
        
    end

    
    N_1 = length(PRMeasured_Observed_1);
    
   % UserPos_1 = ApproxPos_1;
    %UserPos_1(4) = 0;  %assume clock bias is 0 to start with
    
%fix the x y z at the approx position to estimate the clock bias and other
%terms

 [SolutionVec_Observed_1(i,1:4) VarSolutionVec_Observed_1 NumIterations_Observed_1(i) ResVec_Observed_1 M_Observed_1 LSQ_Fail_Observed_1(i) limit_Observed_1(i) DOP_Observed_1(i,1:5)] = GARD_LSQ(UserPos_1,N_1,PRMeasured_Observed_1,SVPos_1);

 %use this one if you only want to estimate the clock bias (assume you know
 %the position exactly of the receiver.
 %[SolutionVec_Observed_1(i,1:4) VarSolutionVec_Observed_1 NumIterations_Observed_1(i) ResVec_Observed_1 M_Observed_1 LSQ_Fail_Observed_1(i) limit_Observed_1(i) DOP_Observed_1(i,1:5)] = GARD_LSQ_ClockBias(UserPos_1,N_1,PRMeasured_Observed_1,SVPos_1);
  
   UserPos_1(1:3) = ApproxPos_1;
  UserPos_1(4) = SolutionVec_Observed_1(i,4) ;
 
end


%DARWIN

UserPos_3 = ApproxPos_3;
    UserPos_3(4) = 0; 
    
for i = 1:length(GPSTime_Sec_3)-5

    k = 1;
    for SV = 1:32

        if  (ValidData_Satellite_1(i,SV) == 1) && (DATASTRUCT_3.C1(SV,i) ~= 0)

            SVPos_3(k,1) = SV_X_Data_1(i,SV);
            SVPos_3(k,2) = SV_Y_Data_1(i,SV);
            SVPos_3(k,3) = SV_Z_Data_1(i,SV);
            SVPos_3(k,4) = c*SV_T_Data_1(i,SV);
            
            
            PRMeasured_Observed_3(k) = DATASTRUCT_3.C1(SV,i);
            
            k = k+1;
            
            
        end         
        
        
    end

    
    N_3 = length(PRMeasured_Observed_3);
    
   % UserPos_1 = ApproxPos_1;
    %UserPos_1(4) = 0;  %assume clock bias is 0 to start with
    
%fix the x y z at the approx position to estimate the clock bias and other
%terms

 [SolutionVec_Observed_3(i,1:4) VarSolutionVec_Observed_3 NumIterations_Observed_3(i) ResVec_Observed_3 M_Observed_3 LSQ_Fail_Observed_3(i) limit_Observed_3(i) DOP_Observed_3(i,1:5)] = GARD_LSQ(UserPos_3,N_3,PRMeasured_Observed_3,SVPos_3);

 %use this one if you only want to estimate the clock bias (assume you know
 %the position exactly of the receiver.
 %[SolutionVec_Observed_1(i,1:4) VarSolutionVec_Observed_1 NumIterations_Observed_1(i) ResVec_Observed_1 M_Observed_1 LSQ_Fail_Observed_1(i) limit_Observed_1(i) DOP_Observed_1(i,1:5)] = GARD_LSQ_ClockBias(UserPos_1,N_1,PRMeasured_Observed_1,SVPos_1);
  
   UserPos_3(1:3) = ApproxPos_3;
  UserPos_3(4) = SolutionVec_Observed_3(i,4) ;
 
end

 
%TOWNSVILLE

UserPos_9 = ApproxPos_9;
    UserPos_9(4) = 0; 
    
for i = 1:length(GPSTime_Sec_9)-5

    k = 1;
    for SV = 1:32

        if  (ValidData_Satellite_1(i,SV) == 1) && (DATASTRUCT_9.C1(SV,i) ~= 0)

            SVPos_9(k,1) = SV_X_Data_1(i,SV);
            SVPos_9(k,2) = SV_Y_Data_1(i,SV);
            SVPos_9(k,3) = SV_Z_Data_1(i,SV);
            SVPos_9(k,4) = c*SV_T_Data_1(i,SV);
            
            
            PRMeasured_Observed_9(k) = DATASTRUCT_9.C1(SV,i);
            
            k = k+1;
            
            
        end         
        
        
    end

    
    N_9 = length(PRMeasured_Observed_9);
    
   % UserPos_1 = ApproxPos_1;
    %UserPos_1(4) = 0;  %assume clock bias is 0 to start with
    
%fix the x y z at the approx position to estimate the clock bias and other
%terms

 [SolutionVec_Observed_9(i,1:4) VarSolutionVec_Observed_9 NumIterations_Observed_9(i) ResVec_Observed_9 M_Observed_9 LSQ_Fail_Observed_9(i) limit_Observed_9(i) DOP_Observed_9(i,1:5)] = GARD_LSQ(UserPos_9,N_9,PRMeasured_Observed_9,SVPos_9);

 %use this one if you only want to estimate the clock bias (assume you know
 %the position exactly of the receiver.
 %[SolutionVec_Observed_1(i,1:4) VarSolutionVec_Observed_1 NumIterations_Observed_1(i) ResVec_Observed_1 M_Observed_1 LSQ_Fail_Observed_1(i) limit_Observed_1(i) DOP_Observed_1(i,1:5)] = GARD_LSQ_ClockBias(UserPos_1,N_1,PRMeasured_Observed_1,SVPos_1);
  
   UserPos_9(1:3) = ApproxPos_9;
  UserPos_9(4) = SolutionVec_Observed_9(i,4) ;
 
end

 

%SYDNEY

UserPos_8 = ApproxPos_8;
    UserPos_8(4) = 0; 
    
for i = 1:length(GPSTime_Sec_8)-5

    k = 1;
    for SV = 1:32

        if  (ValidData_Satellite_1(i,SV) == 1) && (DATASTRUCT_8.C1(SV,i) ~= 0)

            SVPos_8(k,1) = SV_X_Data_1(i,SV);
            SVPos_8(k,2) = SV_Y_Data_1(i,SV);
            SVPos_8(k,3) = SV_Z_Data_1(i,SV);
            SVPos_8(k,4) = c*SV_T_Data_1(i,SV);
            
            
            PRMeasured_Observed_8(k) = DATASTRUCT_8.C1(SV,i);
            
            k = k+1;
            
            
        end         
        
        
    end

    
    N_8 = length(PRMeasured_Observed_8);
    
   % UserPos_1 = ApproxPos_1;
    %UserPos_1(4) = 0;  %assume clock bias is 0 to start with
    
%fix the x y z at the approx position to estimate the clock bias and other
%terms

 [SolutionVec_Observed_8(i,1:4) VarSolutionVec_Observed_8 NumIterations_Observed_8(i) ResVec_Observed_8 M_Observed_8 LSQ_Fail_Observed_8(i) limit_Observed_8(i) DOP_Observed_8(i,1:5)] = GARD_LSQ(UserPos_8,N_8,PRMeasured_Observed_8,SVPos_8);

 %use this one if you only want to estimate the clock bias (assume you know
 %the position exactly of the receiver.
 %[SolutionVec_Observed_1(i,1:4) VarSolutionVec_Observed_1 NumIterations_Observed_1(i) ResVec_Observed_1 M_Observed_1 LSQ_Fail_Observed_1(i) limit_Observed_1(i) DOP_Observed_1(i,1:5)] = GARD_LSQ_ClockBias(UserPos_1,N_1,PRMeasured_Observed_1,SVPos_1);
  
   UserPos_8(1:3) = ApproxPos_8;
  UserPos_8(4) = SolutionVec_Observed_8(i,4) ;
 
end

 

%MELBOURNE

UserPos_7 = ApproxPos_7;
    UserPos_7(4) = 0; 
    
for i = 1:length(GPSTime_Sec_7)-5

    k = 1;
    for SV = 1:32

        if  (ValidData_Satellite_1(i,SV) == 1) && (DATASTRUCT_7.C1(SV,i) ~= 0)

            SVPos_7(k,1) = SV_X_Data_1(i,SV);
            SVPos_7(k,2) = SV_Y_Data_1(i,SV);
            SVPos_7(k,3) = SV_Z_Data_1(i,SV);
            SVPos_7(k,4) = c*SV_T_Data_1(i,SV);
            
            
            PRMeasured_Observed_7(k) = DATASTRUCT_7.C1(SV,i);
            
            k = k+1;
            
            
        end         
        
        
    end

    
    N_7 = length(PRMeasured_Observed_7);
    
   % UserPos_1 = ApproxPos_1;
    %UserPos_1(4) = 0;  %assume clock bias is 0 to start with
    
%fix the x y z at the approx position to estimate the clock bias and other
%terms

 [SolutionVec_Observed_7(i,1:4) VarSolutionVec_Observed_7 NumIterations_Observed_7(i) ResVec_Observed_7 M_Observed_7 LSQ_Fail_Observed_7(i) limit_Observed_7(i) DOP_Observed_7(i,1:5)] = GARD_LSQ(UserPos_7,N_7,PRMeasured_Observed_7,SVPos_7);

 %use this one if you only want to estimate the clock bias (assume you know
 %the position exactly of the receiver.
 %[SolutionVec_Observed_1(i,1:4) VarSolutionVec_Observed_1 NumIterations_Observed_1(i) ResVec_Observed_1 M_Observed_1 LSQ_Fail_Observed_1(i) limit_Observed_1(i) DOP_Observed_1(i,1:5)] = GARD_LSQ_ClockBias(UserPos_1,N_1,PRMeasured_Observed_1,SVPos_1);
  
   UserPos_7(1:3) = ApproxPos_7;
  UserPos_7(4) = SolutionVec_Observed_7(i,4) ;
 
end

 
 

%THERES SOMETHING WRONG WITH SYDNEY DATA, POSITOIN SOLUTION HAD LARGE
%ERRORS



%Simulated GMS function


%Using position data of 4 GVS's, determine the position of the satellites. 



%for alice springs:

max_length = length(GPSTime_Week_1);

for i = 1:

    
    i = 1
    
 N is number of ground stations
  for each satellite do this
     
      
      
      
 N = 4
 

     
     GRSPos(1,1:3) = ApproxPos_1; %alice
     GRSPos(2,1:3) = ApproxPos_3; %darwin
     GRSPos(3,1:3) = ApproxPos_9; %townsville
     GRSPos(4,1:3) = ApproxPos_7; %melbourne
     
     %the receiver clock biases:
     GRSPos(1,4) = SolutionVec_Observed_1(i,4)
       GRSPos(2,4) = SolutionVec_Observed_3(i,4)
         GRSPos(3,4) = SolutionVec_Observed_9(i,4)
           GRSPos(4,4) = SolutionVec_Observed_7(i,4)
     
     
     %for SV 11 
     
     PRMeasured(1) = DATASTRUCT_1.C1(11,i);
       PRMeasured(2) = DATASTRUCT_3.C1(11,i);
        PRMeasured(3) = DATASTRUCT_9.C1(11,i);
          PRMeasured(4) = DATASTRUCT_7.C1(11,i);
           
          
          %GPSTime_Sec_1(1) 
          
     
          SVPosStart(1) = SV_X_Data_1(i,11) ;
          SVPosStart(2) = SV_Y_Data_1(i,11);
          SVPosStart(3) = SV_Z_Data_1(i,11) ;
          SVPosStart(4) = c*SV_T_Data_1(i,11);
          
     
       

[SolutionVec, VarSolutionVec, NumIterations, ResidualVector, M, LSQ_Fail, limit, DOP] = GARD_LSQ_GMS(GRSPos,N,PRMeasured,SVPosStart);

















%Simulate GVS function


















































