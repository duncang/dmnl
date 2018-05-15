
%Main function for GARD Sim *outside GUI



%Version 1.00
%Troy Bruggemann  11 August 2005

%===============================================================
%CONSTANTS
%All constants used should be in this section
%===============================================================

%hold;

%Define constants as global variables

GPSConstants; %All generic constants put in this script, as global variables;


 
%===============================================================
%===============================================================
%PLANNING MODULE
%This module takes inputs from the user, reads in data from files
%===============================================================
%===============================================================

%the number of satellite positions to calculate for


%turn JIT accelerator on (speeds up code processing)
 
%feature accel on %only works in matlab 7

%always make NumSatEpochs greater than NumEpochs
NumSatEpochs = 300;


%-----------------------------------------------------------------
%READ IN DATA FROM FILES
%-----------------------------------------------------------------

%READ RINEX
%Data from Ashtech receiver:

%Observation File
%Read      ashtech micro z 24 hour

%generic filename and path:


filepathObs = 'data\AshtechuZ\EESE2931_mod.05O';



%for Ashtech uZ-CGRS 24 hours test
%filepathObs = 'data\24hrAshtechMicro29.11\EESE334*_mod.04O';

%Check if files exist
IsFilesObs = dir(filepathObs);  %this returns results in order of what the OS gives and it seems to give it in order of filename alphabetically

%check the time returned in gpstime to see if theyre in order or not*=

%Get the total number of files existing
NumFilesObs = size(IsFilesObs);

if NumFilesObs(1) == 0
    error('No such files existing')
else
%Get the file info
[pathstrObs,nameObs,extObs,versnObs] = fileparts(filepathObs) 



%initialise vectors and arrays

%initilaising with zeros makes it run a lot quicker



for filecounterObs = 1:NumFilesObs(1)

%clear temp variables    
clear GPStimeAMicroTemp ValidDataRinexObsTemp L1_PRNAMicroTemp L2_PRNAMicroTemp C1_PRNAMicroTemp P1_PRNAMicroTemp P2_PRNAMicroTemp;

%File to read in  
FilenameTempObs = IsFilesObs(filecounterObs).name;

FilenameObs = fullfile(pathstrObs,FilenameTempObs)

[GPStimeAMicroWeekTemp GPStimeAMicroTemp, ValidDataRinexObsTemp, L1_PRNAMicroTemp, L2_PRNAMicroTemp, C1_PRNAMicroTemp, P1_PRNAMicroTemp, P2_PRNAMicroTemp] = ReadRinexSblockRoof(FilenameObs);




%initialise 

if (filecounterObs == 1)
    
%add the temp data to the data to be processed
%initialise with zeros
GPStimeAMicroWeek = zeros(size(GPStimeAMicroWeekTemp));  %GPS week
GPStimeAMicro = zeros(size(GPStimeAMicroTemp));          % GPS seconds into week
ValidDataRinexObs = zeros(size(ValidDataRinexObsTemp));
L1_PRNAMicro  = zeros(size(L1_PRNAMicroTemp));
L2_PRNAMicro  = zeros(size(L2_PRNAMicroTemp));
C1_PRNAMicro = zeros(size(C1_PRNAMicroTemp));
P1_PRNAMicro = zeros(size(P1_PRNAMicroTemp));
P2_PRNAMicro = zeros(size(P2_PRNAMicroTemp));    
    
GPStimeAMicroWeek = GPStimeAMicroWeekTemp;     
GPStimeAMicro = GPStimeAMicroTemp;
ValidDataRinexObs = ValidDataRinexObsTemp;
L1_PRNAMicro = L1_PRNAMicroTemp;
L2_PRNAMicro = L2_PRNAMicroTemp;
C1_PRNAMicro = C1_PRNAMicroTemp;
P1_PRNAMicro = P1_PRNAMicroTemp;
P2_PRNAMicro = P2_PRNAMicroTemp;
end

if (filecounterObs ~= 1)

    %cocatenate the previous set of data with the next set
    
GPStimeAMicroWeek = cat(2, GPStimeAMicroWeek, GPStimeAMicroWeekTemp);    
GPStimeAMicro = cat(2, GPStimeAMicro, GPStimeAMicroTemp);
ValidDataRinexObs = cat(2, ValidDataRinexObs, ValidDataRinexObsTemp);
L1_PRNAMicro  = cat(2, L1_PRNAMicro, L1_PRNAMicroTemp);
L2_PRNAMicro  = cat(2, L2_PRNAMicro, L2_PRNAMicroTemp);
C1_PRNAMicro  = cat(2, C1_PRNAMicro, C1_PRNAMicroTemp);
P1_PRNAMicro  = cat(2, P1_PRNAMicro, P1_PRNAMicroTemp);
P2_PRNAMicro  = cat(2, P2_PRNAMicro, P2_PRNAMicroTemp);

end %if (filecounterObs ~= 1)


end; %end for filecounterObs = 1:NumFilesObs(1)


end; %end else



%----------------------------------------------------------------------


%READ RINEX Navigation File
%Data from Ashtech receiver:


%generic filename and path:


%for repeater Test

filepathNav = 'data\AshtechuZ\EESE2931_mod.05N';


%for Ashtech uZ-CGRS 24 hours test
%filepathNav = 'data\24hrAshtechMicro29.11\EESE334*_mod.04N';

%Check if files exist
IsFilesNav = dir(filepathNav);

%Get the total number of files existing
NumFilesNav = size(IsFilesNav);

if NumFilesNav(1) == 0
    error('No such files existing')
else
%Get the file info
[pathstrNav,nameNav,extNav,versnNav] = fileparts(filepathNav); 


%initialise vectors and arrays

%initilaising with zeros makes it run a lot quicker

%initialise

for filecounterNav = 1:NumFilesNav(1)

%clear temp variables    
clear NavigationDataTemp;

%File to read in  
FilenameTempNav = IsFilesNav(filecounterNav).name;
FilenameNav = fullfile(pathstrNav,FilenameTempNav)
NavigationDataTemp = freadnav(FilenameNav);

%add the temp data to the data to be processed

%initialise 

if (filecounterNav == 1)
    
NavData = zeros(size(NavigationDataTemp))
NavData =  NavigationDataTemp;

end

if (filecounterNav ~= 1)
    %cocatenate the previous set of data with the next set
NavData  = cat(1, NavData,NavigationDataTemp); %this is the matrix with all the nav data in it

end %if (filecounterNav ~= 1)


end; %end for filecounterNav = 1:NumFilesNav(1)


end; %end else


% 
% 
% %READ RINEX for single frequency measurements
% %Data for 835 Assignment
% 
% %Observation File
% %Read Orion data
% 
% %generic filename and path:
% 
% %filepathObs = 'data\835Assignment\Orion835.O';
% 
% 
% filepathObs = 'data\OrionFiles\Arch.O';
% 
% 
%  %for repeater Test
% %filepathObs = 'data\RepeaterTest8Sept\RoomSess1.05O';
% 
% 
% %Check if files exist
% IsFilesObs = dir(filepathObs);  %this returns results in order of what the OS gives and it seems to give it in order of filename alphabetically
% 
% %check the time returned in gpstime to see if theyre in order or not*=
% 
% %Get the total number of files existing
% NumFilesObs = size(IsFilesObs);
% 
% if NumFilesObs(1) == 0
%     error('No such files existing')
% else
% %Get the file info
% [pathstrObs,nameObs,extObs,versnObs] = fileparts(filepathObs) 
% 
% 
% 
% %initialise vectors and arrays
% 
% %initilaising with zeros makes it run a lot quicker
% 
% 
% 
% for filecounterObs = 1:NumFilesObs(1)
% 
% %clear temp variables    
% clear GPStimeAMicroTemp ValidDataRinexObsTemp L1_PRNAMicroTemp L2_PRNAMicroTemp C1_PRNAMicroTemp P1_PRNAMicroTemp P2_PRNAMicroTemp;
% 
% %File to read in  
% FilenameTempObs = IsFilesObs(filecounterObs).name;
% 
% FilenameObs = fullfile(pathstrObs,FilenameTempObs)
% 
% %[GPStimeAMicroWeekTemp GPStimeAMicroTemp, ValidDataRinexObsTemp, L1_PRNAMicroTemp, L2_PRNAMicroTemp, C1_PRNAMicroTemp, P1_PRNAMicroTemp, P2_PRNAMicroTemp] = ReadRinexSblockRoof(FilenameObs);
% [GPStimeAMicroWeekTemp,GPStimeAMicroTemp,ValidDataRinexObsTemp, L1_PRNAMicroTemp, C1_PRNAMicroTemp] = ReadRinexObsSuperStar(FilenameObs);
% 
% 
% 
% %initialise 
% 
% if (filecounterObs == 1)
%     
% %add the temp data to the data to be processed
% %initialise with zeros
% GPStimeAMicroWeek = zeros(size(GPStimeAMicroWeekTemp));  %GPS week
% GPStimeAMicro = zeros(size(GPStimeAMicroTemp));          % GPS seconds into week
% ValidDataRinexObs = zeros(size(ValidDataRinexObsTemp));
% L1_PRNAMicro  = zeros(size(L1_PRNAMicroTemp));
% 
% C1_PRNAMicro = zeros(size(C1_PRNAMicroTemp));
%  
%     
% GPStimeAMicroWeek = GPStimeAMicroWeekTemp;     
% GPStimeAMicro = GPStimeAMicroTemp;
% ValidDataRinexObs = ValidDataRinexObsTemp;
% L1_PRNAMicro = L1_PRNAMicroTemp;
% 
% C1_PRNAMicro = C1_PRNAMicroTemp;
% end
% 
% if (filecounterObs ~= 1)
% 
%     %cocatenate the previous set of data with the next set
%     
% GPStimeAMicroWeek = cat(2, GPStimeAMicroWeek, GPStimeAMicroWeekTemp);    
% GPStimeAMicro = cat(2, GPStimeAMicro, GPStimeAMicroTemp);
% ValidDataRinexObs = cat(2, ValidDataRinexObs, ValidDataRinexObsTemp);
% L1_PRNAMicro  = cat(2, L1_PRNAMicro, L1_PRNAMicroTemp);
% C1_PRNAMicro  = cat(2, C1_PRNAMicro, C1_PRNAMicroTemp);
% 
% end %if (filecounterObs ~= 1)
% 
% 
% end; %end for filecounterObs = 1:NumFilesObs(1)
% 
% 
% end; %end else
% 
% 
% 
% %----------------------------------------------------------------------
% 
% 
% %READ RINEX Navigation File
% %Data from Ashtech receiver:
% 
% 
% %generic filename and path:
% 
% %filepathNav = 'data\835Assignment\Orion835.N';
% 
% 
% filepathNav = 'data\OrionFiles\Arch.N';
% 
% %for repeater test
% %filepathNav = 'data\RepeaterTest8Sept\RoofSess1.05N';
% 
% 
% %Check if files exist
% IsFilesNav = dir(filepathNav);
% 
% %Get the total number of files existing
% NumFilesNav = size(IsFilesNav);
% 
% if NumFilesNav(1) == 0
%     error('No such files existing')
% else
% %Get the file info
% [pathstrNav,nameNav,extNav,versnNav] = fileparts(filepathNav); 
% 
% 
% %initialise vectors and arrays
% 
% %initilaising with zeros makes it run a lot quicker
% 
% %initialise
% 
% for filecounterNav = 1:NumFilesNav(1)
% 
% %clear temp variables    
% clear NavigationDataTemp;
% 
% %File to read in  
% FilenameTempNav = IsFilesNav(filecounterNav).name;
% FilenameNav = fullfile(pathstrNav,FilenameTempNav)
% NavigationDataTemp = freadnav(FilenameNav);
% 
% %add the temp data to the data to be processed
% 
% %initialise 
% 
% if (filecounterNav == 1)
%     
% NavData = zeros(size(NavigationDataTemp))
% NavData =  NavigationDataTemp;
% 
% end
% 
% if (filecounterNav ~= 1)
%     %cocatenate the previous set of data with the next set
% NavData  = cat(1, NavData,NavigationDataTemp); %this is the matrix with all the nav data in it
% 
% end %if (filecounterNav ~= 1)
% 
% 
% end; %end for filecounterNav = 1:NumFilesNav(1)
% 
% 
% end; %end else
% 















%=============================================================
%OBTAIN SATELLITE POSITIONS (x,y,z,T)
%=============================================================


%-----------------------------------------------------------------
%GET SATELLITE POSITIONS FROM SP3:
%-----------------------------------------------------------------
% READ SP3 File (Reads in all the data at once)
%DataFileName = 'NG05JA08.sp3';

%DataFileName = 'NG05JA05.SP3';

%[NumberSVs, ValidData, NumberEpochs, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data] = readsp3(DataFileName);

%===================================================================

%OR

%-------------------------------------------------------------------
%CALCULATE SATELLITE POSITIONS FROM Rinex Nav file
%-------------------------------------------------------------------

%Using GPSOrbitPropogator to propogate for all the satellites you can
%possibly calculate for in the .nav file. Have to select the specific PRN's
%to use at each epoch based on ValidDataRinexObs


SVWeek = GPStimeAMicroWeek;
SVTime_GPSSecs = GPStimeAMicro;


%determine which set of nav data to use, within 2 hours of the SVTime_GPSSecs






%PR's in Rinex file from Orion are at the UTC second not GPS time so add 13
%seconds 
%SVTime_GPSSecs = SVTime_GPSSecs+13;
% % % 

%initialise with zeros (makes it run faster in matlab)


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


for i = 1:NumSatEpochs  %epoch counter in seconds , number of points to calculate satellite positions for.
    
    for SV = 1:32
        
        PR = C1_PRNAMicro(SV,i);
        
        SVTime_GPSSecsOrb(i) = SVTime_GPSSecs(i) - PR/c; %pass the time of transmission to the GPSOrbitPropagator

        [SV_X_Data(i,SV) SV_Y_Data(i,SV) SV_Z_Data(i,SV) SV_T_Data(i,SV) ValidData_Satellite(i,SV)] = GPSOrbitPropagator(SVWeek(i), SVTime_GPSSecsOrb(i), SV, NavData);

        [SV_Xvel_Data(i,SV) SV_Yvel_Data(i,SV) SV_Zvel_Data(i,SV) SV_Tvel_Data(i,SV) SV_Xacc_Data(i,SV) SV_Yacc_Data(i,SV) SV_Zacc_Data(i,SV) SV_Tacc_Data(i,SV) ValidDataSatVels(i,SV)] = GPSOrbitPropagatorVelocities(SVWeek(i), SVTime_GPSSecsOrb(i), SV, NavData);
       
    end
end



%=======================================================================
%USER DEFINED INPUTS
%=======================================================================



%----------------------------------------------------------------
%Read in User preferences re: how the simulation will run
%Various flags for settings related to operation of simulation
%eg parameter to run solution from observed measurements, simulated or both






%READ IN Flight Planning Data
%This section reads in parameters for flight planning, such as
%Can read in longitude, latitude and height and time and this is converted to ECEF
%for processing.


% NumEpochs = 8380; %number of epochs to use %note if Sp3 used then each epoch is 15 minutes.
NumEpochs = 200; %number of epochs to use %note if Sp3 used then each epoch is 15 minutes.

%Standard duration is in seconds


%for i = 1: NumEpochs
%user x pos =
%user y pos =
%user z pos =


%--------------------------------------------
%RAIM and Least Squares-Specific Parameters
%--------------------------------------------

%MAximum allowable alarm rate = 0.002 per hour. %Constant alarm rate.
%Probability with any independante sample will be no greater than 1/15000.
PFalseAlarm = 1/15000; %From RTCA MOPS:
SigmaS = 20; %metres , Standard deviation of pseudorange measurements, this is receiver noise?
Alarm_Limit = 300;%300; %metres ..for particular phase of flight.
ElevationMask = 7.5 ; %7.5 degrees Satellite elevation mask to use for satellite accept or reject criteria





%=======================================================================
%ERROR MODELS FOR SIMULATED MEASUREMENTS
%These are parameters and calculations for generating simulated Pseudorange
%measurementss:
%All parameters relating to models for simulated measurements should be
%included in this section
%=======================================================================

%------------------------------------------------------------------------
%Atmospheric Models
%-----------------------------------------------------------------------

%IONOSPHERIC PARAMETERS for IONO Model
ALPHA = [0.1118e-07  -0.7451e-08  -0.5960e-07   .01192e-06];%          ION ALPHA
BETA = [0.1167e+06  -0.2294e+06  -0.1311e+06   .01049e+07]; %         ION BETA


%------------------------------------------------------------------------
%Receiver Clock model
%------------------------------------------------------------------------


% Generates clock bias and drift for one epoch only. Updates the previous clock bias with the new one. 

%start at 0
% UserPos(4) = 30/c;


%for i = 2:1000

%[ClockBias_t(i) ClockDrift_t(i)] = GARD_ClockModel(ClockBias_t(i-1)); %time in seconds, frequencies in hertz


%end


%integrate ClockDrift

%ClockDriftInt = cumtrapz(ClockDrift_t);


%calculate allan variance of the simulated clock



%work out allan variance of estimated receiver clock bias


%data is in 1 second sample length

% time = UserPos(4); %tau  = 1 second.
% 
% %average within each sample
% 
% difft = diff(time);
% 
% squared = difft.^2;
% 
% 
% %allan variance
% allanvar = 0.5*mean(squared);
% 
% %square root of allan variance
% allanvarsqrt = sqrt(allanvar);
% 
% %allan variance in ppb 
% 
% allanvarsqrtPPB = allanvarsqrt*1e9;





%------------------------------------------------------------------------
%Receiver Antenna model
%------------------------------------------------------------------------

%blah blah blah



%=======================================================================
%------------------------------------------------------------------------
%Write configuration file from Planning Module
%------------------------------------------------------------------------







%========================================================================
%========================================================================
% Variable Initialization
% General Variables to be initialised should be put in this section
%========================================================================
%========================================================================


PreviousSmoothed(1) = 0;
R_Ex = 0;



dTpos(1) = 30/c;


%===============================================================
%==========================================================================
%
%START CORE RUN MODULE
%This module does the number crunching
%==========================================================================
%===============================================================
% i = epoch counter


%initilaise some of the variables with zero (makes it run faster)

FinalSolution_Observed = zeros(NumEpochs,4);
FinalSolution_Simulated = zeros(NumEpochs,4);

SolutionVec_Observed = zeros(NumEpochs,4);
SolutionVec_Observed_Vel = zeros(NumEpochs,4);



EstimatedClockBias(1) = 0;


for i = 1:NumEpochs


       
    %  variables to be cleared from workspace at each iteration
    clear  PRMeasured_Observed SVPos PRMeasured_FDE_Observed...  %PR Observed
           PRMeasured_Simulated SVPos PRMeasured_FDE_Simulated... %PR Simulated
           PRRate_Observed PRRate SVVel PRRate_Simulated GeometricRange ICPMeasured_Simulated CPMeasured_Simulated ICPMeasured_Observed                 %PR rates observed & carrier phase.
       
          
    
           
    %-------------------------------------
    %ERROR CHECKING FLAGS

    RUNFLAG_Solution = 0;  %Solution can be made, (number of satellites is okay, run the code)
    RUNFLAG_FDE = 0;               %RAIM FDE can be done.
    %-------------------------------------


    %True User Position and clock error. (from previous estimate of clock error from clock model)
    Xpos(i) = -5046773.574;  %put      antenna position in here and can verify that the nav. solution is correct or not.
    Ypos(i) = 2568446.555;
    Zpos(i) = -2925288.451;
    %dTpos(i) = 0;
    
    
    %clock model for predicted PR measurements, calculates one epoch in advance
    [dTpos(i+1) dTposdot(i+1)] = GARD_ClockModel(dTpos(i)); %time in seconds, frequencies in hertz
    
    
    
    
    Xvel(i) = 0;
    Yvel(i) = 0;
    Zvel(i) = 0;
    dTvel(i) = 0;
    
    
%     Xpos(i) = -5045790.38107811;
%      Ypos(i) = 2570291.64376626;
%       Zpos(i) = -2922762.10591193;
                           


%     Xpos(i) = 0;  %put      antenna position in here and can verify that the nav. solution is correct or not.
%      Ypos(i) = 0;
%      Zpos(i) = 0;
    %Xpos(i) = 0;  %put      antenna position in here and can verify that the nav. solution is correct or not.
    %Ypos(i) = 0;
    %Zpos(i) = 0;


    %dTpos(i) = ClockBias_t(i)*Speedoflight;  %in metres
  


    %Update next Clock bias estimate with previous estimate from observed
    %pseudorange solution
    %dTpos(i) = FinalSolution_Observed(i,4);


    %this is the real, True user position and clock bias  based on flight planning inputs

    UserPos = ([Xpos(i),Ypos(i),Zpos(i),c*dTpos(i)]); %this will be a copy of the values set by the user in the flight planning section


    UserVel = ([Xvel(i),Yvel(i),Zvel(i),c*dTvel(i)]);

    %----------------------------------------------------
    %To Test FDE
    %Put error on a pseudorange measurement

    %--------------------------------------------------------


    %Vector of satellite PRNS to use in the solution, from the Rinex Obs
    %file


    %if using simulated measurements then   ValidData_SatelliteVector = ValidData_Satellite
    %else

    %else, sort which satellites to use at this epoch based on what PR
    %measurements are available.
    %Note:not using ValidData_Satellite because taking the satellites to
    %use from Rinex observation file

    
    %original (uncomment the following if not doing comparison of 2 sets of
    %observations
    
%     Valid = size(ValidDataRinexObs);
% 
%     for PRN = 1:Valid(1)
%         if ValidDataRinexObs(PRN,i) == 1
%             ValidData_SatelliteVector(i,PRN) = 1;
%         end
%     end
    
    
    
    
    
    %for comparing two sets of measurements from observations.
        %check for which satellties exist in the      receiver's file
    Valid = size(ValidDataRinexObs);

    for PRN = 1:Valid(1)
        if ValidDataRinexObs(PRN,i) == 1
            ValidData_SatelliteVector(i,PRN) = 1;
        end
    end
    
    
    %check for which satellties exist in the Room receivers file
%    Valid = size(ValidDataRinexObsRoom);
% 
%     for PRN = 1:Valid(1)
%         if ValidDataRinexObsRoom(PRN,i) == 1
%             ValidData_SatelliteVectorRoom(i,PRN) = 1;
%         end
%     end


    
    %Choose which satellites to use which are common to both      and room
    %receivers for the repeater test
    
%     
%     if ValidData_SatelliteVectorRoof(i,Counter) == 1 && ValidData_SatelliteVectorRoom(i,Counter) == 1 ValidData_Satellite(i,Counter) == 1
%     
    
    
    
    
    
    
    
    
    


    %======================================================================
    %Check data Read in and select satellites to use in solution.

    %is below elevation mask ??

    XApos = Xpos(i);
    YApos = Ypos(i);
    ZApos = Zpos(i);


    el = atan2(ZApos,sqrt(XApos^2 + YApos^2)); %lat, lon, hgt in ECEF
    az = atan2(YApos,XApos);
    hgt = sqrt(XApos^2 + YApos^2 +ZApos^2);

    %direction cosines (page 202 of Kathon this is the transformation
    %matrix for from ECEF to LTP. (East , North, Up)

    
     TMatrixECEFtoENU = T_ECEF2ENU(az,el);
%     sel = sin(el);
%     cel = cos(el);
%     saz = sin(az);
%     caz = cos(az);
%     Telev(1,1) = -sel*caz;
%     Telev(1,2) = -sel*saz;
%     Telev(1,3) = cel;
%     Telev(2,1) = -saz;
%     Telev(2,2) = caz;
%     Telev(2,3) = 0.0;
%     Telev(3,1) = cel*caz;
%     Telev(3,2) = cel*saz;
%     Telev(3,3) = sel;

    nn = 0;

    %transform xyz to spherical coords and NEU

    
        
    Valid = size(ValidData_SatelliteVector);
    
    %Valid = size(ValidData_SatelliteVector);

    for Counter = 1:Valid(2)

        %if the satellite data is in the Observation file AND there is data
        %from the navigation file then:
        if ValidData_SatelliteVector(i,Counter) == 1 & ValidData_Satellite(i,Counter) == 1  %the original

        %checks      observation file, room observation file, and     
        %navigation file (use a common nav file)
        %if ValidData_SatelliteVectorRoof(i,Counter) == 1 && ValidData_SatelliteVectorRoom(i,Counter) == 1 && ValidData_Satellite(i,Counter) == 1 %for repeater testing
            
            GPS_Vec(Counter,1) = SV_X_Data(i,Counter) - XApos;
            GPS_Vec(Counter,2)= SV_Y_Data(i,Counter) - YApos;
            GPS_Vec(Counter,3) = SV_Z_Data(i,Counter) - ZApos;
            
            
            
            
           %[topox(Counter); topoy(Counter); topoz(Counter)] =  
           
           topo = TMatrixECEFtoENU*[GPS_Vec(Counter,1); GPS_Vec(Counter,2); GPS_Vec(Counter,3)];
            
           
          topox(Counter) = topo(1);
          topoy(Counter) = topo(2);
          topoz(Counter) = topo(3);
            
            
%             topox(Counter) = Telev(1,1)*GPS_Vec(Counter,1) + Telev(1,2)*GPS_Vec(Counter,2) + Telev(1,3)*GPS_Vec(Counter,3);
%             topoy(Counter) = Telev(2,1)*GPS_Vec(Counter,1) + Telev(2,2)*GPS_Vec(Counter,2) + Telev(2,3)*GPS_Vec(Counter,3);
%             topoz(Counter) = Telev(3,1)*GPS_Vec(Counter,1) + Telev(3,2)*GPS_Vec(Counter,2) + Telev(3,3)*GPS_Vec(Counter,3);
            
%             topox(Counter) = Telev(1,1)*SV_X_Data(i,Counter) + Telev(1,2)*SV_Y_Data(i,Counter) + Telev(1,3)*SV_Z_Data(i,Counter);
%             topoy(Counter) = Telev(2,1)*SV_X_Data(i,Counter) + Telev(2,2)*SV_Y_Data(i,Counter) + Telev(2,3)*SV_Z_Data(i,Counter);
%             topoz(Counter) = Telev(3,1)*SV_X_Data(i,Counter) + Telev(3,2)*SV_Y_Data(i,Counter) + Telev(3,3)*SV_Z_Data(i,Counter);

            % Calculate the SV elevation

            rangevecxy = sqrt(topox(Counter)^2 + topoy(Counter)^2);
            elevation(i,Counter) = (180/pi)*atan2(topoz(Counter),rangevecxy) ;
            
            %Alternative calculation of elevation
            %rangeveccheck = sqrt(topox(Counter)^2 + topoy(Counter)^2 +topoz(Counter)^2);
%             elevationcheck(i,Counter) = (180/pi)*asin(topoz(Counter)/rangeveccheck) ;
%             
            %calculate azimuth - not sure if this is right or if it should
            %be atan2(topox(Counter),topoy(Counter))
            azimuth(i,Counter) = (180/pi)*atan2(topoy(Counter),topox(Counter)); 

            %satellites above elevation mask at user location are:

            if elevation(i,Counter) > ElevationMask   %7.5 degree elevation mask
                nn = nn+1;
                SV_AboveElevationMask(nn) = Counter ;

            end

        end

    end

        
    
    %check there's enough satellites for RAIM and Position Solution
    NumSatsForSolution(i) = nn;  %number of satellites in view to use in the solution
    if NumSatsForSolution(i) < 4; %Can't do solution
        ALARM_NotEnoughSatsSolution(i) = 100;
        %error('Not enough satellites for navigation solution')
        RUNFLAG_Solution = 1;  %If RUNFLAG = 1 then don't attempt a solution for this epoch, cannot be done.
    else if NumSatsForSolution(i) < 6;
            ALARM_NotEnoughSatsRAIMFDE(i) = 100;
            RUNFLAG_FDE = 1;
            %error('Not enough satellites for RAIM FDE solution')
        end
    end



    %=====================================================================

     if RUNFLAG_Solution == 0 % If number of satellites is OK for a least squares solution then continue



        N = NumSatsForSolution(i);
        
        %N = 6; %Use for now for the purpose of doing underdetermined velocity solutions 
        

        for uu = 1:N
            Sats(uu) = SV_AboveElevationMask(uu);
        end
        
        
       % 6    15    16    18    21    22    26    29
        
        
%         Sats = [21 16 6 22];
%         N = 4;
%         
        
              %--------------------------------------------------------------------------
%         %Generate Matrix of satellite positions and Predicted and Observed Pseudorange Vectors.
%         %----------------------------------------------------------------
%         -----



          
%         
% 
%              ClockBias_t(1) = 84982.380068232/c;
%             ClockBias_t(2) =   84982.3846606345/c;
%           ClockBias_t(3) =  84982.4053022903/c;
%          ClockBias_t(4) =   84982.3926327852/c;
%           ClockBias_t(5) =   84982.4258549901/c;



           
    if i == 1 %have to initialise these:
        PRMeasured_SimulatedPrevious = zeros(NumEpochs,N) ;
        GeometricRange_SimulatedPrevious = zeros(NumEpochs,N);
        CPMeasured_SimulatedPrevious = zeros(NumEpochs,N);
        ICPMeasured_SimulatedPrevious = zeros(NumEpochs,N);
    end      
          
    
    
         
         aa = PRMeasured_SimulatedPrevious(i,:);
         qq = ICPMeasured_SimulatedPrevious(i,:);      %IntegratedCarrierPhase
         bb = CPMeasured_SimulatedPrevious(i,:);
         cc = GeometricRange_SimulatedPrevious(i,:);
                
        
        
         [SVPos,PRMeasured_Observed,PRMeasured_Simulated,PRRate_Simulated,ICPMeasured_Observed,ICPMeasured_Simulated,CPMeasured_Simulated,NAmb, GeometricRange, SVVel, PRRate_Observed] = GARDSim_GenerateMeasurements(i,GPStimeAMicro(i),UserPos,N, Sats,C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data, SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data, ALPHA, BETA,aa,qq,bb,cc,EstimatedClockBias(i));
        
         
         
        

        %=================================================================
        %LEAST SQUARES and RAIM Solution from REAL Measurements
        %=================================================================

        %Do Least squares and RAIM FDE solution For Observed Pseudoranges
        %from Rinex Files

        [SolutionVec_Observed(i,1:4) VarSolutionVec_Observed NumIterations_Observed_Vel(i) ResVec_Observed M_Observed LSQ_Fail_Observed(i) limit_Observed(i) DOP_Observed(i,1:4)] = GARD_LSQ(UserPos,N,PRMeasured_Observed,SVPos);

      
        %Calculate GPS velocity from range rates which were calculated from
        %the carrier phase
        
        %note that the velocity calculation here is done BEFORE any RAIM checking is done.  You could put this vel. calculation after the FD or FDE is done on the pseudorange measurements.
                
        [SolutionVec_Observed_Vel(i,1:4) VarSolutionVec_Observed_Vel NumIterations_Observed_Vel(i) ResVec_Observed_Vel M_Observed_Vel LSQ_Fail_Observed_Vel(i) limit_Observed_Vel(i)] = GARD_LSQVel(SolutionVec_Observed(i,1:4),UserVel,N,PRRate_Observed,SVPos,SVVel);

        
        %calculate GPS velocity from range rates but for underdetermined
        %case (receiver clock drift is estimated by differentiating
        %receiver clock bias
%         SolutionVec_ObservedUDS(i,1:4) = SolutionVec_Observed(i,1:4);
%         
%         SolutionVec_ObservedUDS(i,4) = SolutionVec_Observed(i,4) - SolutionVec_Observed(i-1,4);
%         
%       
%         
%         
%         [SolutionVec_Observed_VelUDS(i,1:4) VarSolutionVec_Observed_VelUDS NumIterations_Observed_VelUDS(i) ResVec_Observed_VelUDS M_Observed_VelUDS LSQ_Fail_Observed_VelUDS(i) limit_VelUDS(i)] = GARD_LSQVelUDS(SolutionVec_ObservedUDS(i,1:4),UserVel,N,PRRate_Observed,SVPos,SVVel);
% 
%         
%         
%           %make it calculate velocity with only 2 unknowns ie 2
%           %measurements required
%           %estimate z by differentiation z and dt
%         SolutionVec_ObservedUDSTwoMeas(i,1:4) = SolutionVec_ObservedUDS(i,1:4);
%         
%         %SolutionVec_ObservedUDSTwoMeas(i,3) = SolutionVec_Observed(i,3) - SolutionVec_Observed(i-1,3);
%         Zterm = SolutionVec_Observed(i,3) - SolutionVec_Observed(i-1,3);
%         
%         [SolutionVec_Observed_VelUDSTwoMeas(i,1:4) VarSolutionVec_Observed_VelUDSTwoMeas NumIterations_Observed_VelUDSTwoMeas(i) ResVec_Observed_VelUDSTwoMeas M_Observed_VelUDSTwoMeas LSQ_Fail_Observed_VelUDSTwoMeas(i) limit_VelUDSTwoMeas(i)] = GARD_LSQVelUDSTwoMeas(SolutionVec_ObservedUDSTwoMeas(i,1:4),UserVel,N,PRRate_Observed,SVPos,SVVel,Zterm);

        
        if LSQ_Fail_Observed(i) == 0  %If solution has converged
            %careful with this because if theres
            %no LSQ then it wont increment
            %RAIM_ALARM_Observed and others, so an else
            %statement has been put below

            FinalSolution_Observed(i,1) = SolutionVec_Observed(i,1);
            FinalSolution_Observed(i,2) = SolutionVec_Observed(i,2);
            FinalSolution_Observed(i,3) = SolutionVec_Observed(i,3);
            FinalSolution_Observed(i,4) = SolutionVec_Observed(i,4);

            %Do RAIM FD solution using least squares residuals method

            %RAIM_ALARM == 100 = true
            %RAIM_ALARM == Anything else = False

            %need to implement parity method, and parity isolation method. cf
            %parkinson blue book.

            
            
            
            [BadGeometry_Observed(i) RAIM_ALARM_Observed(i) SLOPE_Max_Observed(i) r_Observed(i) r_Observed_Threshold(i) ARP_Observed(i)] = GARD_RAIM(N,PFalseAlarm,SigmaS,Alarm_Limit,ResVec_Observed,M_Observed);

            
            
            %Using the parity method
            [BadGeometry_ObservedPar(i) RAIM_ALARM_ObservedPar(i) SLOPE_Max_ObservedPar(i) r_ObservedPar(i) r_Observed_ThresholdPar(i) ARP_ObservedPar(i) FaultySatFDIPar(i)] = GARDSim_RAIMParity(N,PFalseAlarm,SigmaS,Alarm_Limit,ResVec_Observed,M_Observed,PRMeasured_Observed);

                
            %The faulty sat is identified as:
            
            if RAIM_ALARM_ObservedPar(i) == 100
            
            FaultySat(i) = Sats(FaultySatFDIPar(i));  %PRN number of the faulty satellite
            
        else
            FaultySat(i) = 0;
        end
        
            
            
            
            
        else
            BadGeometry_Observed(i) = 0;
            RAIM_ALARM_Observed(i)= 0;
            SLOPE_Max_Observed(i)= 0;
            r_Observed(i)= 0;

            %don't need this for these if initialising with zeros before
            %the loop starts
            FinalSolution_Observed(i,1) = 0;
            FinalSolution_Observed(i,2) = 0;
            FinalSolution_Observed(i,3) = 0;
            FinalSolution_Observed(i,4) = 0;
            
%             SolutionVec_Observed(i,1) = 0;
%             
%             
%             
%             SolutionVec_Observed_Vel(i,1) = 0;


        end
        
        
        
        %output results from velocity calculated from observed measurements.
        
        
        if LSQ_Fail_Observed_Vel(i) == 0  %If solution has converged
            %careful with this because if theres
            %no LSQ then it wont increment
            %RAIM_ALARM_Observed and others, so an else
            %statement has been put below

            FinalSolution_ObservedVel(i,1) = SolutionVec_Observed_Vel(i,1);
            FinalSolution_ObservedVel(i,2) = SolutionVec_Observed_Vel(i,2);
            FinalSolution_ObservedVel(i,3) = SolutionVec_Observed_Vel(i,3);
            FinalSolution_ObservedVel(i,4) = SolutionVec_Observed_Vel(i,4);
        
        
        
        end
        
        
        
              
        
          
        
        
        
        
        

        %==================================================================

        %Do RAIM FDE solution for Observed PRs if required

        %FDE taken from Parkinson Volume 2

        %if failure...go through all the available subsets...
        % if 6 satellites, can only detect one error at a time.

        % if error occurs on more than one subset then have to assume that more
        % than one satellite failure has occurred.


        %Do FDE solution
        %if failure:

        if RUNFLAG_FDE == 0  %Run only if there's enough satellites for FDE (6 or more).

            if RAIM_ALARM_Observed(i) == 100 %If there's a RAIM alarm attempt FDE solution

                %Do FDE - ONLY DETECTS ONE SATELLITE FAILURE AT A TIME AT THE MOMENT

                %Form N-1 Subsets
           
              SimulatedOrObserved_Flag = 0; %Use observed PR measurements
              
             
                
                
             [FinalSolution_FDE_Obs RAIM_ALARM_FDE_Obs GoodSatsSats_FDE_Obs Sats_FDE_Obs DOP_Observed(i,1:4)] = GARDSim_FDE(SimulatedOrObserved_Flag,i,GPStimeAMicro(i),UserPos,N, Sats,C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data,SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data,PFalseAlarm,SigmaS,Alarm_Limit,ALPHA,BETA,aa,qq,bb,cc,EstimatedClockBias(i));
            
             
         
         
             
            FinalSolution_Observed(i,1) = FinalSolution_FDE_Obs(1);
            FinalSolution_Observed(i,2) = FinalSolution_FDE_Obs(2);
            FinalSolution_Observed(i,3) = FinalSolution_FDE_Obs(3);
            FinalSolution_Observed(i,4) = FinalSolution_FDE_Obs(4);
               
      


            else   %if there's no RAIM alarm
%             GoodSatsSats_FDE_Obs(i,:) = 0;
%             RAIM_ALARM_FDE_Obs(i,:) = 0;
            
            end  %If there's a RAIM alarm attempt FDE solution


        end %RUNFLAG_FDE == 0   %End FDE Section for observed measurements




      
  
        %=================================================================
        %LEAST SQUARES and RAIM Solution from SIMULATED Measurements
        %=================================================================
        
       
        


        [SolutionVec_Simulated(i,1:4) VarSolutionVec_Simulated NumIterations_Simulated(i) ResVec_Simulated M_Simulated LSQ_Fail_Simulated(i) limit_Simulated(i) DOP_Simulated(i,1:4)] = GARD_LSQ(UserPos,N,PRMeasured_Simulated,SVPos);

        
                   
        [SolutionVec_Simulated_Vel(i,1:4) VarSolutionVec_Simulated_Vel NumIterations_Simulated_Vel(i) ResVec_Simulated_Vel M_Simulated_Vel LSQ_Fail_Simulated_Vel(i) limit_Simulated_Vel(i)] = GARD_LSQVel(SolutionVec_Simulated(i,1:4),UserVel,N,PRRate_Simulated,SVPos,SVVel);

        
        
        
       
        if LSQ_Fail_Simulated(i) == 0   %If solution has converged

            FinalSolution_Simulated(i,1) = SolutionVec_Simulated(i,1);
            FinalSolution_Simulated(i,2) = SolutionVec_Simulated(i,2);
            FinalSolution_Simulated(i,3) = SolutionVec_Simulated(i,3);
            FinalSolution_Simulated(i,4) = SolutionVec_Simulated(i,4);
            
            EstimatedClockBias(i+1) = FinalSolution_Simulated(i,4); 

            %Do RAIM FD solution using least squares residuals method

            %RAIM_ALARM == 100 = true
            %RAIM_ALARM == Anything else = False

            %need to implement parity method, and parity isolation method. cf
            %parkinson blue book.

            [BadGeometry_Simulated(i) RAIM_ALARM_Simulated(i) SLOPE_Max_Simulated(i) r_Simulated(i) r_Simulated_Threshold(i) ARP_Simulated(i)] = GARD_RAIM(N,PFalseAlarm,SigmaS,Alarm_Limit,ResVec_Simulated,M_Simulated);

            
        else

            BadGeometry_Simulated(i) = 0;
            RAIM_ALARM_Simulated(i)= 0;
            SLOPE_Max_Simulated(i)= 0;
            r_Simulated(i)= 0;

            FinalSolution_Simulated(i,1) = 0;
            FinalSolution_Simulated(i,2) = 0;
            FinalSolution_Simulated(i,3) = 0;
            FinalSolution_Simulated(i,4) = 0;
            
            


        end
        
        
        
         %output results from velocity calculated from simulated measurements.
        
        
        if LSQ_Fail_Simulated_Vel(i) == 0  %If solution has converged
            %careful with this because if theres
            %no LSQ then it wont increment
            %RAIM_ALARM_Observed and others, so an else
            %statement has been put below

            FinalSolution_SimulatedVel(i,1) = SolutionVec_Simulated_Vel(i,1);
            FinalSolution_SimulatedVel(i,2) = SolutionVec_Simulated_Vel(i,2);
            FinalSolution_SimulatedVel(i,3) = SolutionVec_Simulated_Vel(i,3);
            FinalSolution_SimulatedVel(i,4) = SolutionVec_Simulated_Vel(i,4);
        
        
        
        end
        
        
        
        
        
        
        
        
        

        %==================================================================

     %Do FDE solution
        %if failure:

        if RUNFLAG_FDE == 0  %Run only if there's enough satellites for FDE (6 or more).

            if RAIM_ALARM_Simulated(i) == 100 %If there's a RAIM alarm attempt FDE solution

                %Do FDE - ONLY DETECTS ONE SATELLITE FAILURE AT A TIME AT THE MOMENT

                %Form N-1 Subsets

             SimulatedOrObserved_Flag = 1; %Use simulated PR measurements
                
             [FinalSolution_FDE_Sim RAIM_ALARM_FDE_Sim GoodSatsSats_FDE_Sim Sats_FDE_Sim DOP_Simulated(i,1:4)] = GARDSim_FDE(SimulatedOrObserved_Flag,i,GPStimeAMicro(i),UserPos,N, Sats,C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data,SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data,PFalseAlarm,SigmaS,Alarm_Limit,ALPHA,BETA,aa,qq,bb,cc,EstimatedClockBias(i));
               
             
           
                             
            FinalSolution_Simulated(i,1) = FinalSolution_FDE_Sim(1);
            FinalSolution_Simulated(i,2) = FinalSolution_FDE_Sim(2);
            FinalSolution_Simulated(i,3) = FinalSolution_FDE_Sim(3);
            FinalSolution_Simulated(i,4) = FinalSolution_FDE_Sim(4);
               
             
            
             
            else   %if there's no RAIM alarm
%             GoodSatsSats_FDE_Sim(i,:) = 0;   
%             RAIM_ALARM_FDE_Sim(i,:) = 0;
            
            end  %If there's a RAIM alarm attempt FDE solution


        end %RUNFLAG_FDE == 0   %End FDE Section for observed measurements



        %Update next Clock bias estimate with previous estimate from observed
        %pseudorange solution
        %dTpos(i+1) = FinalSolution_Observed(i,4);


        %Take copy of PR's
        jjjjj = size(PRMeasured_Simulated);  %can use simulated or observed here..just to get the array size


        %Copy any local 'single epoch' variables into Global Variables for later analysis/plotting
        for pp = 1:jjjjj(2)
            
            %Satellites used at each iteration
                       
            SatsGlobal(i,pp) = Sats(pp);
            
            
            %Pseudoranges
             PRMeasured_ObservedGlobal(i,pp) = PRMeasured_Observed(pp);
             PRMeasured_SimulatedGlobal(i,pp) = PRMeasured_Simulated(pp);
            
            %Pseudorange Rates
            PRRate_ObservedGlobal(i,pp) = PRRate_Observed(pp);
            %Simulated Pseduorange rates to go in here
            
            PRRate_SimulatedGlobal(i,pp) = PRRate_Simulated(pp);
            
            
            
            
            
            %Carrier phases
            ICPMeasured_SimulatedGlobal(i,pp) = ICPMeasured_Simulated(pp); 
            ICPMeasured_ObservedGlobal(i,pp) = ICPMeasured_Observed(pp); 
                         
            
            %needed for simulated PR rate and carrier phase calculations
            PRMeasured_SimulatedPrevious(i+1,pp) = PRMeasured_Simulated(pp);
            GeometricRange_SimulatedPrevious(i+1,pp) = GeometricRange(pp);
            CPMeasured_SimulatedPrevious(i+1,pp) = CPMeasured_Simulated(pp);    
            ICPMeasured_SimulatedPrevious(i+1,pp) = ICPMeasured_Simulated(pp);
            
            
          
            
            
            %Smoothed psuedoranges
            %PreviousSmoothedGlobal(i,pp) = PreviousSmoothed(pp); %for smoothed PR solution
             %global copy of satellites used
             
            
            
            
            
        end

    end %for i = 1:NumEpochs

end %if RUNFLAG_Solution ==1






%===============================================================
%===============================================================
%OUTPUT MODULE
% This module is the module where outputs are displayed to the user
%===============================================================
%===============================================================
%XYZ Position error from true solution
% PosErr(i) = sqrt( (SolutionVec(1)-UserPos(1))^2 + (SolutionVec(2)-UserPos(2))^2 + (SolutionVec(3)-UserPos(3))^2);
%
% if RAIM_ALARM(i) == 100
% plot(r(i),PosErr(i),'or');
% else
% plot(r(i),PosErr(i),'o');
% end
%
%
%
% xlabel('range residual/sqrt(n-4) (m)');
% ylabel('radial position error (m)');

%Check Accuracy

%True position (TruthdT is unknown). 

%This is the ashtech uZ surveyed antenna position on S block roof.
TruthX = -5046773.57;
TruthY = 2568446.555;
TruthZ = -2925288.451;
TruthdT = 0;

TruthXdot = 0;
TruthYdot = 0;
TruthZdot = 0;
TruthdTdot = 0;




%compare positions


Xerr_Observed = TruthX - FinalSolution_Observed(:,1);
Yerr_Observed = TruthY - FinalSolution_Observed(:,2);
Zerr_Observed = TruthZ - FinalSolution_Observed(:,3);
dTerr_Observed = TruthdT - FinalSolution_Observed(:,4);


Xerr_Simulated = TruthX - FinalSolution_Simulated(:,1);
Yerr_Simulated = TruthY - FinalSolution_Simulated(:,2);
Zerr_Simulated = TruthZ - FinalSolution_Simulated(:,3);
dTerr_Simulated = TruthdT - FinalSolution_Simulated(:,4);


subplot(2,1,1); plot(Xerr_Observed); title 'X Error Obs (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(Xerr_Simulated); title 'X Error Sim (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;

%Difference between observed and simulated. 

plot(Xerr_Observed - Xerr_Simulated);  title 'X Error Obs - X Error Sim (m)'; xlabel('Time (Secs)'); 

pause;

subplot(2,1,1); plot(Yerr_Observed); title 'Y Error Obs (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(Yerr_Simulated); title 'Y Error Sim (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;

plot(Yerr_Observed - Yerr_Simulated);  title 'Y Error Obs - Y Error Sim (m)'; xlabel('Time (Secs)'); 
 pause;

subplot(2,1,1); plot(Zerr_Observed); title 'Z Error Obs (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(Zerr_Simulated); title 'Z Error Sim (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;

plot(Zerr_Observed - Zerr_Simulated);  title 'Z Error Obs - Z Error Sim (m)'; xlabel('Time (Secs)'); 
 pause;
subplot(2,1,1); plot(dTerr_Observed); title 'dT Error Obs (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(dTerr_Simulated); title 'dT Error Sim (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;

plot(dTerr_Observed - dTerr_Simulated);  title 'dT Error Obs - dT Error Sim (m)'; xlabel('Time (Secs)'); 
pause;





%compare velocities



Xerr_Observeddot = TruthXdot - FinalSolution_ObservedVel(:,1);
Yerr_Observeddot = TruthYdot - FinalSolution_ObservedVel(:,2);
Zerr_Observeddot = TruthZdot - FinalSolution_ObservedVel(:,3);
dTerr_Observeddot = TruthdTdot - FinalSolution_ObservedVel(:,4);


Xerr_Simulateddot = TruthXdot - FinalSolution_SimulatedVel(2:NumEpochs,1);
Yerr_Simulateddot = TruthYdot - FinalSolution_SimulatedVel(2:NumEpochs,2);
Zerr_Simulateddot = TruthZdot - FinalSolution_SimulatedVel(2:NumEpochs,3);
dTerr_Simulateddot = TruthdTdot - FinalSolution_SimulatedVel(2:NumEpochs,4);


%compare velocity errors by differentiating the position errors. 


Xerr_Observeddotdiff = diff(Xerr_Observed);
Yerr_Observeddotdiff = diff(Yerr_Observed);
Zerr_Observeddotdiff = diff(Zerr_Observed);
dTerr_Observeddotdiff = diff(dTerr_Observed);


%There is a 1 second delay due to the differentiation...hence i have to adjust these..
Xerr_Simulateddotdiff = diff(Xerr_Simulated);
Yerr_Simulateddotdiff = diff(Yerr_Simulated);
Zerr_Simulateddotdiff = diff(Zerr_Simulated);
dTerr_Simulateddotdiff = diff(dTerr_Simulated);



subplot(2,1,1); plot(Xerr_Observeddot); title 'X dot Error Obs (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(Xerr_Simulateddot); title 'X dot Error Sim (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;

%Difference between observed and simulated. 

plot(Xerr_Observeddot(2:NumEpochs) - Xerr_Simulateddot);  title 'X dot Error Obs - X dot Error Sim (m)'; xlabel('Time (Secs)'); 

pause;

subplot(2,1,1); plot(Yerr_Observeddot); title 'Y dot Error Obs (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(Yerr_Simulateddot); title 'Y dot Error Sim (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;

plot(Yerr_Observeddot(2:NumEpochs) - Yerr_Simulateddot);  title 'Y dot Error Obs - Y dot Error Sim (m)'; xlabel('Time (Secs)'); 
 pause;

subplot(2,1,1); plot(Zerr_Observeddot); title 'Z dot Error Obs (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(Zerr_Simulateddot); title 'Z dot Error Sim (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;

plot(Zerr_Observeddot(2:NumEpochs) - Zerr_Simulateddot);  title 'Z Error Obs - Z Error Sim (m)'; xlabel('Time (Secs)'); 
 pause;
subplot(2,1,1); plot(dTerr_Observeddot); title 'dT dot Error Obs (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(dTerr_Simulateddot); title 'dT dot Error Sim (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;

plot(dTerr_Observeddot(2:NumEpochs) - dTerr_Simulateddot);  title 'dT dot Error Obs - dT dot Error Sim (m)'; xlabel('Time (Secs)'); 
pause;


%compare differentiated positions with velocities


hh = size(Xerr_Observeddot);

plot(Xerr_Observeddot(2:hh(1)) - Xerr_Observeddotdiff);  title 'X dot Error Obs - X dot Error Obs Diff (m)'; xlabel('Time (Secs)'); 

pause;


plot(Yerr_Observeddot(2:hh(1)) - Yerr_Observeddotdiff);  title 'Y dot Error Obs - Y dot Error Obs Diff (m)'; xlabel('Time (Secs)'); 

pause;


plot(Zerr_Observeddot(2:hh(1)) - Zerr_Observeddotdiff);  title 'Z dot Error Obs - Z dot Error Obs Diff (m)'; xlabel('Time (Secs)'); 

pause;


plot(dTerr_Observeddot(2:hh(1)) - dTerr_Observeddotdiff);  title 'dT dot Error Obs - dT dot Error Obs Diff (m)'; xlabel('Time (Secs)'); 

pause;


%Compare simulated clock model 
%dTpos dTposdot


subplot(2,1,1); plot(dTpos*c); title 'dT Simulated Clock Model (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(dTerr_Observed); title 'dT Error Observed (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;



subplot(2,1,1); plot(dTposdot*c); title 'dT dot Simulated Clock Model (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(dTerr_Observeddot); title 'dT dot Error Observed (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);

pause;
clf;







%compare pseudoranges


%difference between measured and predicted PR's
PRresult = PRMeasured_SimulatedGlobal-PRMeasured_ObservedGlobal;



%difference between measured and predicted PR's without receiver clock biases


q = size(PRMeasured_SimulatedGlobal);

for q = 1:q(2)
    

PRresultnoclock(:,q) = (PRMeasured_SimulatedGlobal(:,q) - FinalSolution_Simulated(:,4)) - (PRMeasured_ObservedGlobal(:,q)  - FinalSolution_Observed(:,4));


end




for g = 1:N


plot(PRresultnoclock(:,g));  title (['PR Sim - PR Obs No Rx Clock Bias (m) SV ', num2str(Sats(g))]); xlabel('Time (Secs)'); 

pause;


end









%compare carrier phase


%can't really compare CP because they both have different starting points, have to compare the difference between them. 
%difference between measured and predicted PR's
CPresult = diff(ICPMeasured_SimulatedGlobal)-diff(ICPMeasured_ObservedGlobal);





%quality analysis of  ICP in METRES

ICPSim = (ICPMeasured_SimulatedGlobal(10:200,1))*L1_Wavelength;

ICPObs = (ICPMeasured_ObservedGlobal(10:200))*L1_Wavelength;




[CurveFitSim, BestFitSim] = LeastSquaresBestFit(ICPSim', 5, 1);
[CurveFitObs, BestFitSim] = LeastSquaresBestFit(ICPObs, 5, 1);



pause;


subplot(2,1,1); plot(CurveFitSim - ICPSim'); title 'ICP Sim Residuals (m)'; %axis([1 3999+10 axesvaluemin axesvaluemax ]);
subplot(2,1,2); plot(CurveFitObs - ICPObs); title 'ICP Obs Residuals (m)'; xlabel('Time (Secs)'); %axis([1 3999+10 axesvaluemin axesvaluemax ]);


pause;

