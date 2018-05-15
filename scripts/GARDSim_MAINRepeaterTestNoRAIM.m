
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
 
%feature accel on

%always make NumSatEpochs greater than NumEpochs
NumSatEpochs = 4050 ;


%-----------------------------------------------------------------
%READ IN DATA FROM FILES
%-----------------------------------------------------------------

%READ RINEX
%Data from Ashtech receiver:

%Observation File
%Read roof ashtech micro z 24 hour

%Read Dual frequency data
%generic filename and path:

%8Sept
%for repeater Test
%filepathObs = 'data\RepeaterTest8Sept\RoofSess1.05O';
%filepathObs = 'data\RepeaterTest8Sept\RoofSess2.05O';



%15 Sept
%filepathObs = 'data\RepeaterTest15Sept\RoofSess1.05O';
%filepathObs = 'data\RepeaterTest15Sept\RoofSess2.05O';
filepathObs = 'data\RepeaterTest15Sept\RoofSess3.05O';




%filepathObs = 'data\RepeaterTest15Sept\RoofSess3ForCompareSatPositionstoCheckGPSPropagator.05O';

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
clear GPStimeAMicroTempRoof ValidDataRinexObsTempRoof L1_PRNAMicroTempRoof L2_PRNAMicroTempRoof C1_PRNAMicroTempRoof P1_PRNAMicroTempRoof P2_PRNAMicroTempRoof;

%File to read in  
FilenameTempObs = IsFilesObs(filecounterObs).name;

FilenameObs = fullfile(pathstrObs,FilenameTempObs)

[GPStimeAMicroWeekTempRoof, GPStimeAMicroTempRoof, ValidDataRinexObsTempRoof, L1_PRNAMicroTempRoof, L2_PRNAMicroTempRoof, C1_PRNAMicroTempRoof, P1_PRNAMicroTempRoof, P2_PRNAMicroTempRoof] = ReadRinexSblockRoof(FilenameObs);




%initialise 

if (filecounterObs == 1)
    
%add the temp data to the data to be processed
%initialise with zeros
GPStimeAMicroWeekRoof = zeros(size(GPStimeAMicroWeekTempRoof));  %GPS week
GPStimeAMicroRoof = zeros(size(GPStimeAMicroTempRoof));          % GPS seconds into week
ValidDataRinexObsRoof = zeros(size(ValidDataRinexObsTempRoof));
L1_PRNAMicroRoof  = zeros(size(L1_PRNAMicroTempRoof));
L2_PRNAMicroRoof  = zeros(size(L2_PRNAMicroTempRoof));
C1_PRNAMicroRoof = zeros(size(C1_PRNAMicroTempRoof));
P1_PRNAMicroRoof = zeros(size(P1_PRNAMicroTempRoof));
P2_PRNAMicroRoof = zeros(size(P2_PRNAMicroTempRoof));    
    
GPStimeAMicroWeekRoof = GPStimeAMicroWeekTempRoof;     
GPStimeAMicroRoof = GPStimeAMicroTempRoof;
ValidDataRinexObsRoof = ValidDataRinexObsTempRoof;
L1_PRNAMicroRoof = L1_PRNAMicroTempRoof;
L2_PRNAMicroRoof = L2_PRNAMicroTempRoof;
C1_PRNAMicroRoof = C1_PRNAMicroTempRoof;
P1_PRNAMicroRoof = P1_PRNAMicroTempRoof;
P2_PRNAMicroRoof = P2_PRNAMicroTempRoof;
end

if (filecounterObs ~= 1)

    %cocatenate the previous set of data with the next set
    
GPStimeAMicroWeekRoof = cat(2, GPStimeAMicroWeekRoof, GPStimeAMicroWeekTempRoof);    
GPStimeAMicroRoof = cat(2, GPStimeAMicroRoof, GPStimeAMicroTempRoof);
ValidDataRinexObsRoof = cat(2, ValidDataRinexObsRoof, ValidDataRinexObsTempRoof);
L1_PRNAMicroRoof  = cat(2, L1_PRNAMicroRoof, L1_PRNAMicroTempRoof);
L2_PRNAMicroRoof  = cat(2, L2_PRNAMicroRoof, L2_PRNAMicroTempRoof);
C1_PRNAMicroRoof  = cat(2, C1_PRNAMicroRoof, C1_PRNAMicroTempRoof);
P1_PRNAMicroRoof  = cat(2, P1_PRNAMicroRoof, P1_PRNAMicroTempRoof);
P2_PRNAMicroRoof  = cat(2, P2_PRNAMicroRoof, P2_PRNAMicroTempRoof);

end %if (filecounterObs ~= 1)


end; %end for filecounterObs = 1:NumFilesObs(1)


end; %end else



%----------------------------------------------------------------------


%READ RINEX Navigation File
%Data from Ashtech receiver:


%generic filename and path:


%for repeater Test

%8 Sept
%filepathNav = 'data\RepeaterTest8Sept\RoofSess1.05N';
%filepathNav = 'data\RepeaterTest8Sept\RoofSess2.05N';

%15 Sept
%filepathNav = 'data\RepeaterTest15Sept\RoofSess1.05N';
%filepathNav = 'data\RepeaterTest15Sept\RoofSess2.05N';
filepathNav = 'data\RepeaterTest15Sept\RoofSess3.05N';




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




%READ RINEX for single frequency measurements
%Data for 835 Assignment

%Observation File
%Read Orion data

%generic filename and path:

%filepathObs = 'data\835Assignment\Orion835.O';


 %for repeater Test
 %8 Sept
%filepathObs = 'data\RepeaterTest8Sept\RoomSess1.05O';
%filepathObs = 'data\RepeaterTest8Sept\RoomSess2.05O';



%15 Sept
%filepathObs = 'data\RepeaterTest15Sept\RoomSess1.05O';
%filepathObs = 'data\RepeaterTest15Sept\RoomSess2.05O';
filepathObs = 'data\RepeaterTest15Sept\RoomSess3.05O';




%filepathObs = 'data\RepeaterTest15Sept\RoomSess3ForCompareSatPositionstoCheckGPSPropagator.05O';

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
clear GPStimeAMicroTempRoom ValidDataRinexObsTempRoom L1_PRNAMicroTempRoom L2_PRNAMicroTempRoom C1_PRNAMicroTempRoom P1_PRNAMicroTempRoom P2_PRNAMicroTempRoom;

%File to read in  
FilenameTempObs = IsFilesObs(filecounterObs).name;

FilenameObs = fullfile(pathstrObs,FilenameTempObs)

%[GPStimeAMicroWeekTemp GPStimeAMicroTemp, ValidDataRinexObsTemp, L1_PRNAMicroTemp, L2_PRNAMicroTemp, C1_PRNAMicroTemp, P1_PRNAMicroTemp, P2_PRNAMicroTemp] = ReadRinexSblockRoof(FilenameObs);
[GPStimeAMicroWeekTempRoom,GPStimeAMicroTempRoom,ValidDataRinexObsTempRoom, L1_PRNAMicroTempRoom, C1_PRNAMicroTempRoom] = ReadRinexObsSuperStar(FilenameObs);



%initialise 

if (filecounterObs == 1)
    
%add the temp data to the data to be processed
%initialise with zeros
GPStimeAMicroWeekRoom = zeros(size(GPStimeAMicroWeekTempRoom));  %GPS week
GPStimeAMicroRoom = zeros(size(GPStimeAMicroTempRoom));          % GPS seconds into week
ValidDataRinexObsRoom = zeros(size(ValidDataRinexObsTempRoom));
L1_PRNAMicroRoom  = zeros(size(L1_PRNAMicroTempRoom));

C1_PRNAMicroRoom = zeros(size(C1_PRNAMicroTempRoom));
 
    
GPStimeAMicroWeekRoom = GPStimeAMicroWeekTempRoom;     
GPStimeAMicroRoom = GPStimeAMicroTempRoom;
ValidDataRinexObsRoom = ValidDataRinexObsTempRoom;
L1_PRNAMicroRoom = L1_PRNAMicroTempRoom;

C1_PRNAMicroRoom = C1_PRNAMicroTempRoom;
end

if (filecounterObs ~= 1)

    %cocatenate the previous set of data with the next set
    
GPStimeAMicroWeekRoom = cat(2, GPStimeAMicroWeekRoom, GPStimeAMicroWeekTempRoom);    
GPStimeAMicroRoom = cat(2, GPStimeAMicroRoom, GPStimeAMicroTempRoom);
ValidDataRinexObsRoom = cat(2, ValidDataRinexObsRoom, ValidDataRinexObsTempRoom);
L1_PRNAMicroRoom  = cat(2, L1_PRNAMicroRoom, L1_PRNAMicroTempRoom);
C1_PRNAMicroRoom  = cat(2, C1_PRNAMicroRoom, C1_PRNAMicroTempRoom);

end %if (filecounterObs ~= 1)


end; %end for filecounterObs = 1:NumFilesObs(1)


end; %end else



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
% %for repeater test
% filepathNav = 'data\RepeaterTest8Sept\RoofSess1.05N';
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


SVWeek = GPStimeAMicroWeekRoof;
SVTime_GPSSecs = GPStimeAMicroRoof;


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



%Use Roof or Repeater measurements to calculate position

RoofFlag = 0; 

%Switch between using observed measurements from roof receiver or room
%receiver

if RoofFlag == 1  %if is 1 then use roof measurements. if is 0 use room ones.

    
GPStimeAMicroWeek = GPStimeAMicroWeekRoof;
GPStimeAMicro = GPStimeAMicroRoof;
ValidDataRinexObs = ValidDataRinexObsRoof;
L1_PRNAMicro = L1_PRNAMicroRoof;
L2_PRNAMicro = L2_PRNAMicroRoof;
C1_PRNAMicro = C1_PRNAMicroRoof;
P1_PRNAMicro = P1_PRNAMicroRoof;
P2_PRNAMicro = P2_PRNAMicroRoof;    
    
    
end



if RoofFlag == 0  % use room receiver measurements.
    
    
GPStimeAMicroWeek = GPStimeAMicroWeekRoom;
GPStimeAMicro = GPStimeAMicroRoom;
ValidDataRinexObs = ValidDataRinexObsRoom;
L1_PRNAMicro  = L1_PRNAMicroRoom;
C1_PRNAMicro = C1_PRNAMicroRoom;   
    
    
    
    
end








for i = 1:NumSatEpochs  %epoch counter in seconds , number of points to calculate satellite positions for.
    
    for SV = 1:32
        
        %the pseudorange so that time of transmission can be calculated
        PR = C1_PRNAMicro(SV,i);
        
        SVTime_GPSSecsOrb(i) = SVTime_GPSSecs(i) - PR/c; %pass the time of transmission to the GPSOrbitPropagator

        [SV_X_Data(i,SV) SV_Y_Data(i,SV) SV_Z_Data(i,SV) SV_T_Data(i,SV) ValidData_Satellite(i,SV)] = GPSOrbitPropagator(SVWeek(i),  SVTime_GPSSecsOrb(i), SV, NavData);

        [SV_Xvel_Data(i,SV) SV_Yvel_Data(i,SV) SV_Zvel_Data(i,SV) SV_Tvel_Data(i,SV) SV_Xacc_Data(i,SV) SV_Yacc_Data(i,SV) SV_Zacc_Data(i,SV) SV_Tacc_Data(i,SV) ValidDataSatVels(i,SV)] = GPSOrbitPropagatorVelocities(SVWeek(i),  SVTime_GPSSecsOrb(i), SV, NavData);
       
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


NumEpochs = 4000; %number of epochs to use %note if Sp3 used then each epoch is 15 minutes.

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
Alarm_Limit = 100000%300; %metres ..for particular phase of flight.
ElevationMask = 0 ; %7.5 degrees Satellite elevation mask to use for satellite accept or reject criteria





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


% Generates clock bias and %drift for 1000 seconds. (this can be made longer or smaller in GARD_CLockModel())

[ClockBias_t ClockDrift_t ClockBias_f ClockDrift_f] = GARD_ClockModel; %time in seconds, frequencies in hertz


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



dTpos = 0; %initialise

dTpos  = zeros(1,1:NumEpochs);

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


for i = 1:NumEpochs


       
    %  variables to be cleared from workspace at each iteration
    clear  PRMeasured_Observed SVPos PRMeasured_FDE_Observed...  %PR Observed
           PRMeasured_Simulated SVPos PRMeasured_FDE_Simulated... %PR Simulated
           PRRate_Observed SVVel                                   %PR rates observed

    %-------------------------------------
    %ERROR CHECKING FLAGS

    RUNFLAG_Solution = 0;  %Solution can be made, (number of satellites is okay, run the code)
    RUNFLAG_FDE = 0;               %RAIM FDE can be done.
    %-------------------------------------


    %True User Position and clock error. (from previous estimate of clock error from clock model)
    Xpos(i) = -5046773.574;  %put roof antenna position in here and can verify that the nav. solution is correct or not.
    Ypos(i) = 2568446.555;
    Zpos(i) = -2925288.451;
    %dTpos(i) = 0;  %put this back in if dont' want to update the next solution with the clock
    
   
    
    
    
    Xvel(i) = 0;
    Yvel(i) = 0;
    Zvel(i) = 0;
    dTvel(i) = 0;
    
    
%     Xpos(i) = -5045790.38107811;
%      Ypos(i) = 2570291.64376626;
%       Zpos(i) = -2922762.10591193;
                           


%     Xpos(i) = 0;  %put roof antenna position in here and can verify that the nav. solution is correct or not.
%      Ypos(i) = 0;
%      Zpos(i) = 0;
    %Xpos(i) = 0;  %put roof antenna position in here and can verify that the nav. solution is correct or not.
    %Ypos(i) = 0;
    %Zpos(i) = 0;


    %dTpos(i) = ClockBias_t(i)*Speedoflight;  %in metres
  


    %Update next Clock bias estimate with previous estimate from observed
    %pseudorange solution
    %dTpos(i) = FinalSolution_Observed(i,4);


    %this is the real, True user position and clock bias  based on flight planning inputs

    UserPos = ([Xpos(i),Ypos(i),Zpos(i),dTpos(i)]); %this will be a copy of the values set by the user in the flight planning section


    UserVel = ([Xvel(i),Yvel(i),Zvel(i),dTvel(i)]);

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
        %check for which satellties exist in the roof receiver's file
    Valid = size(ValidDataRinexObsRoof);

    for PRN = 1:Valid(1)
        if ValidDataRinexObsRoof(PRN,i) == 1
            ValidData_SatelliteVectorRoof(i,PRN) = 1;
        else
            ValidData_SatelliteVectorRoof(i,PRN) = 0;
        end
        
    end
   
     
    
    %check for which satellties exist in the Room receivers file
   Valid = size(ValidDataRinexObsRoom);

    for PRN = 1:Valid(1)
        if ValidDataRinexObsRoom(PRN,i) == 1
            ValidData_SatelliteVectorRoom(i,PRN) = 1;
        else
            ValidData_SatelliteVectorRoom(i,PRN) = 0;
        end
    end
    


    
    %Choose which satellites to use which are common to both roof and room
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

    sel = sin(el);
    cel = cos(el);
    saz = sin(az);
    caz = cos(az);
    Telev(1,1) = -sel*caz;
    Telev(1,2) = -sel*saz;
    Telev(1,3) = cel;
    Telev(2,1) = -saz;
    Telev(2,2) = caz;
    Telev(2,3) = 0.0;
    Telev(3,1) = cel*caz;
    Telev(3,2) = cel*saz;
    Telev(3,3) = sel;

    nn = 0;

    %transform xyz to spherical coords and NEU

    
        
    Valid = size(ValidData_SatelliteVectorRoof);
    
    %Valid = size(ValidData_SatelliteVector);

    for Counter = 1:Valid(2)

        %if the satellite data is in the Observation file AND there is data
        %from the navigation file then:
        %if ValidData_SatelliteVector(i,Counter) == 1 && ValidData_Satellite(i,Counter) == 1  %the original

        %checks roof observation file, room observation file, and roof
        %navigation file (use a common nav file)
        if ValidData_SatelliteVectorRoof(i,Counter) == 1 & ValidData_SatelliteVectorRoom(i,Counter) == 1 & ValidData_Satellite(i,Counter) == 1 %for repeater testing
            
            GPS_Vec(Counter,1) = SV_X_Data(i,Counter) - XApos;
            GPS_Vec(Counter,2)= SV_Y_Data(i,Counter) - YApos;
            GPS_Vec(Counter,3) =	SV_Z_Data(i,Counter) - ZApos;
            
            topox(Counter) = Telev(1,1)*GPS_Vec(Counter,1) + Telev(1,2)*GPS_Vec(Counter,2) + Telev(1,3)*GPS_Vec(Counter,3);
            topoy(Counter) = Telev(2,1)*GPS_Vec(Counter,1) + Telev(2,2)*GPS_Vec(Counter,2) + Telev(2,3)*GPS_Vec(Counter,3);
            topoz(Counter) = Telev(3,1)*GPS_Vec(Counter,1) + Telev(3,2)*GPS_Vec(Counter,2) + Telev(3,3)*GPS_Vec(Counter,3);
            
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
        
        
        
         [SVPos,PRMeasured_Observed,PRMeasured_Simulated,CPMeasPredict,NAmb, GeometricRange, SVVel, PRRate_Observed] = GARDSim_GenerateMeasurements(i,GPStimeAMicro(i),UserPos,N, Sats,C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data, SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data, ALPHA, BETA);
        
                         
        

        %=================================================================
        %LEAST SQUARES and RAIM Solution from REAL Measurements
        %=================================================================

        %Do Least squares and RAIM FDE solution For Observed Pseudoranges
        %from Rinex Files

        [SolutionVec_Observed(i,1:4) VarSolutionVec_Observed NumIterations_Observed(i) ResVec_Observed M_Observed LSQ_Fail_Observed(i) limit_Observed(i)] = GARD_LSQ(UserPos,N,PRMeasured_Observed,SVPos);

      
        %Calculate GPS velocity from range rates which were calculated from
        %the carrier phase
                
        [SolutionVec_Observed_Vel(i,1:4) VarSolutionVec_Observed_Vel NumIterations_Observed_Vel(i) ResVec_Observed_Vel M_Observed_Vel LSQ_Fail_Observed_Vel(i) limit_Vel(i)] = GARD_LSQVel(SolutionVec_Observed(i,1:4),UserVel,N,PRRate_Observed,SVPos,SVVel);

        
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

           % [BadGeometry_Observed(i) RAIM_ALARM_Observed(i) SLOPE_Max_Observed(i) r_Observed(i) r_Observed_Threshold(i) ARP_Observed(i)] = GARD_RAIM(N,PFalseAlarm,SigmaS,Alarm_Limit,ResVec_Observed,M_Observed);

        else
%             BadGeometry_Observed(i) = 0;
%             RAIM_ALARM_Observed(i)= 0;
%             SLOPE_Max_Observed(i)= 0;
%             r_Observed(i)= 0;

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

        %==================================================================

        %Do RAIM FDE solution for Observed PRs

        %FDE taken from Parkinson Volume 2

        %if failure...go through all the available subsets...
        % if 6 satellites, can only detect one error at a time.

        % if error occurs on more than one subset then have to assume that more
        % than one satellite failure has occurred.


        %Do FDE solution
        %if failure:

%         if RUNFLAG_FDE == 0  %Run only if there's enough satellites for FDE (6 or more).
% 
%             if RAIM_ALARM_Observed(i) == 100 %If there's a RAIM alarm attempt FDE solution
% 
%                 %Do FDE - ONLY DETECTS ONE SATELLITE FAILURE AT A TIME AT THE MOMENT
% 
%                 %Form N-1 Subsets
%            
%               SimulatedOrObserved_Flag = 0; %Use observed PR measurements
%                 
%                 
%              [FinalSolution_FDE_Obs RAIM_ALARM_FDE_Obs GoodSatsSats_FDE_Obs Sats_FDE_Obs] = GARDSim_FDE(SimulatedOrObserved_Flag,i,GPStimeAMicro(i),UserPos,N, Sats,C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data,SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data,PFalseAlarm,SigmaS,Alarm_Limit,ALPHA,BETA);
%             
%              
%          
%              
%             FinalSolution_Observed(i,1) = FinalSolution_FDE_Obs(1);
%             FinalSolution_Observed(i,2) = FinalSolution_FDE_Obs(2);
%             FinalSolution_Observed(i,3) = FinalSolution_FDE_Obs(3);
%             FinalSolution_Observed(i,4) = FinalSolution_FDE_Obs(4);
%                
%       
% 
% 
%             else   %if there's no RAIM alarm
% %             GoodSatsSats_FDE_Obs(i,:) = 0;
% %             RAIM_ALARM_FDE_Obs(i,:) = 0;
%             
%             end  %If there's a RAIM alarm attempt FDE solution
% 
% 
%         end %RUNFLAG_FDE == 0   %End FDE Section for observed measurements






        %=================================================================
        %LEAST SQUARES and RAIM Solution from SIMULATED Measurements
        %=================================================================


        [SolutionVec_Simulated(i,1:4) VarSolutionVec_Simulated NumIterations_Simulated(i) ResVec_Simulated M_Simulated LSQ_Fail_Simulated(i) limit_Simulated(i)] = GARD_LSQ(UserPos,N,PRMeasured_Simulated,SVPos);

       
        if LSQ_Fail_Simulated(i) == 0   %If solution has converged

            FinalSolution_Simulated(i,1) = SolutionVec_Simulated(i,1);
            FinalSolution_Simulated(i,2) = SolutionVec_Simulated(i,2);
            FinalSolution_Simulated(i,3) = SolutionVec_Simulated(i,3);
            FinalSolution_Simulated(i,4) = SolutionVec_Simulated(i,4);

            %Do RAIM FD solution using least squares residuals method

            %RAIM_ALARM == 100 = true
            %RAIM_ALARM == Anything else = False

            %need to implement parity method, and parity isolation method. cf
            %parkinson blue book.

            %[BadGeometry_Simulated(i) RAIM_ALARM_Simulated(i) SLOPE_Max_Simulated(i) r_Simulated(i) r_Simulated_Threshold(i) ARP_Simulated(i)] = GARD_RAIM(N,PFalseAlarm,SigmaS,Alarm_Limit,ResVec_Simulated,M_Simulated);

            
        else

%             BadGeometry_Simulated(i) = 0;
%             RAIM_ALARM_Simulated(i)= 0;
%             SLOPE_Max_Simulated(i)= 0;
%             r_Simulated(i)= 0;

            FinalSolution_Simulated(i,1) = 0;
            FinalSolution_Simulated(i,2) = 0;
            FinalSolution_Simulated(i,3) = 0;
            FinalSolution_Simulated(i,4) = 0;
            
            


        end

        %==================================================================

     %Do FDE solution
        %if failure:

%         if RUNFLAG_FDE == 0  %Run only if there's enough satellites for FDE (6 or more).
% 
%             if RAIM_ALARM_Simulated(i) == 100 %If there's a RAIM alarm attempt FDE solution
% 
%                 %Do FDE - ONLY DETECTS ONE SATELLITE FAILURE AT A TIME AT THE MOMENT
% 
%                 %Form N-1 Subsets
% 
%              SimulatedOrObserved_Flag = 1; %Use simulated PR measurements
%                 
%              [FinalSolution_FDE_Sim RAIM_ALARM_FDE_Sim GoodSatsSats_FDE_Sim Sats_FDE_Sim] = GARDSim_FDE(SimulatedOrObserved_Flag,i,GPStimeAMicro(i),UserPos,N, Sats,C1_PRNAMicro,L1_PRNAMicro, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data,SV_Xvel_Data, SV_Yvel_Data, SV_Zvel_Data, SV_Tvel_Data,PFalseAlarm,SigmaS,Alarm_Limit,ALPHA,BETA);
%                
%                              
%             FinalSolution_Simulated(i,1) = FinalSolution_FDE_Sim(1);
%             FinalSolution_Simulated(i,2) = FinalSolution_FDE_Sim(2);
%             FinalSolution_Simulated(i,3) = FinalSolution_FDE_Sim(3);
%             FinalSolution_Simulated(i,4) = FinalSolution_FDE_Sim(4);
%                
%              
%              
%              
%             else   %if there's no RAIM alarm
% %             GoodSatsSats_FDE_Sim(i,:) = 0;   
% %             RAIM_ALARM_FDE_Sim(i,:) = 0;
%             
%             end  %If there's a RAIM alarm attempt FDE solution
% 
% 
%         end %RUNFLAG_FDE == 0   %End FDE Section for observed measurements



        %Update next Clock bias estimate with previous estimate from observed
        %pseudorange solution
        %dTpos(i+1) = FinalSolution_Observed(i,4);


        %Take copy of PR's
        jjjjj = size(PRMeasured_Simulated);  %can use simulated or observed here..just to get the array size


        %Copy any local 'single epoch' variables into Global Variables for later analysis/plotting
        for pp = 1:jjjjj(2)
            PRMeasured_SimulatedGlobal(i,pp) = PRMeasured_Simulated(pp);
            PRMeasured_ObservedGlobal(i,pp) = PRMeasured_Observed(pp);
            PRRate_ObservedGlobal(i,pp) = PRRate_Observed(pp);
            
            %global copy of satellites used
            SatsGlobal(i,pp) = Sats(pp);
           
            %PreviousSmoothedGlobal(i,pp) = PreviousSmoothed(pp); %for smoothed PR solution
        end
        
        
        
        
        
        dTpos(i+1) = dTpos(i) + FinalSolution_Observed(i,4); %(updates the clock model at each epoch, take this out if want to calculate on epoch by epoch basis
        
        
            
        
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

TruthX = -5046773.57;
TruthY = 2568446.555;
TruthZ = -2925288.451;
TruthdT = 0;


Xerr_Observed = TruthX - FinalSolution_Observed(:,1);
Yerr_Observed = TruthY - FinalSolution_Observed(:,2);
Zerr_Observed = TruthZ - FinalSolution_Observed(:,3);
dTerr_Observed = TruthdT - FinalSolution_Observed(:,4);


Xerr_Simulated = TruthX - FinalSolution_Simulated(:,1);
Yerr_Simulated = TruthY - FinalSolution_Simulated(:,2);
Zerr_Simulated = TruthZ - FinalSolution_Simulated(:,3);
dTerr_Simulated = TruthdT - FinalSolution_Simulated(:,4);



plot(Zerr_Observed)
hold
plot(Zerr_Simulated,'r')

%difference between measured and predicted PR's
PRresult = PRMeasured_SimulatedGlobal-PRMeasured_ObservedGlobal;

results = PRresult(:,1) - ClockBias_t(1:1000)'  ;

aaaa = results - dTer_Simulated;
%Subtract clock bias



blahblah = PRresult(:,1) + FinalSolution_Observed(:,4);


%see p201 Kayton, HDOP or VDOP needs to be converted to local tangent plane
%(NEU) before use.




