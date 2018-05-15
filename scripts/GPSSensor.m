
%======================================================================
%GPS SENSOR
%======================================================================

%generate TCXO clock errors for the GPS.

TimeInterval = 1; %1 second

NumEpochs = (endepoch - startepoch)+10;

Delta_t = TimeInterval; %this is the time between samples in seconds
%Numberpoints = floor(NumEpochs/Delta_t); %this is the number of sample points %Num epochs is really in seconds.
%Numberpoints = floor(NumEpochs*TimeInterval/Delta_t);
Numberpoints = NumEpochs;

load 'data/WhiteNoiseForClockCoast/WhiteNoise.mat';  %load 100,000 points of white noise.

U1t = U1(1:Numberpoints);
U2t = U2(1:Numberpoints);
%[dTpos(i,:) dTposdot(i,:) dTposflicker(i,:) dTposdotflicker(i,:)] = GARD_ClockModelKF(Numberpoints,Delta_t);
[dTposTCXO dTposdotTCXO dTposCSAC dTposdotCSAC] = GARD_ClockModelKF(Numberpoints,Delta_t,U1t,U2t); %note this is in seconds!!


% %set to zero, no clock error, for now.
% dTposTCXO = zeros(1,NumEpochs);
% dTposdotTCXO = zeros(1,NumEpochs);


%calculates GPS measurements at 1 Hz

%dont forget to run GPSConstants above before running this code by itself.
FilenameOptimised = 'GPSConstellationOptimised.txt';
%read optimised constellation data from file
[NavData, TruePosVec] = freadOptimisedNav(FilenameOptimised);

GPSWeek = 829;
GPSSec = 432000;



%pos truth for GPS

%simulate rate in hz , can't be greater than 50, must be a whole number
Rate_GPS = 1;

%dt_GPS = 50/Rate_GPS;

a = 1; %change the rate for this next section

Pseudoranges = zeros(1,16); %allow for up to 16 satellites
PseudorangeRates = zeros(1,16); %allow for up to 16 satellites

%This bit represents the GPS receiver
for i = startepoch-1:endepoch  %start the GPS at epoch 1 so the pr rates and satellite vels are not a large number or zero for epoch 2

    %Xpos_truth(i*a-(a-1))
    %initial guess start point
    PosTruthGPS(i,1:4) = [Xpos_truth1Hz(i), Ypos_truth1Hz(i), Zpos_truth1Hz(i), c*dTposTCXO(i)];
    VelTruthGPS(i,1:4) = [ Xvel_truth1Hz(i), Yvel_truth1Hz(i), Zvel_truth1Hz(i), c*dTposdotTCXO(i)];

    %this function spits out pseudoranges and calculated satellites at each
    %epoch and the receivers calculated position.

    %starting point DECEMBER 1, 1995 AT 00:00 UTC
    
    ElevationMask = 0; %degrees
    
%     if i >5
%         
%            ElevationMask = 0; %degrees
%            
%     end
%     
%     if i > 10
%         
%         ElevationMask = 20;
%     end
%         
%         if i > 12
%             ElevationMask = 5;
%             
%         end
                
        
        %calculate centre of antenna to satellite vector
        
        %rotate plane with the true roll pitch and yaw of aircraft
                
        %if the vector cuts through the plane then signal is blocked, so
        %block that particular measurement.    
        
                
    
    %put code here for antenna mask angle algorithm which simulates
    %satellites dropping in and out with aircraft bank angle.

    %this function generates 

    [GPSPosLSQTemp, SatPosTemp,SatVelTemp,SatellitesTemp,PseudorangesTemp, PseudorangeRatesTemp,NTemp, Sat_PRN_VecTemp] = SimulatedGPS(PosTruthGPS(i,1:4), VelTruthGPS(i,1:4), GPSWeek,GPSSec,ElevationMask,NavData,Roll_truth1Hz(i), Pitch_truth1Hz(i), Yaw_truth1Hz(i) );
            
     
    
    
    
    

    GPSPosLSQ(i,:) = GPSPosLSQTemp;    
    SatPos(i,:,:) = SatPosTemp;  %where i is the epoch, : is the index corresponding to the satellites PRN number and : is the x, y, z or dT
    SatVel(i,:,:) = SatVelTemp;
       
    
    
    uuuu = size(PseudorangesTemp);
    PRVecSize = uuuu(2);
    
    Satellites(i,1:NTemp) = SatellitesTemp;   
    Pseudoranges(i,1:PRVecSize)  = PseudorangesTemp;
    PseudorangeRates(i,1:PRVecSize) = PseudorangeRatesTemp;
    N_save(i) = NTemp;   
    
    Sat_PRN_Vec(i,:) = Sat_PRN_VecTemp;
   
    
    
   
    


       

    Xpos_GPS(i) = GPSPosLSQ(i,1);
    Ypos_GPS(i) = GPSPosLSQ(i,2);
    Zpos_GPS(i) = GPSPosLSQ(i,3);
     dt_GPS(i) = GPSPosLSQ(i,4);

    %GPS velocities
    if i == startepoch-1  %BE CAREFUL, i didn't have this as -1, so it was giving hte wrong velocity. Because its startepoch-1 above.

        Xvel_GPS(i) = (Xpos_GPS(i) - Xpos_truth(i*50))/1;  %51 at 50 Hz is i = 1 at 1 Hz
        Yvel_GPS(i) = (Ypos_GPS(i) - Ypos_truth(i*50))/1;
        Zvel_GPS(i) = (Zpos_GPS(i) - Zpos_truth(i*50))/1;
         dtvel_GPS(i) = 0;

    else
        Xvel_GPS(i) = (Xpos_GPS(i) - Xpos_GPS(i-1))/1;
        Yvel_GPS(i) = (Ypos_GPS(i) - Ypos_GPS(i-1))/1;
        Zvel_GPS(i) = (Zpos_GPS(i) - Zpos_GPS(i-1))/1;
         dtvel_GPS(i) = (dt_GPS(i) - dt_GPS(i-1))/TimeInterval;

    end
    
    
    %NOTE
   % note that diff(Xpos_GPS(49:99)) equals Xvel_GPS(50:100)



    GPSSec = GPSSec+1;

    %check for GPS week rollover

    if GPSSec > 604800
        GPSWeek = GPSWeek+1;
        GPSSec = 0;
    end
end

%estimate satellite velocities by differentiating (theres ap roblem
%with the sat. vel propagator code.


drift = 0;


%======================================================================
%END GPS SENSOR
%======================================================================