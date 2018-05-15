%Use state 0 of the random number generator

randn('state',0);




i = startepoch; % to start it

Temp_noise(i-1) = Sigma_Temp_noise*randn(1);
Alpha_noise(i-1) = Sigma_Alpha_noise*randn(1);
Beta_noise(i-1) = Sigma_Beta_noise*randn(1);
Elevator_noise(i-1) = Sigma_Elevator_noise*randn(1);
Aileron_noise(i-1) = Sigma_Aileron_noise*randn(1);
Rudder_noise(i-1) = Sigma_Rudder_noise*randn(1);
Throttle_noise(i-1) = Sigma_Throttle_noise*randn(1);
Altimeter_noise(i-1) = Sigma_Altimeter_noise*randn(1);


Airspeed_noise(i-1) = Sigma_Airspeed_noise*randn(1);


    CL_noise(i-1) = Sigma_CL_noise*randn(1);
    CD_noise(i-1) = Sigma_CD_noise*randn(1);
    CY_noise(i-1) = Sigma_CY_noise*randn(1);
    Cl_noise(i-1) = Sigma_Cl_noise*randn(1);
    Cm_noise(i-1) = Sigma_Cm_noise*randn(1);
    Cn_noise(i-1) = Sigma_Cn_noise*randn(1);
            
    
  
    %add satellite orbit error
%     
%     for isat_count = 1:32       
%         
%      %GPSPr_noise(isat_count,i) = Sigma_PR_noise*randn(1); 
%      
%      
%       SatOrbError(isat_count,1,i-1) = 30*randn(1);   %metres
%       SatOrbError(isat_count,2,i-1) = 30*randn(1);
%        SatOrbError(isat_count,3,i-1) = 30*randn(1);
% %     
%    
%     end
%     
           
    
   
      
    
    %this is 1 sigma standard deviations all in metres
    Sigma_Ephem_noise = 3;
    
%     
%     
    Sigma_Iono_noise = 5;    
    Sigma_Tropo_noise = 2;    
    
    
%     %for L1-L5
%     Sigma_Iono_noise = 0;    
%     Sigma_Tropo_noise = 0;    
    
    
    
    Sigma_Multipath_noise = 5;     
    Sigma_ReceiverMeasurement_noise = 3;
    
    
    
    Sigma_Multipath_noise_carr = 0.048; 
    Sigma_ReceiverMeasurement_noise_carr = 0.0019;
    
    %Beta = 1/tau (1/seconds)
    
    Beta_Ephem_noise = 1/1800;
    Beta_Iono_noise = 1/1800;    
    Beta_Tropo_noise = 1/3600;    
    Beta_Multipath_noise = 1/600;     
    Beta_ReceiverMeasurement_noise = 1000000;  %large because its white noise
     
    
    Beta_Multipath_noise_carr = 1/600;
    Beta_ReceiverMeasurement_noise_carr = 1000000;

    
%     
%     
%     Sigma_Multipath_noise = 0;     
%     Sigma_Multipath_noise_carr = 0 ;
%     Sigma_Ephem_noise = 0;
    
     
     
     %should MP be absolute or +ve or -ve? I thought +ve or -ve..
     

     
      %GPS Noises
      
      for isat_count = 1:32
    i = startepoch;
    GPS_Ephem_noise(isat_count,i-1) = abs(Sigma_Ephem_noise*randn(1));
    GPS_Iono_noise(isat_count,i-1) = abs(Sigma_Iono_noise*randn(1));
    GPS_Tropo_noise(isat_count,i-1) = abs(Sigma_Tropo_noise*randn(1));
    GPS_Multipath_noise(isat_count,i-1) = abs(Sigma_Multipath_noise*randn(1));
    GPS_ReceiverMeasurement_noise(isat_count,i-1) = abs(Sigma_ReceiverMeasurement_noise*randn(1));
           
    
    GPS_Multipath_noise_carr(isat_count,i-1) = abs(Sigma_Multipath_noise_carr*randn(1));
    GPS_ReceiverMeasurement_noise_carr(isat_count,i-1) = Sigma_ReceiverMeasurement_noise_carr*randn(1);
     
    PrCorrNoiseCheck(isat_count,i-1) = 7.94*randn(1); 
    
%     
%      GPS_Ephem_noise(isat_count,i-1) = 0;
%     GPS_Iono_noise(isat_count,i-1) = 0;
%     GPS_Tropo_noise(isat_count,i-1) = 0;
%     GPS_Multipath_noise(isat_count,i-1) = 0;
%     GPS_ReceiverMeasurement_noise(isat_count,i-1) = 0;
%            
%     
%     GPS_Multipath_noise_carr(isat_count,i-1) = 0;
%     GPS_ReceiverMeasurement_noise_carr(isat_count,i-1) = 0;
    
    
     
    %GPSPr_noise(isat_count,i-1) = 0;
     
      end
      
      
    
    for isat_count = 1:32
        
        
        for i = startepoch:endepoch
                        
         %for i = startepoch:10000                           
            
            %GPS C/A model for PR's from James Rankin "GPS and differential
            %GPS: An Error Model for Sensor Simulation
            %assumed it is standard C/A correlator (not narrow)
            
            %also see brown and hwang section 11.3 gps error models
            
            %=====================
            %CORRELATED NOISE
            %=====================
            
            
            % ephemeris noise 
            [GPS_Ephem_noise(isat_count,i),variance,randU] = GaussMarkov_Process(GPS_Ephem_noise(isat_count,i-1),Beta_Ephem_noise,Sigma_Ephem_noise,1);
           
            
            
                % iono noise            
            [GPS_Iono_noise(isat_count,i),variance,randU] = GaussMarkov_Process(GPS_Iono_noise(isat_count,i-1),Beta_Iono_noise,Sigma_Iono_noise,1);
            
              % tropo noise            
           [GPS_Tropo_noise(isat_count,i),variance,randU] = GaussMarkov_Process(GPS_Tropo_noise(isat_count,i-1),Beta_Tropo_noise,Sigma_Tropo_noise,1);             
                                       
            
              %Multipath C/A standard
              [GPS_Multipath_noise(isat_count,i),variance,randU] = GaussMarkov_Process(GPS_Multipath_noise(isat_count,i-1),Beta_Multipath_noise,Sigma_Multipath_noise,1);
        
               
              %Multipath L1 Carrier
              [GPS_Multipath_noise_carr(isat_count,i),variance,randU] = GaussMarkov_Process(GPS_Multipath_noise_carr(isat_count,i-1),Beta_Multipath_noise_carr,Sigma_Multipath_noise_carr,1);
                        
                  
              
              
              [PrCorrNoiseCheck(isat_count,i),variance,randU] = GaussMarkov_Process(PrCorrNoiseCheck(isat_count,i-1),1/1100,7.94,1);
                        
                  
                          
                  
     %=======================     
       %MEASUREMENT NOISE
       %=======================   
       %i'm assuming this is different for every PR measurement instead of
       %common to all PR measurements. 
                
       
       %this is gaussian white noise, so Beta is set large 
              [GPS_ReceiverMeasurement_noise(isat_count,i),variance,randU] = GaussMarkov_Process(GPS_ReceiverMeasurement_noise(isat_count,i-1),Beta_ReceiverMeasurement_noise,Sigma_ReceiverMeasurement_noise,1);
    
    
              %for the carrier phase
               [GPS_ReceiverMeasurement_noise_carr(isat_count,i),variance,randU] = GaussMarkov_Process(GPS_ReceiverMeasurement_noise_carr(isat_count,i-1),Beta_ReceiverMeasurement_noise_carr,Sigma_ReceiverMeasurement_noise_carr,1);
    
              
                            
         %add up the errors      
          %GPSPr_noise(i) =  sqrt(GPS_Ephem_noise(i)^2 + GPS_Iono_noise(i)^2 + GPS_Tropo_noise(i)^2 + GPS_Multipath_noise(i)^2);% + GPS_ReceiverMeasurement_noise(i);
              
            GPSPr_noise(isat_count,i) =  abs((GPS_Ephem_noise(isat_count,i)) + abs(GPS_Iono_noise(isat_count,i)) + abs(GPS_Tropo_noise(isat_count,i)) + abs(GPS_Multipath_noise(isat_count,i))) + GPS_ReceiverMeasurement_noise(isat_count,i);
                  
              
             GPSCarrier_noise(isat_count,i) =  abs((GPS_Ephem_noise(isat_count,i)) - abs(GPS_Iono_noise(isat_count,i)) +  abs(GPS_Tropo_noise(isat_count,i)) + abs(GPS_Multipath_noise_carr(isat_count,i))) + GPS_ReceiverMeasurement_noise_carr(isat_count,i);
             
             
             %plotting this will give standard deviation of about 20cm
             %which is typical of GPS PR Rate I think
             GPSPr_Rate_noise(isat_count,i) =  GPSCarrier_noise(isat_count,i) -  GPSCarrier_noise(isat_count,i-1);
      
                
             
              CorrCheck(isat_count,i) =  abs((GPS_Ephem_noise(isat_count,i)) + abs(GPS_Iono_noise(isat_count,i)) + abs(GPS_Tropo_noise(isat_count,i)) + abs(GPS_Multipath_noise(isat_count,i)));% + GPS_ReceiverMeasurement_noise(isat_count,i);
                  
              
             
%               
%        %generate one noise for each pseudorange
%        
%        
%        
%        
%    
%      [GPSPr_noise(isat_count,i+1)] = AutoRegressive_Process(GPSPr_noise(isat_count,i),0.998,Sigma_PR_noise);% + (5+rand(1)*5); %this bit on the end represents iono error %metres
%      
% %      SatOrbError(isat_count,1,i) = 3.5*randn(1);   %metres
% %       SatOrbError(isat_count,2,i) = 3.5*randn(1);
% %        SatOrbError(isat_count,3,i) = 3.5*randn(1);
% 


       end
           
     
    end
    
    
    
    
    
       
    
%     
%     
%     
% 
% maxlags = 10000;
% 
% [autocorr1,lags] = xcorr(CorrCheck(6,:),maxlags);
% 
% autocorr1p = autocorr1/length(xcorr(CorrCheck(6,:)));
% 
% %autocorr1p = autocorr1;
% 
% plot(lags, autocorr1p);
% 
% 
% 
% 
% maxlags = 10000;
% [autocorr1,lags] = xcorr(PrCorrNoiseCheck(6,:),maxlags);
% 
% autocorr1p = autocorr1/length(xcorr(PrCorrNoiseCheck(6,:)));
% 
% %autocorr1p = autocorr1;
% 
% 
% plot(lags, autocorr1p,'r');
%     
    

% for i = startepoch:endepoch
%     %========================================================================
%     %AIR DATA SENSORS, at 1 Hz
%     %========================================================================
% 
% 
%     %-----------------------------------------------------
%     %Temperature Sensor
%     %-----------------------------------------------------
%     %in Degrees C
%     
     dT = 1;
%     
%    [Temp_noise(i),variance,randU] = GaussMarkov_Process(Temp_noise(i-1),Beta_Temp_noise,Sigma_Temp_noise,dT );
%    Temp_SensorTruth(i)  = AtmosTruth(2,i*100+1);
%    Temp_Sensor(i)  = Temp_SensorTruth(i) + Temp_noise(i);
%         
%        
%    %GPS PR Noise
%   % for i = startepoch:endepoch+8000
%    
%       %[GPSPr_noise(i)] = GaussMarkov_Process(GPSPr_noise(i-1),Beta_PR_noise,Sigma_PR_noise);
%     %[GPSPr_noise(i)] = GaussMarkov_Process(GPSPr_noise(i-1),5,Sigma_PR_noise);
%    %end
%    
%    
%    
%       
%    
%    
% 
%     
% 
%    
%    
%    
% 
%     %-----------------------------------------------------
%     %Static Pressure Sensor
%     %-----------------------------------------------------
% 
% 
%     %-----------------------------------------------------
%     %Dynamic Pressure Sensor
%     %-----------------------------------------------------
% 
% 
%     %-----------------------------------------------------
%     %Airspeed Sensor
%     %-----------------------------------------------------  
%     
%     [Airspeed_noise(i),variance,randU] =  GaussMarkov_Process(Airspeed_noise(i-1),Beta_Airspeed_noise,Sigma_Airspeed_noise,dT);
%    Airspeer_SensorTruth(i) = Airspeed_truth(i*100+1);
%     Airspeed_Sensor(i) = Airspeer_SensorTruth(i) + Airspeed_noise(i) ;   %round quantizes it to a metre
%     
%     
%     
%       
% 
% 
%   
%     
%     %-----------------------------------------------------
%     %Altimeter Sensor
%     %-----------------------------------------------------
% 
%     %in Metres
%     [Altimeter_noise(i),variance,randU] =  GaussMarkov_Process(Altimeter_noise(i-1),Beta_Altimeter_noise,Sigma_Altimeter_noise,dT);
%     Altimeter_SensorTruth(i) = Hgt_truth1Hz(i);
%     Altimeter_Sensor(i) = round(Altimeter_SensorTruth(i) + Altimeter_noise(i)) ;   %round quantizes it to a metre
%     
%   
%    
%     
% 
%     %-----------------------------------------------------
%     %Angle of Attack Sensor
%     %-----------------------------------------------------
%     %in radians
%     
%     [Alpha_noise(i),variance,randU] =  GaussMarkov_Process(Alpha_noise(i-1),Beta_Alpha_noise,Sigma_Alpha_noise,dT);
%     Alpha_SensorTruth(i) = Alpha_truth(i*100+1);
%     Alpha_Sensor(i) = Alpha_SensorTruth(i) + Alpha_noise(i) ;
% 
%     %-----------------------------------------------------
%     %Angle of Sideslip Sensor
%     %-----------------------------------------------------
%     [Beta_noise(i),variance,randU] =  GaussMarkov_Process(Beta_noise(i-1),Beta_Beta_noise,Sigma_Beta_noise,dT);
%     Beta_SensorTruth(i) = Beta_truth(i*100+1);
%     Beta_Sensor(i) = Beta_SensorTruth(i) + Beta_noise(i);
% 
% 
%     %========================================================================
%     %Control Surface Sensors
%     %========================================================================
% 
%     [Elevator_noise(i),variance,randU] = GaussMarkov_Process(Elevator_noise(i-1),Beta_Elevator_noise,Sigma_Elevator_noise,dT);
%     [Aileron_noise(i),variance,randU] = GaussMarkov_Process(Aileron_noise(i-1),Beta_Aileron_noise,Sigma_Aileron_noise,dT);
%     [Rudder_noise(i),variance,randU] = GaussMarkov_Process(Rudder_noise(i-1),Beta_Rudder_noise,Sigma_Rudder_noise,dT);
%     [Throttle_noise(i),variance,randU] = GaussMarkov_Process(Throttle_noise(i-1),Beta_Throttle_noise,Sigma_Throttle_noise,dT);
%     
%     
%     
%     %Aircraft parameter noise in %
%     
%     [CL_noise(i),variance,randU] = GaussMarkov_Process(CL_noise(i-1),Beta_CL_noise,Sigma_CL_noise,dT);     
%       [CD_noise(i),variance,randU] = GaussMarkov_Process(CD_noise(i-1),Beta_CD_noise,Sigma_CD_noise,dT); 
%        [CY_noise(i),variance,randU] = GaussMarkov_Process(CY_noise(i-1),Beta_CY_noise,Sigma_CY_noise,dT); 
%        [Cl_noise(i),variance,randU] = GaussMarkov_Process(Cl_noise(i-1),Beta_Cl_noise,Sigma_Cl_noise,dT); 
%         [Cm_noise(i),variance,randU] = GaussMarkov_Process(Cm_noise(i-1),Beta_Cm_noise,Sigma_Cm_noise,dT); 
%         [ Cn_noise(i),variance,randU] = GaussMarkov_Process(Cn_noise(i-1),Beta_Cn_noise,Sigma_Cn_noise,dT); 
%     
%     
% 
%     
%     % 
% % %add errors
% %  CL = CL + 2*CL;% +  rand(1);
% %  
% %  CD = CD - 0.9*CD +  0.001*rand(1);
% %  CY = CY - 0.9*CY +  0.001*rand(1);
% %  Cm = Cm + 0.1*Cm +  0.001*rand(1);
% %  Cl = Cl - 0.1*Cl +  0.001*rand(1);
% %  Cn = Cn + 0.1*Cn +  0.001*rand(1);
% 
%                                                     
% 
% end

% Sigma_CL_noise = 0.068;   %this is to give roughly +/- 10% of true value , 0.068 is standard deviation of course for gaussian distribution
% Beta_CL_noise = 1/120;

Beta_Param_noise = 1/120;

Sigma_Param_noise =  0.068; %1 sigma value for 10% ie 0.1

 ParamNoise(1:33,startepoch-1) = Sigma_Param_noise*randn(1);

%for individual 10% uncertainty 



for i = startepoch:endepoch


for k = 1:33

[ParamNoise(k,i),variance,randU] = GaussMarkov_Process(ParamNoise(k,i-1),Beta_Param_noise,Sigma_Param_noise,dT);     
     
end    
        

end











% 
% XGyroBiasTruth = 0.1;   %in rad/s
% YGyroBiasTruth = -0.15;  %in rad/s
% ZGyroBiasTruth = 0.2;   %in rad/s
% 
% XAccelBiasTruth = 1;   %in m/s^2
% YAccelBiasTruth = 2;  %in m/s^2
% ZAccelBiasTruth = 0.5;  %in m/s^2


i = startepochHighRate;



An1Bias = 0.5;    %m/s^2
An2Bias = 0.45;   %m/s^2
An3Bias = -0.2;  %m/s^2


Gn1Bias = 0.57*pi/180;  %rad/s
Gn2Bias = 0.40*pi/180;   %rad/s
Gn3Bias = 0.38*pi/180;   %rad/s


% 
% An1Bias = 0;    %m/s^2
% An2Bias = 0;   %m/s^2
% An3Bias = 0;  %m/s^2
% 
% 
% Gn1Bias = 0;  %rad/s
% Gn2Bias = 0;   %rad/s
% Gn3Bias = 0;   %rad/s




%for putting ins fault on:

% 
% Gn1Bias = 0.1;  %rad/s
% Gn2Bias = -0.15;   %rad/s
% 
% Gn3Bias = 0.2;   %rad/s

% 
% 
% Gn1Bias = 0;  %rad/s
% Gn2Bias = 0;   %rad/s
% Gn3Bias = 0;   %rad/s

%0.0124

%starting point, I add the large initial bias. 
%dtINS = 0.01



Sigma_p = 0.00925;
Sigma_q = 0.00785;
Sigma_r = 0.00768;

Sigma_Ax = 0.013;
Sigma_Ay = 0.018;
Sigma_Az = 0.01;


% 
% 
% 
% Sigma_p = 0.00925/100;
% Sigma_q = 0.00785/100;
% Sigma_r = 0.00768/100;
% 
% Sigma_Ax = 0.013/100;
% Sigma_Ay = 0.018/100;
% Sigma_Az = 0.01/100;


dtINS = 0.01;

% if i > 2000  && i < 5000
% AN(1,i-1) = 0.0124*randn(1)*dtINS + An1Bias;
% AN(2,i-1) = 0.0124*randn(1)*dtINS + An2Bias;
% AN(3,i-1) = 0.0124*randn(1)*dtINS + An3Bias;
% 
% % AN(1,i-1) = 1*randn(1)*dtINS + An1Bias;
% % AN(2,i-1) = 1*randn(1)*dtINS + An2Bias;
% % AN(3,i-1) = 1*randn(1)*dtINS + An3Bias;
% 
% GN(1,i-1) = 0.005*randn(1)*dtINS + Gn1Bias; 
% GN(2,i-1) = 0.005*randn(1)*dtINS + Gn2Bias;
% GN(3,i-1) = 0.005*randn(1)*dtINS + Gn3Bias;
% else
AN(1,i-1) = Sigma_Ax*randn(1)*dtINS ;
AN(2,i-1) = Sigma_Ay*randn(1)*dtINS  ;
AN(3,i-1) = Sigma_Az*randn(1)*dtINS  ;

% AN(1,i-1) = 1*randn(1)*dtINS + An1Bias;
% AN(2,i-1) = 1*randn(1)*dtINS + An2Bias;
% AN(3,i-1) = 1*randn(1)*dtINS + An3Bias;

GN(1,i-1) = Sigma_p*randn(1)*dtINS  ; 
GN(2,i-1) = Sigma_q*randn(1)*dtINS  ;
GN(3,i-1) = Sigma_r*randn(1)*dtINS  ;



betaINSNoise = 100000;
% end
    
    

%  
% AN(1,i-1) = 0;
% AN(2,i-1) = 0;
% AN(3,i-1) = 0;
% 
% GN(1,i-1) = 0;
% GN(2,i-1) = 0;
% GN(3,i-1) = 0;


for i = startepochHighRate:endepochHighRate
    
    
    %for i = 2:200
%for i = 2:100
        
   %INS Sensor Noise 
    %Accel    
    %dtINS = 0.01;
 
    %NOTE THAT I INPUT BETA TO THE FUNCTION WHICH IS 1/tau = 1/300 sec
    
%     [AN(1,i), hhh] = GaussMarkov_Processfortesting(AN(1,i-1),0,1/300,0.0124,dtINS );
%     [AN(2,i),hhh] = GaussMarkov_Processfortesting(AN(2,i-1),0,1/300,0.0124,dtINS );
%     [AN(3,i),hhh] = GaussMarkov_Processfortesting(AN(3,i-1),0,1/300,0.0124,dtINS );
%     
%     [GN(1,i),hhh] = GaussMarkov_Processfortesting(GN(1,i-1),0,1/300,0.005,dtINS );
%     [GN(2,i),hhh] = GaussMarkov_Processfortesting(GN(2,i-1),0,1/300,0.005,dtINS );
%     [GN(3,i),hhh] = GaussMarkov_Processfortesting(GN(3,i-1),0,1/300,0.005,dtINS );
%     
%            
     
        
    [AN(1,i) ,variance,randU] = GaussMarkov_Process(AN(1,i-1),betaINSNoise,Sigma_Ax,dtINS );
    [AN(2,i),variance,randU] = GaussMarkov_Process(AN(2,i-1),betaINSNoise,Sigma_Ay,dtINS );
    [AN(3,i),variance,randU] = GaussMarkov_Process(AN(3,i-1),betaINSNoise,Sigma_Az,dtINS );
    
    [GN(1,i),variance,randU] = GaussMarkov_Process(GN(1,i-1),betaINSNoise,Sigma_p,dtINS );
    [GN(2,i),variance,randU] = GaussMarkov_Process(GN(2,i-1),betaINSNoise,Sigma_q,dtINS );
    [GN(3,i),variance,randU] = GaussMarkov_Process(GN(3,i-1),betaINSNoise,Sigma_r,dtINS );
    
    
%     [AN(1,i) ,variance,randU] = GaussMarkov_Process(AN(1,i-1),1/300,0,dtINS );
%     [AN(2,i),variance,randU] = GaussMarkov_Process(AN(2,i-1),1/300,0,dtINS );
%     [AN(3,i),variance,randU] = GaussMarkov_Process(AN(3,i-1),1/300,0,dtINS );
%     
%     [GN(1,i),variance,randU] = GaussMarkov_Process(GN(1,i-1),1/300,0,dtINS );
%     [GN(2,i),variance,randU] = GaussMarkov_Process(GN(2,i-1),1/300,0,dtINS );
%     [GN(3,i),variance,randU] = GaussMarkov_Process(GN(3,i-1),1/300,0,dtINS );
%     
           
        
    %add spikes in data, due to sudden jolts in the aircraff
          
    %need to do this on the second though to see any effect. 
%     if i == 650 || i == 660 || i == 800
%         
%      p_INS_50Hz(i) =  p_INS_50Hz(i) + 1;
%      
%       ax_b_INS_50Hz(i) =  ax_b_INS_50Hz(i) + 1;
%      
%     end
    
        
    %add slowly variation depending upon temperature, to simulate real INS.
             
                    
    
end




%need to scale noise by 1/dtINS so thatts in rad/s not 100ths of rad/s
 AN = AN*100;
GN = GN*100;
    




% 
Gn1Fault = 0.1;  %rad/s
Gn2Fault = -0.15;   %rad/s

Gn3Fault = 0.2;   %rad/s

for i = startepochHighRate:endepochHighRate


    
    
%%use this code to add a INS fault  
% 
% %     
% if i > 7500  && i < 8500
% % 
% %   %turn ins noise off while i debug aero model
%     ax_b_INS_50Hz(i) = AccelTruth100Hz(1,i) + AN(1,i) + An1Bias;
%     ay_b_INS_50Hz(i) = AccelTruth100Hz(2,i) + AN(2,i) + An2Bias;
%     az_b_INS_50Hz(i) = AccelTruth100Hz(3,i) + AN(3,i) + An3Bias;
% 
%     p_INS_50Hz(i) = GyroTruth100Hz(1,i) + GN(1,i) + Gn1Bias +  Gn1Fault*(i-7498) ;  
%     q_INS_50Hz(i) = GyroTruth100Hz(2,i) + GN(2,i) +  Gn2Bias +  Gn2Fault*(i-7498) ;  
%     r_INS_50Hz(i) = GyroTruth100Hz(3,i) + GN(3,i) + Gn3Bias +  Gn3Fault ; 
% 
% 
% 
% 
% %     
%  else
  
%     ax_b_INS_50Hz(i) = AccelTruth100Hz(1,i);% + AN(1,i) + An1Bias ; 
%     ay_b_INS_50Hz(i) = AccelTruth100Hz(2,i);% + AN(2,i) + An2Bias;
%     az_b_INS_50Hz(i) = AccelTruth100Hz(3,i);% + AN(3,i) + An3Bias;
% 
%     p_INS_50Hz(i) = GyroTruth100Hz(1,i);% + GN(1,i) + Gn1Bias;
%     q_INS_50Hz(i) = GyroTruth100Hz(2,i);% + GN(2,i) + Gn2Bias ;
%     r_INS_50Hz(i) = GyroTruth100Hz(3,i);% + GN(3,i) +Gn3Bias;
    
 %end

  
end

 



for i = startepoch-1:endepoch
    
     
    
   
   ax_b_INS1Hz(i) = ax_b_INS_50Hz(i*100+1);
    ay_b_INS1Hz(i) = ay_b_INS_50Hz(i*100+1);
    az_b_INS1Hz(i) = az_b_INS_50Hz(i*100+1);    

    p_INS1Hz(i) = p_INS_50Hz(i*100+1);
    q_INS1Hz(i) = q_INS_50Hz(i*100+1);
    r_INS1Hz(i) = r_INS_50Hz(i*100+1);       
       
    
    AN1Hz(1,i) =  AN(1,i*100+1) ;
      AN1Hz(2,i) =  AN(2,i*100+1) ;
        AN1Hz(3,i) =  AN(3,i*100+1) ;
        
          GN1Hz(1,i) =  GN(1,i*100+1) ;
            GN1Hz(2,i) =  GN(2,i*100+1) ;
              GN1Hz(3,i) =  GN(3,i*100+1) ;
    
    

end




% 
% 
% 
% Signal1 = GN(3,:)*10000;
% [CurveFit1, BestFit] = LeastSquaresBestFit(Signal1, 10, 1);
% 
% 
% ma = Signal1-CurveFit1;
% 
% 
% 
% [CurveFit2, BestFit] = LeastSquaresBestFit(ma, 10, 1);
% 
% 
% ma2 = ma-CurveFit2;
% 
% 
% 
% 
%    
% TruthINSNoise(1,:) = Signal1-CurveFit1;





