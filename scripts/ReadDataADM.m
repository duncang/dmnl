feature accel on

%read in truth data

%START THE READING IN AT 1 , EVEN THOUGH 101 IS THE START POINT for the simulation, (in 50 Hz intervals), BECAUSE
%IT USES i-1 IN SOME PLACES AND THIS MAKES SURE THESE VALUES AREN'T ZERO

percent = 0.01 ;%start at 1 percent
for i = startepochHighRate:endepochHighRate
    
    
      %this is aerosim output from its airplane block. Note this is the same as
    %XposCalc.
    
    
    
    
   
        
        
      %velocity truth
% 
%     V_n_truth(i) = vel_truth(2,i);
%       V_e_truth(i) = vel_truth(3,i);
%         V_d_truth(i) = vel_truth(4,i);
              
        
%I was using the vel_truthECEF.mat which is calculated from NED velocities
%in Aerosim but when the lat and lon parameters were switched around in the
%Tned2Ecef.m function i forgot to change it in the simulink. So I won't use
%vel_truthECEF anymore since some datasets will be in error for
%vel_truthECEF. 
%         Xvel_truth(i) =   vel_truthECEF(2,i);
%         Yvel_truth(i) =   vel_truthECEF(3,i);
%         Zvel_truth(i) =   vel_truthECEF(4,i);

        
%         TMatrix = T_ECEF2NED( Lat_truth(i), Lon_truth(i));        
%         
%         
%         Veceftemp = TMatrix'*[V_n_truth(i), V_e_truth(i), V_d_truth(i)]';
%         
%         
%          Xvel_truth(i) =   Veceftemp(1);
%          Yvel_truth(i) =   Veceftemp(2);
%          Zvel_truth(i) =   Veceftemp(3);
        
        
        
        

    %sensor measuremements
    %Read in INS sensor data

    ax_b_INS_50Hz(i) = ax_b_MODEL100Hz(i);
    ay_b_INS_50Hz(i) = ay_b_MODEL100Hz(i);
    az_b_INS_50Hz(i) = az_b_MODEL100Hz(i);

    p_INS_50Hz(i) = omega_xMODEL100Hz(i);
    q_INS_50Hz(i) =  omega_yMODEL100Hz(i);
    r_INS_50Hz(i) =  omega_zMODEL100Hz(i);
    
    
    
    
    GyroTruth100Hz(1,i) = omega_xMODEL100Hz(i);
        GyroTruth100Hz(2,i) = omega_yMODEL100Hz(i);
            GyroTruth100Hz(3,i) = omega_zMODEL100Hz(i);
 
    AccelTruth100Hz(1,i) = ax_b_MODEL100Hz(i);
    AccelTruth100Hz(2,i) = ay_b_MODEL100Hz(i);
    AccelTruth100Hz(3,i) = az_b_MODEL100Hz(i);
    
%     
%     
%     if i == 100*50+1
%         ax_b_INS_50Hz(i) =  ax_b_INS_50Hz(i) + 50;
%     end
    
    
%     %add noise on the IMU
%     
%     if i >= 2  %this is 50 seconds in 1 HZ rate
%         
%         percent = percent + 0.02; %add two percent at a time
%         %start at 5 percent
%         
% %        ax_b_INS(i) =  ax_b_INS(i) +  ax_b_INS(i)*percent*randn(1); 
% %         ay_b_INS(i) =  ay_b_INS(i) +  ay_b_INS(i)*percent*randn(1); 
% %          az_b_INS(i) =  az_b_INS(i) +  az_b_INS(i)*percent*randn(1); 
%          
%           p_INS(i) =  p_INS(i) +  p_INS(i)*percent*abs(randn(1)); %put abs so it only add its each time 
% %            q_INS(i) =  q_INS(i) +  q_INS(i)*percent*randn(1); 
% %             r_INS(i) =  r_INS(i) +  r_INS(i)*percent*randn(1); 
% %     
%        
%        %if up to 30 percent, start again
%        if percent >= 0.50
%            percent = 0.02;
%     
%        end
%        
%     end

    
  
        

    %Attitude truth

%     Roll_truth(i) = att_truth(2,i);   % bank angle(rad)  PHI
%     Pitch_truth(i) = att_truth(3,i);    %pitch angle (rad) THETA
%     Yaw_truth(i) = att_truth(4,i);    %heading (rad) PSI


    %Angular Accelerations

%     p_truth(i) = AngRate(2,i);   % p
%     q_truth(i) = AngRate(3,i);    %q
%     r_truth(i) = AngRate(4,i);    %r
% 
% 
%     p_dot_truth(i) = AngAcc(2,i);   % p dot
%     q_dot_truth(i) = AngAcc(3,i);    %q dot
%     r_dot_truth(i) = AngAcc(4,i);    %r dot


    %Air data
    %CORRECT

%     Airspeed_truth(i) = AirData1(2,i);    %Airspeed (m/s)
%     Beta_truth(i) = AirData1(3,i);     %Sideslip (rad)
%     Alpha_truth(i) = AirData1(4,i);              %AOA(rad)

    % flight controls
    %CORRECT

    Flap_truth(i) = flight_controls(2,i);          % Flap (rad)
    Elevator_truth(i) = flight_controls(3,i);       %Elevator (rad)
    Aileron_truth(i) = flight_controls(4,i);        % Aileron (rad)
    Rudder_truth(i) = flight_controls(5,i);         % Rudder (rad)
    Throttle_truth(i) = flight_controls(6,i);          %Throttle (frac (0 to 1??))
    Mixture_truth(i) = flight_controls(7,i);       %Mixture (ratio)
    Ignition_truth(i) = flight_controls(8,i);      %Ignition (bool)


    %coefficients aerosim calculates
    %aerodynamic coeffs.

    CD_truth(i) = AeroCoeff(2,i);    %rad^-1
    CY_truth(i) = AeroCoeff(3,i);     %ditto
    CL_truth(i) = AeroCoeff(4,i);
    Cl_truth(i) = AeroCoeff(5,i);
    Cm_truth(i) = AeroCoeff(6,i);
    Cn_truth(i) = AeroCoeff(7,i);


    %prop coeffs.

    J_truth(i) = PropCoeff(2,i);   % J
    CT_truth(i) = PropCoeff(3,i);    %CT
    CP_truth(i) = PropCoeff(4,i);    %CP

    %Eng. coeffs.

    MAP_truthEng(i) = EngCoeff(2,i);   % MAP  Current manifold pressure for current throttle setting and altitude, in kPa
    m_dotAir_truthEng(i) = EngCoeff(3,i);    %m_dotAir current instantaneous mass air flow
    m_dotFuel_truthEng(i) = EngCoeff(4,i);    %m_dotFuel current instantaneous mass fuel flow
    BSFC_truthEng(i) = EngCoeff(5,i);    %BSFC  Brake specific fuel consumption
    P_truthEng(i) = EngCoeff(6,i);    %P  Current engine power


    OMEGA_truthEng(i) = EngParameters(2,i);

    Faero_truth(1,i) = Forces(2,i);   % all forces In body axes
    Faero_truth(2,i) = Forces(3,i);   % all forces In body axes
    Faero_truth(3,i) = Forces(4,i);   % all forces In body axes

    Fprop_truth(1,i) = Forces(5,i);    %
    Fprop_truth(2,i) = Forces(6,i);    %
    Fprop_truth(3,i) = Forces(7,i);    %


    Mass_truth(i) = Forces(8,i);    %mass used in total force calculation

    Maero_truth(1,i) = Moments(2,i);   % all moments are in body axes
    Maero_truth(2,i) = Moments(3,i);   % all moments are in body axes
    Maero_truth(3,i) = Moments(4,i);   % all moments are in body axes


    Mprop_truth(1,i) = Moments(5,i);    %
    Mprop_truth(2,i) = Moments(6,i);    %
    Mprop_truth(3,i) = Moments(7,i);    %

    Mcg_truth(1,i) = Moments(8,i);    %Moment about the c of g
    Mcg_truth(2,i) = Moments(9,i);    %Moment about the c of g
    Mcg_truth(3,i) = Moments(10,i);    %Moment about the c of g


    Jx_truth(i) = InertiaCG(2,i);
    Jy_truth(i)= InertiaCG(3,i);
    Jz_truth(i)= InertiaCG(4,i);
    Jxz_truth(i) = InertiaCG(5,i);
    CGxpos_truth(i) = InertiaCG(6,i);
    CGypos_truth(i) = InertiaCG(7,i);
    CGzpos_truth(i) = InertiaCG(8,i);

    AtmosTruth(1,i)  = AtmosGrav(2,i);   %this is P
    AtmosTruth(2,i)  = AtmosGrav(3,i);                          %this is T
    AtmosTruth(3,i)  = AtmosGrav(4,i);                 %this is rho
    AtmosTruth(4,i)  = AtmosGrav(5,i);              %this is a (speed of sound)
    GravityTruth(i) =    AtmosGrav(6,i);

    Mach_truth(i) = Mach(2,i);

    %Quaternions

%     Quaternions_truth(1,i) = Quaternions(2,i); %q0
%     Quaternions_truth(2,i) = Quaternions(3,i); %q1
%     Quaternions_truth(3,i) = Quaternions(4,i); %q2
%     Quaternions_truth(4,i) = Quaternions(5,i); %q3
% 
%     
    
    
    

%This is wind speed in body axes, not u v w body velocities

    %Load wind speeds
    u_wind_truth(i) = WindB(2,i);
    v_wind_truth(i) = WindB(3,i);
    w_wind_truth(i) = WindB(4,i);
    
    
    %this is wind component on p q r

    p_wind_truth(i) = WindRates(2,i);
    q_wind_truth(i) = WindRates(3,i);
    r_wind_truth(i) = WindRates(4,i);
    
    


%     if i == 1
%         Xvel_truth(i) = 0;
%         Yvel_truth(i) = 0;
%         Zvel_truth(i) = 0;
%     else
% 
%         Xvel_truth(i) = (Xpos_truth(i) - Xpos_truth(i-1))/0.02;
%         Yvel_truth(i) = (Ypos_truth(i) - Ypos_truth(i-1))/0.02;
%         Zvel_truth(i) = (Zpos_truth(i) - Zpos_truth(i-1))/0.02;
%     end
% 
% 
%     TMatrix_ECEF2NED = T_ECEF2NED( Lat_truth(i),Lon_truth(i));
%     VelocityNED = TMatrix_ECEF2NED*[Xvel_truth(i),Yvel_truth(i),Zvel_truth(i)]';
%     
%     V_n_truth(i) = VelocityNED(1);
%       V_e_truth(i) = VelocityNED(2);
%         V_d_truth(i) = VelocityNED(3);



% if UseNoisy == 1
%     load('data\TestFlight4.7.07\sensors_noisy.mat');
%     sensors = sensors_noisy;
%     clear 'sensors_noisy';
% else
%     load('data\TestFlight4.7.07\sensors_clean.mat');
%     sensors = sensors_clean;
%     clear 'sensors_clean';
% end


% 
% 
% INSSensorNoise(1,i) = sensorsNoisy(2,i) - sensorsClean(2,i);            %accel.
% INSSensorNoise(2,i) = sensorsNoisy(3,i) - sensorsClean(3,i);            %accel.
% INSSensorNoise(3,i) = sensorsNoisy(4,i) - sensorsClean(4,i);            %accel.
% INSSensorNoise(4,i) = sensorsNoisy(5,i) - sensorsClean(5,i);            %gyro.
% INSSensorNoise(5,i) = sensorsNoisy(6,i) - sensorsClean(6,i);            %gyro.
% INSSensorNoise(6,i) = sensorsNoisy(7,i) - sensorsClean(7,i);            %gyro.
% 

%  AccelTruth100Hz(1:3,i) = sensorsClean(2:4,i);
%  GyroTruth100Hz(1:3,i) = sensorsClean(5:7,i);


%AccelTruth100Hz(1:3,i) = sensorsNoisy(2:4,i);
%GyroTruth100Hz(1:3,i) = sensorsNoisy(5:7,i);


end;



%truth data at 1 Hz


%for i = startepoch - 1


i = startepoch - 1; 









for i = startepoch:endepoch

    %read in the data and perform necessary conversions

    %this should be the GPS inputs but use the truth for now.

    %ECEF positions

    %Note that first row is simulation time data thats why the index starts at
    %2 in the following

    Lat_truth1Hz(i) = x_state_total_plusMODEL(8,i);
    Lon_truth1Hz(i) = x_state_total_plusMODEL(9,i);
    Hgt_truth1Hz(i) = x_state_total_plusMODEL(10,i);
    
    
        
    
    
    
     Lat_truth(i) = pos_truth_llh(2,i);
    Lon_truth(i) = pos_truth_llh(3,i);
    Hgt_truth(i) = pos_truth_llh(4,i);
    
    
    
    xtemp = LLH2ECEF(Lat_truth1Hz(i), Lon_truth1Hz(i),   Hgt_truth1Hz(i));
    
    
    Xpos_truth1Hz(i) = xtemp(1);
    Ypos_truth1Hz(i) = xtemp(2);
    Zpos_truth1Hz(i) = xtemp(3);     
    


    
      V_n_truth1Hz(i) = x_state_total_plusMODEL(5,i);
      V_e_truth1Hz(i) = x_state_total_plusMODEL(6,i);
      V_d_truth1Hz(i) = x_state_total_plusMODEL(7,i);
      
      
      
      TMatrix = T_ECEF2NED(Lat_truth1Hz(i), Lon_truth1Hz(i));
      
    
      
          
    
    xveltemp = TMatrix'*[ V_n_truth1Hz(i), V_e_truth1Hz(i), V_d_truth1Hz(i)]';
    

    Xvel_truth1Hz(i) = xveltemp(1);
    Yvel_truth1Hz(i) = xveltemp(2);
    Zvel_truth1Hz(i) = xveltemp(3);
    
      
      
    
    %note this is using uncorrected roll, pitch and yaw estimates from theINS
    %these are also corrupted by the wind I think?

    % Xvel_GPS(i) = Xvel_truth1Hz(i);
    % Yvel_GPS(i) = Yvel_truth1Hz(i);
    % Zvel_GPS(i) = Zvel_truth1Hz(i);


     Quaternions_truth1Hz(1,i) = x_state_total_plusMODEL(1,i);
        Quaternions_truth1Hz(2,i) = x_state_total_plusMODEL(2,i);
            Quaternions_truth1Hz(3,i) = x_state_total_plusMODEL(3,i);
                Quaternions_truth1Hz(4,i) = x_state_total_plusMODEL(4,i);
                
                
                
           
                  
    quat = [   Quaternions_truth1Hz(1,i) ,  Quaternions_truth1Hz(2,i) ,  Quaternions_truth1Hz(3,i) ,  Quaternions_truth1Hz(4,i) ];     
      
    [euler] = QuatToEuler(quat);
    
    
          
    Roll_truth1Hz(i) = euler(1);
    Pitch_truth1Hz(i) = euler(2);
    Yaw_truth1Hz(i) = euler(3);
       
    
    
    
    %calculate u v w truth
    
    TMatrix = T_Body2NED(Roll_truth1Hz(i),Pitch_truth1Hz(i), Yaw_truth1Hz(i));
    
    uvwtruth1hz = TMatrix'*[V_n_truth1Hz(i), V_e_truth1Hz(i), V_d_truth1Hz(i)]';
    u_truth1Hz(i) = uvwtruth1hz(1);
     v_truth1Hz(i) = uvwtruth1hz(2);
      w_truth1Hz(i) = uvwtruth1hz(3);
    
    
      
      
      
          
    Airspeed_truth1Hz(i) =  sqrt( u_truth1Hz(i)^2 +  v_truth1Hz(i)^2 +  w_truth1Hz(i)^2);
    Beta_truth1Hz(i) =  asin(v_truth1Hz(i)/Airspeed_truth1Hz(i) );
    Alpha_truth1Hz(i) =   atan2(w_truth1Hz(i),u_truth1Hz(i));
    
      
      
      
    
    
    
%     Roll_truth1Hz(i) = mean(Roll_truth((i-1)*50+1:i*50+1));   
%     Pitch_truth1Hz(i) = mean(Pitch_truth((i-1)*50+1:i*50+1));    
%     Yaw_truth1Hz(i) = mean(Yaw_truth((i-1)*50+1:i*50+1));    
%     
       

                
               %convert from quat to euler 
                
                
GyroTruth1Hz(1,i) = GyroTruth100Hz(1,i*100+1) ;       
  GyroTruth1Hz(2,i) = GyroTruth100Hz(2,i*100+1)   ;               
     GyroTruth1Hz(3,i) = GyroTruth100Hz(3,i*100+1)  ; 
     
     
  AccelTruth1Hz(1,i) =  AccelTruth100Hz(1,i*100+1);
AccelTruth1Hz(2,i) =  AccelTruth100Hz(2,i*100+1);
AccelTruth1Hz(3,i) =  AccelTruth100Hz(3,i*100+1);

     
             
                
  %sensor data at 1 Hz
   
   ax_b_INS1Hz(i) = ax_b_INS_50Hz(i*100+1);
    ay_b_INS1Hz(i) = ay_b_INS_50Hz(i*100+1);
    az_b_INS1Hz(i) = az_b_INS_50Hz(i*100+1);
    

    p_INS1Hz(i) = p_INS_50Hz(i*100+1);
    q_INS1Hz(i) = q_INS_50Hz(i*100+1);
       r_INS1Hz(i) = r_INS_50Hz(i*100+1);       



% 
% 
%     ax_b_INS_50Hz(i) = ax_b_MODEL100Hz(i);
%     ay_b_INS_50Hz(i) = ay_b_MODEL100Hz(i);
%     az_b_INS_50Hz(i) = az_b_MODEL100Hz(i);
% 
%     p_INS_50Hz(i) = omega_xMODEL100Hz(i);
%     q_INS_50Hz(i) =  omega_yMODEL100Hz(i);
%     r_INS_50Hz(i) =  omega_zMODEL100Hz(i);

%        
       
       
    CD_truth1Hz(i) = CD_truth(i*100+1);    %rad^-1
    CY_truth1Hz(i) = CY_truth(i*100+1);     %ditto
    CL_truth1Hz(i) = CL_truth(i*100+1);
    Cl_truth1Hz(i) = Cl_truth(i*100+1);
    Cm_truth1Hz(i) = Cm_truth(i*100+1);
    Cn_truth1Hz(i) = Cn_truth(i*100+1);
                
                
%                 
%     
%     
%     
%     INSSensorNoise1Hz(1,i) = INSSensorNoise(1,i*100+1);
%         INSSensorNoise1Hz(2,i) = INSSensorNoise(2,i*100+1);
%             INSSensorNoise1Hz(3,i) = INSSensorNoise(3,i*100+1);
%                 INSSensorNoise1Hz(4,i) = INSSensorNoise(4,i*100+1);
%                     INSSensorNoise1Hz(5,i) = INSSensorNoise(5,i*100+1);
%                         INSSensorNoise1Hz(6,i) = INSSensorNoise(6,i*100+1);
%     
    
    
                
   %Roll_truth1HzTest(i) = mean(Roll_truth((i-1)*50+1:i*50+1));   
   
   %sensor data at 1 Hz
%    
%    ax_b_INS1Hz(i) = mean(ax_b_INS((i-1)*50+1:i*50+1));
%     ay_b_INS1Hz(i) = mean(ay_b_INS((i-1)*50+1:i*50+1));
%     az_b_INS1Hz(i) = mean(az_b_INS((i-1)*50+1:i*50+1));
% 
%     p_INS1Hz(i) = mean(p_INS((i-1)*50+1:i*50+1));
%     q_INS1Hz(i) = mean(q_INS((i-1)*50+1:i*50+1));
%        r_INS1Hz(i) = mean(r_INS((i-1)*50+1:i*50+1));






% 
% 
% AccelTruth1Hz(1:3,i) = sensorsClean(2:4,i*100+1);
% GyroTruth1Hz(1:3,i) = sensorsClean(5:7,i*100+1);


%This is wind speed in body axes, not u v w body velocities

%     u_wind_truth1Hz(i) = u_wind_truth(i*100+1);
%     v_wind_truth1Hz(i) = v_wind_truth(i*100+1);
%     w_wind_truth1Hz(i) = w_wind_truth(i*100+1);
% 
%     p_wind_truth1Hz(i) = p_wind_truth(i*100+1);
%     q_wind_truth1Hz(i) = q_wind_truth(i*100+1);
%     r_wind_truth1Hz(i) = r_wind_truth(i*100+1);
%     
                

end





i = startepoch -1;


    Lat_truth1Hz(i) = x_state_total_plusMODEL(8,i+1);
    Lon_truth1Hz(i) = x_state_total_plusMODEL(9,i+1);
    Hgt_truth1Hz(i) = x_state_total_plusMODEL(10,i+1);
    
    
        
    
    
    
     Lat_truth(i) = pos_truth_llh(2,i+1);
    Lon_truth(i) = pos_truth_llh(3,i+1);
    Hgt_truth(i) = pos_truth_llh(4,i+1);
    
    
    
    xtemp = LLH2ECEF(Lat_truth1Hz(i+1), Lon_truth1Hz(i+1),   Hgt_truth1Hz(i+1));
    
    
    Xpos_truth1Hz(i) = xtemp(1);
    Ypos_truth1Hz(i) = xtemp(2);
    Zpos_truth1Hz(i) = xtemp(3);     
    


    
      V_n_truth1Hz(i) = x_state_total_plusMODEL(5,i+1);
      V_e_truth1Hz(i) = x_state_total_plusMODEL(6,i+1);
      V_d_truth1Hz(i) = x_state_total_plusMODEL(7,i+1);
      
      
      
      TMatrix = T_ECEF2NED(Lat_truth1Hz(i+1), Lon_truth1Hz(i+1));
      
    
      
          
    
    xveltemp = TMatrix'*[ V_n_truth1Hz(i+1), V_e_truth1Hz(i+1), V_d_truth1Hz(i+1)]';
    

    Xvel_truth1Hz(i) = xveltemp(1);
    Yvel_truth1Hz(i) = xveltemp(2);
    Zvel_truth1Hz(i) = xveltemp(3);
    
      
      
    
    %note this is using uncorrected roll, pitch and yaw estimates from theINS
    %these are also corrupted by the wind I think?

    % Xvel_GPS(i) = Xvel_truth1Hz(i);
    % Yvel_GPS(i) = Yvel_truth1Hz(i);
    % Zvel_GPS(i) = Zvel_truth1Hz(i);


     Quaternions_truth1Hz(1,i) = x_state_total_plusMODEL(1,i+1);
        Quaternions_truth1Hz(2,i) = x_state_total_plusMODEL(2,i+1);
            Quaternions_truth1Hz(3,i) = x_state_total_plusMODEL(3,i+1);
                Quaternions_truth1Hz(4,i) = x_state_total_plusMODEL(4,i+1);
                
                
                
           
                  
    quat = [   Quaternions_truth1Hz(1,i+1) ,  Quaternions_truth1Hz(2,i+1) ,  Quaternions_truth1Hz(3,i+1) ,  Quaternions_truth1Hz(4,i+1) ];     
      
    [euler] = QuatToEuler(quat);
    
    
          
    Roll_truth1Hz(i) = euler(1);
    Pitch_truth1Hz(i) = euler(2);
    Yaw_truth1Hz(i) = euler(3);
       
    
    
    
    %calculate u v w truth
    
    TMatrix = T_Body2NED(Roll_truth1Hz(i+1),Pitch_truth1Hz(i+1), Yaw_truth1Hz(i+1));
    
    uvwtruth1hz = TMatrix'*[V_n_truth1Hz(i+1), V_e_truth1Hz(i+1), V_d_truth1Hz(i+1)]';
    u_truth1Hz(i) = uvwtruth1hz(1);
     v_truth1Hz(i) = uvwtruth1hz(2);
      w_truth1Hz(i) = uvwtruth1hz(3);
    
    
      
      
      
          
    Airspeed_truth1Hz(i) =  sqrt( u_truth1Hz(i+1)^2 +  v_truth1Hz(i+1)^2 +  w_truth1Hz(i+1)^2);
    Beta_truth1Hz(i) =  asin(v_truth1Hz(i+1)/Airspeed_truth1Hz(i+1) );
    Alpha_truth1Hz(i) =   atan2(w_truth1Hz(i+1),u_truth1Hz(i+1));
    
      
      
      
    
    
    
%     Roll_truth1Hz(i) = mean(Roll_truth((i-1)*50+1:i*50+1));   
%     Pitch_truth1Hz(i) = mean(Pitch_truth((i-1)*50+1:i*50+1));    
%     Yaw_truth1Hz(i) = mean(Yaw_truth((i-1)*50+1:i*50+1));    
%     
    
    

                
               %convert from quat to euler 
                
            
             
                
  %sensor data at 1 Hz
   
   ax_b_INS1Hz(i) =  ax_b_INS1Hz(i+1);
    ay_b_INS1Hz(i) =  ay_b_INS1Hz(i+1);
    az_b_INS1Hz(i) =  az_b_INS1Hz(i+1) ;
    

    p_INS1Hz(i) = p_INS1Hz(i+1);
    q_INS1Hz(i) = q_INS1Hz(i+1);
     r_INS1Hz(i) = r_INS1Hz(i+1);  


GyroTruth1Hz(1,i) = GyroTruth1Hz(1,i+1);
GyroTruth1Hz(2,i) = GyroTruth1Hz(2,i+1);
GyroTruth1Hz(3,i) = GyroTruth1Hz(3,i+1);

AccelTruth1Hz(1,i) =  AccelTruth1Hz(1,i+1);
AccelTruth1Hz(2,i) =  AccelTruth1Hz(2,i+1);
AccelTruth1Hz(3,i) =  AccelTruth1Hz(3,i+1);
