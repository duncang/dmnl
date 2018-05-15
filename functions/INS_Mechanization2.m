function [X_state_out] = INS_Mechanization2(X_state_in, A_b, omega_b, gravity, ins_dt,llh_dot);



OMEGAedot = 7.2921151467e-5;

% $Id: INS_Mechanization2.m 1887 2008-07-15 06:11:51Z n2523710 $

%This function takes gyro and accel outputs from inertial sensors, which are assumed to be already corrected 


%note that gravity should be positive ie +9.8 m/s. not negative.


% x_state_in vector is as follows:
% 
% q0
% q1
% q2
% q3
% Vn
% Ve
% Vd
% Lat
% Lon
% Hgt


lat_dot = llh_dot(1);
lon_dot = llh_dot(2);
hgt_dot = llh_dot(3);

%OMEGA_e = 7.292115e-5; % eearth rotation rate (inertial frame)
    
    
    %Attitude
    
    q0 = X_state_in(1);
    q1 = X_state_in(2);
    q2 = X_state_in(3);
    q3 = X_state_in(4);       
    
    
                 %normalise quaternions 
[q0,q1,q2,q3] = Normalise_Quat(q0,q1,q2,q3);

                
    %convert to phi, theta psi
    
%    phi = atan2( 2*(q0*q1 + q2*q3), q0^2 + q3^2 - q1^2 -q2^2);
%    theta = asin(2*(q0*q2 - q1*q3));
%    psi = atan2( 2*(q0*q3 + q1*q2), q0^2 + q1^2 - q2^2 - q3^2);
%     
    
%     InitialAttitude(1) = phi;
%     InitialAttitude(2) = theta;
%     InitialAttitude(3) = psi;
    
     %velocity
        Vn = X_state_in(5);  %north 
           Ve = X_state_in(6);  %east 
               Vd = X_state_in(7); %down
    
    %position
        lat = X_state_in(8);  %lat
            lon = X_state_in(9); %lon
                hgt = X_state_in(10); %hgt           
                          
                
         
    %INS measurements
        omega_x = omega_b(1);
        omega_y = omega_b(2);
        omega_z = omega_b(3);

        A_x = A_b(1);
        A_y = A_b(2);
        A_z = A_b(3);
        
        g = gravity;  % gravity is local gravity as computed by aerosim. Note that g should be positive i.e. +9.8 because A_n is accel in down direction, so
        %gravity must add to this error.
                        
    
        %propagation of attitude
    p = omega_x;
      q = omega_y;
     r = omega_z;
     
     
          
     
%      Q_dot = pqrToQuat([p, q, r], [q0,q1,q2,q3]);
%               
%        %Q_dot = -0.5*[0, p, q, r; -p, 0, -r, q; -q, r, 0, -p; -r, -q, p, 0;]*[q0;q1;q2;q3];     
%     
% 
% %         a_dot(i) = -0.5 * (bq * omega_x + cq * omega_y + dq * omega_z);
% %         b_dot(i) = 0.5 * (aq * omega_x - dq * omega_y + cq * omega_z);
% %         c_dot(i) = 0.5 * (dq * omega_x + aq * omega_y - bq * omega_z);
% %         d_dot(i) = -0.5 * (cq * omega_x - bq * omega_y - aq * omega_z);
% 
%         q0_dot = Q_dot(1);
%          q1_dot = Q_dot(2);
%           q2_dot = Q_dot(3);
%            q3_dot = Q_dot(4);      
%            
%            
% %               %Renormalise quaternions
% %         vector_length = sqrt(q0^2 + q1^2 + q2^2 + q3^2);
% %         
% %         q0 = q0/vector_length;
% %          q1 = q1/vector_length;
% %           q2 = q2/vector_length;
% %             q3 = q3/vector_length;
%          
% 
%          
%            
%         q0_new = q0 + q0_dot * ins_dt;
%         q1_new = q1 + q1_dot * ins_dt;
%         q2_new = q2 + q2_dot * ins_dt;
%         q3_new = q3 + q3_dot * ins_dt;
        
        
        
        
        %quaternion update using formula from 'sigma-point kalman filters
        %for integrated navigation by der Merwe and Wan
        
        
        OMEGA_tilda = [0, p, q, r; -p, 0, -r, q; -q, r, 0, -p; -r, -q, p, 0;];
        
        q_new = expm(-0.5*OMEGA_tilda*ins_dt)*[q0,q1,q2,q3]';
        
        
        q0_new = q_new(1);
        q1_new = q_new(2);
        q2_new = q_new(3);
        q3_new = q_new(4);
        
        
    
        
        %Renormalise quaternions
%         vector_length = sqrt(q0_new^2 + q1_new^2 + q2_new^2 + q3_new^2);
%         
%         q0_new = q0_new/vector_length;
%          q1_new = q1_new/vector_length;
%           q2_new = q2_new/vector_length;
%             q3_new = q3_new/vector_length;
                  
            
            %normalise quaternions 
            %no need to normalise them with the closed form solution i used
            %above
%[q0_new,q1_new,q2_new,q3_new] = Normalise_Quat(q0_new,q1_new,q2_new,q3_new);

            
        
        %note that this DCM is different from the one Aerosim calls the
        %DCM. This is CBN though (transform from body to nav frame). This
        %has been checked
         CBN = GARD_QuatToDCM([q0_new,q1_new,q2_new,q3_new]);
        
         Accel_nav = CBN*[A_x, A_y, A_z]';
        
        an = Accel_nav(1);
            ae = Accel_nav(2);
                ad = Accel_nav(3);

% 
%         an = A_x;
%             ae =  A_y;
%                 ad =  A_z;
        
        


 [Rn, Re] = WGS84_calcRnRe(lat);     
    Rnh = Rn + hgt;
    Reh = Re + hgt;

        
        
%my original implementation, from titterton

       %Vn_dot = an - 2*OMEGAedot*sin(lat)*Ve + ( Vn*Vd - Ve^2*tan(lat))/Reh;  %should this be Reh or Rnh, check 
       
       Vn_dot = an - 2*OMEGAedot*sin(lat)*Ve + ( Vn*Vd - Ve^2*tan(lat))/Rnh;
       Ve_dot = ae + 2*OMEGAedot*(Vn*sin(lat) + Vd*cos(lat)) + (Ve*Vd)/Reh + Ve*Vn*tan(lat)/Reh;
       Vd_dot = ad - 2*OMEGAedot*Ve*cos(lat) - (Ve^2/Reh) - (Vn^2/Rnh)   + g;  %note that g shoudl be +ve here, to add to the down acceleration, because gravity acts downward
        
                
        
        %no coriolis
        
    %    Vn_dot = an + (Vn*Vd - Ve^2*tan(lat))/Reh;
    %    Ve_dot = ae  +  (Ve*Vd)/Reh +  Ve*Vn*tan(lat)/Reh;
    %    Vd_dot = ad  - (Ve^2/Reh) - (Vn^2/Rnh)   + g;  %note that g shoudl be +ve here, to add to the down acceleration, because gravity acts downward
        


%checking to see whether changing this improves anything
%         Vn_dot = an - 2*OMEGAedot*sin(lat)*Ve                      + (Vn*Vd - Ve^2*tan(lat))/Reh;
%         Ve_dot = ae + 2*OMEGAedot*(Vn*sin(lat) + Vd*cos(lat))      + (Ve*Vd)/Reh + Ve*Vn*tan(lat)/Reh;
%         Vd_dot = ad +0                   - (Ve^2/Reh) - (Vn^2/Rnh)   + g;  %note that g shoudl be +ve here, to add to the down acceleration, because gravity acts downward
%         
%         

        
        %this is how aerosim has it, taken from Farrel and Barth (the
        %global positionig system and inertial navigation)
        
%         Vn_dot = -(lon_dot + 2*OMEGAedot)*sin(lat)*Ve + lat_dot*Vd + an;
%         Ve_dot = (lon_dot + 2*OMEGAedot)*sin(lat)*Vn + (lon_dot + 2*OMEGAedot)*cos(lat)*Vd + ae;
%         Vd_dot = -lat_dot*Vn - (lon_dot + 2*OMEGAedot)*cos(lat)*Ve + ad + g; 
        


% 
%         Vn_dot = an + (Vn*Vd - Ve^2*tan(lat))/Reh;
%         Ve_dot = ae  +  (Ve*Vd)/Reh +  Ve*Vn*tan(lat)/Reh;
%         Vd_dot = ad  - (Ve^2/Reh) - (Vn^2/Rnh);   %+ g;  %note that g shoudl be +ve here, to add to the down acceleration, because gravity acts downward
%         
%         
% %     
%           %Renormalise quaternions
%         vector_length = sqrt(q0_new^2 + q1_new^2 + q2_new^2 + q3_new^2);
%         
%         q0_new = q0_new/vector_length;
%          q1_new = q1_new/vector_length;
%           q2_new = q2_new/vector_length;
%             q3_new = q3_new/vector_length;
%          
        
        
   
        %integrate velocities to get new ones
   Vn_new = Vn + Vn_dot*ins_dt;
   Ve_new = Ve + Ve_dot*ins_dt;
   Vd_new = Vd + Vd_dot*ins_dt;
        
                 
  
Lat_dot = Vn_new/Rnh;
Lon_dot = Ve_new/(Reh*cos(lat));
Hgt_dot = -Vd_new;       
        
        
   % propogate position estimate
        
        
   lat_new = lat + Lat_dot*ins_dt;
   lon_new = lon + Lon_dot*ins_dt;
   hgt_new = hgt + Hgt_dot*ins_dt;        
      



X_state_out(1) = q0_new;
X_state_out(2) = q1_new;
X_state_out(3) = q2_new;
X_state_out(4) = q3_new;
X_state_out(5) = Vn_new; %vn
X_state_out(6) = Ve_new; %ve
X_state_out(7) = Vd_new; %vd
X_state_out(8) = lat_new; %lat
X_state_out(9) = lon_new; %lon
X_state_out(10) = hgt_new; %hgt



