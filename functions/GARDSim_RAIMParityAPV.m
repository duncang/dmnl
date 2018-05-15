function [BadGeometry_H,BadGeometry_V, RAIM_ALERT_H,RAIM_ALERT_V, SLOPE_Max_H,SLOPE_Max_V, r_H,r_V, Td_H,Td_V, HPL,VPL, FaultySatFDI] = GARDSim_RAIMParityAPV(a_H,a_V, lambda_H,lambda_V, N,PFalseAlarm,SigmaS,Alert_Limit_HAL,Alert_Limit_VAL,ResVec,M);

  

%$Id: GARDSim_RAIMParity.m 841 2007-06-13 01:52:53Z greerd $
%Troy Bruggemann 19 January 2006
%
%Does RAIM FDI using Parity Scheme - this 
% INPUT
%a is the vector of pre-computed values to be used to calculate the threshold,
%with degrees of freedom range from 2 to 18 (ie number of satellites from 6
%to 24
%
%lambda is the vector of pre-computed noncentrality parameters from chi square distibution to be used to calculate the pbias,
%with degrees of freedom range from 2 to 18 (ie number of satellites from 6
%to 24
%
%
% N - Number of Satellites to do RAIM solution with
% PFalseAlarm - Probability of False Alarm
% SigmaS - Expected Standard Deviation of Pseudorange errors
% Alert_Limit - 3D Protection limit for given phase of flight
% ResVec - [1..N] Vector of Pseudorange residuals from least squares (m) 
% M = M matrix from least squares
%=======================================================================
% OUTPUT
% BadGeometry - Flag indicating poor satellite-user geometry   
% RAIM_ALERT - Flag indicating poor position solution integrity
% SLOPE_Max - Value of Maximum Slope
% r - test statistic


%for testing and verification
%SigmaS = 33;  %metres
%PFalseAlarm = 1e-6;
%Pmiss = 0.001

%5 in view case in RTCA paper by Brown
% M = [-0.8607445743 -0.3446039300 0.3746557209 1; 0.2109370676 0.3502943374 0.9125784518 1; -0.0619331310 -0.4967359072 0.8656891623 1;...
%     -0.7248969588 0.4759681238 0.4979746422 1; -0.4009266538 0.1274180997 0.9072058455 1]

%6 in view case in RTCA paper by Brown
% M = [0.7460266527 -0.4689257437 0.4728137904 1;...
%     -0.8607445743 -0.3446039300 0.3746557209 1; 0.2109370676 0.3502943374 0.9125784518 1; -0.0619331310 -0.4967359072 0.8656891623 1;...
%     -0.7248969588 0.4759681238 0.4979746422 1; -0.4009266538 0.1274180997 0.9072058455 1]
% 



%in the function header pass the a values and the pbias values from offline
%calculations.



%check that Number of satellites enough to do RAIM FDI

 if N < 5;
       %error('Not enough satellites for RAIM solution')  
       

    
    
    
    
    
    
    BadGeometry_H = 100;
    BadGeometry_V = 100;
    RAIM_ALERT_H = 0;
    RAIM_ALERT_V = 0;
    SLOPE_Max_H = 1;
    SLOPE_Max_V = 1;
    r_H= 0;
    r_V= 0;
    Td_H= 0;
    Td_V= 0;
    HPL= 0;
    VPL= 0;
    FaultySatFDI= 0;
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     return;  %exit function.
       
 end
  
 
 
 
 %first values in thhe a and lambda 
 
 
 %Index to the a and lambda array
 
 a_ind = N-4;
 
 %For Horizontal Direction
 %-------------------------------
 %determine threshold value
 %------------------------------ 
%normalised Td
Td_norm_H = sqrt(a_H(a_ind));
%Unnormalised Td
Td_H = SigmaS*Td_norm_H; %metres
  
 %-------------------------
 %calculate pbias values
 %-------------------------
 %normalised pbias 
  
 pbias_norm_H = sqrt(lambda_H(a_ind));
 pbias_H = pbias_norm_H*SigmaS;
 



%determine SLOPE values and SLOPE_Max

%IMPORTANT: M MUST BE CONVERTED TO LTP FIRST BEFORE USE HERE, M STRAIGHT
%OUT OF GARD_LSQ WILL BE IN ECEF, NOT LTP.
AA = inv(M'*M)*M';
BB = M*AA;

%compute Slope for each satellite in view , X, Y and Z (~GDOP)

% for i_slope = 1:N
%     SLOPE(i_slope) = sqrt(  (AA(1,i_slope)^2 + AA(2,i_slope)^2 + AA(3,i_slope)^2)*(N-4)/(1 - BB(i_slope,i_slope)));
% end

%Slope for X and Y only
%  for i_slope = 1:N
%      SLOPE(i_slope) = sqrt(  (AA(1,i_slope)^2 + AA(2,i_slope)^2)*(N-4)/(1 - BB(i_slope,i_slope)));
%    
%  end
 

%See "A measurement domain receiver autonomous integrity monitoring
%algorithm" by Feng, Ochieng, for fairly good explanation on difference
%between N-4 and N. 

%Since Brown's paper (GPS RAIM:Calculation of Thresholds and protection radius etc") doesn't use N-4, I won't use it either. 
 
 for i_slope = 1:N
     SLOPE(i_slope) = sqrt(  (AA(1,i_slope)^2 + AA(2,i_slope)^2))/(1 - BB(i_slope,i_slope));
   
 end
 
 
     
SLOPE_Max_H = max(SLOPE);  %the value of slope_max
maxslopsatellite_H = find(SLOPE == SLOPE_Max_H) ;%i_slope indice to the max satellite ie. which satellite has the max_slope



%Form Parity transformation Matrix
[Q,R] = qr(M);  %M is also called the G matrix.

%get submatrix %sub-matrix with lower n-4 rows of Q' is P matrix:

Qt = Q';
P = Qt(N-(N-4)+1:N,1:N);



%calculate range bias

%The maxslopsatellite'th column of P is what we want:

rbias_H = pbias_H/norm(P(:,maxslopsatellite_H));

%compute the position error vector induced by the bias of rbias in the max
%slope satellite
epsilon = zeros(N,1);  %this is already a column vector, dont need to transpose it
epsilon(maxslopsatellite_H) = rbias_H;
Poserr = AA*epsilon;

HPL = norm([Poserr(1) Poserr(2)]);






 %For Vertical Direction
 %-------------------------------
 %determine threshold value
 %------------------------------ 
%normalised Td
Td_norm_V = sqrt(a_V(a_ind));
%Unnormalised Td
Td_V = SigmaS*Td_norm_V; %metres
  
 %-------------------------
 %calculate pbias values
 %-------------------------
 %normalised pbias 
  
 pbias_norm_V = sqrt(lambda_V(a_ind));
 pbias_V = pbias_norm_V*SigmaS;
 



%% calculate VPL

 %% calculate V-slope
 for i_slope = 1:N
    SLOPE_V(i_slope) = sqrt(  (AA(3,i_slope)^2 ))/(1 - BB(i_slope,i_slope));
 end



SLOPE_Max_V = max(SLOPE_V);  %the value of slope_max
maxslopsatellite_V = find(SLOPE_V == SLOPE_Max_V) ;%i_slope indice to the max satellite ie. which satellite has the max_slope


%note, I think this is norm(P(:,maxslopsatellite_H) not
%norm(P(3,maxslopsatellite_H) because P is not a 3*N matrix but an N-4*N
%matrix, so the rows are not north, east, up. 
rbias_V = pbias_V/norm(P(:,maxslopsatellite_H));  %NOTE, I USE MAXSLOPSATELILTE_H here, as it shouldnt be different for the vertical case

%compute the position error vector induced by the bias of rbias in the max
%slope satellite
epsilon = zeros(N,1);  %this is already a column vector, dont need to transpose it
epsilon(maxslopsatellite_H) = rbias_V;
Poserr = AA*epsilon;

VPL = norm(Poserr(3));


%Check for RAIM not being available:

if HPL > Alert_Limit_HAL
    
BadGeometry_H = 100;
else
    BadGeometry_H = 0;    
end


if VPL > Alert_Limit_VAL
    
BadGeometry_V = 100;
else
    BadGeometry_V = 0;    
end




%Detect any errors


%Calculate parity vector
%horizontal direction

p = P*ResVec; 




%Use magnitude of p as test statistic. Alternatively could use 
%SSE = p'*p i think.

if N>5
r_H = norm([p(1),p(2)]);

else  %special case if N = 5. 
r_H = p;      
end

%SSE = abs(p)

%SSE = norm(p)

%SSE = p;

%compare to threshold 

if r_H > Td_H 
    RAIM_ALERT_H = 100;
else
    RAIM_ALERT_H = 0;
    
end



% 
% 
% if N>5
% r_V = norm([p(3)]);  23.5.08 NEED TO CHECK THIS I DONT THINK ITS RIGHT
% BECAUSE this is N-4, not North east down
% 
% else  %special case if N = 5. 
% r_V = p;      
% end

r_V = 0;  

if r_V > Td_V 
    RAIM_ALERT_V = 100;
else
    RAIM_ALERT_V = 0;
    
end















% 
% 
% 
% %Attempt FDI using parity method - only detects one faulty measurement
% %If a fault is detected
% %Run only if there are 2 or more redundant measurement sources available
% if RAIM_ALERT == 100
%         if N>=6
%             
%          %compute fi^2/Sii quantities
%          
%          for i = 1:N
%              
%              FDItest(i) = (f(i)^2)/S(i,i);
%             
%          end
%          
%          %find maximum
%          FDItestMax = max(FDItest);
%          
%          FaultySatFDI = find(FDItest == FDItestMax);       
%             
%     
%         end
%     else
%         
%         FaultySatFDI = 0;
%         
%     end
%     
% 
% 
% 
FaultySatFDI = 1;


%just for output, output the index to the satellite wiht max slope
SLOPE_Max_H = maxslopsatellite_H;





