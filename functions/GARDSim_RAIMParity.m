function [BadGeometry, RAIM_ALERT, SLOPE_Max, r, Td, HPL,VPL, FaultySatFDI] = GARDSim_RAIMParity(a, lambda, N,PFalseAlarm,SigmaS,Alert_Limit,ResVec,M)
%$Id: GARDSim_RAIMParity.m 3023 2009-10-18 00:37:28Z greerd $
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
       
    BadGeometry = 100;       
    RAIM_ALERT = 0;
    SLOPE_Max = 0;
    r = 0;
    Td = 0;
    HPL =0;
    VPL =0;
    FaultySatFDI=0;
    
    
     return;  %exit function.
       
 end
  
 
 
 
 %first values in thhe a and lambda 
 
 
 %Index to the a and lambda array
 
 a_ind = N-4;
 
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
 



%determine SLOPE values and SLOPE_Max
AA = inv(M'*M)*M';
BB = M*AA;

%compute Slope for each satellite in view , X, Y and Z (~GDOP)

% for i_slope = 1:N
%     SLOPE(i_slope) = sqrt(  (AA(1,i_slope)^2 + AA(2,i_slope)^2 + AA(3,i_slope)^2)*(N-4)/(1 - BB(i_slope,i_slope)));
% end

%Slope for X and Y only
 for i_slope = 1:N
     SLOPE(i_slope) = sqrt(  (AA(1,i_slope)^2 + AA(2,i_slope)^2)*(N-4)/(1 - BB(i_slope,i_slope)));
   
 end
 
 
     
SLOPE_Max = max(SLOPE);  %the value of slope_max
maxslopsatellite = find(SLOPE == SLOPE_Max); %i_slope indice to the max satellite ie. which satellite has the max_slope



%Form Parity transformation Matrix
[Q,R] = qr(M);  %M is also called the G matrix.

%get submatrix %sub-matrix with lower n-4 rows of Q' is P matrix:

Qt = Q';
P = Qt(N-(N-4)+1:N,1:N);



%calculate range bias

%The maxslopsatellite'th column of P is what we want:

rbias = pbias/norm(P(:,maxslopsatellite));

%compute the position error vector induced by the bias of rbias in the max
%slope satellite
epsilon = zeros(N,1);  %this is already a column vector, dont need to transpose it
epsilon(maxslopsatellite) = rbias;
Poserr = AA*epsilon;

HPL = norm([Poserr(1) Poserr(2)]);




%% calculate VPL

 %% calculate V-slope
 for i_slope = 1:N
    V_SLOPE(i_slope) = sqrt(  (AA(3,i_slope)^2 )*(N-4)/(1 - BB(i_slope,i_slope)));
 end



V_SLOPE_Max = max(V_SLOPE);  %the value of slope_max
maxslopsatellite = find(SLOPE == SLOPE_Max) ;%i_slope indice to the max satellite ie. which satellite has the max_slope


rbias = pbias/norm(P(:,maxslopsatellite));

%compute the position error vector induced by the bias of rbias in the max
%slope satellite
epsilon = zeros(N,1);  %this is already a column vector, dont need to transpose it
epsilon(maxslopsatellite) = rbias;
Poserr = AA*epsilon;

VPL = norm(Poserr(3));


%Check for RAIM not being available:

if HPL > Alert_Limit
    
BadGeometry = 100;
else
    BadGeometry = 0;    
end



%Detect any errors


%Calculate parity vector


p = P*ResVec; 




%Use magnitude of p as test statistic. Alternatively could use 
%SSE = p'*p i think.

if N>5
r = norm([p(1),p(2)]);

else  %special case if N = 5. 
r = p;      
end

%SSE = abs(p)

%SSE = norm(p)

%SSE = p;

%compare to threshold 

if r > Td 
    RAIM_ALERT = 100;
else
    RAIM_ALERT = 0;
    
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



















