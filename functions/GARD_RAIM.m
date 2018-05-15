
function [BadGeometry, RAIM_ALARM, SLOPE_Max, r, rd, ARP] = GARD_RAIM(N,PFalseAlarm,SigmaS,Alarm_Limit,ResVec,M);
%Version 1.00
%Troy Bruggemann 30 June 2005

%Does RAIM FDI And FDE using Parkinson's PR residual method.
% INPUT
% N - Number of Satellites to do RAIM solution with
% PFalseAlarm - Probability of False Alarm
% SigmaS - Expected Standard Deviation of Pseudorange errors
% Alarm_Limit - 3D Protection limit for given phase of flight
% ResVec - [1..N] Vector of Pseudorange residuals from least squares (m)
% M = M matrix from least squares
%=======================================================================
% OUTPUT
% BadGeometry - Flag indicating poor satellite-user geometry
% RAIM_ALARM - Flag indicating poor position solution integrity
% SLOPE_Max - Value of Maximum Slope
% r - test statistic
% rd - calculated threshold value from Chi square distribution
% ARP - approximate radial protected


%check that Number of satellites enough to do RAIM FDI

if N < 5;
    error('Not enough satellites for RAIM solution')
end

%RAIM section , using Least-Squares-Residuals Method
%determine the threshold using chi-square statistics.
%Chi-square degrees of freedom , d = n - 4 = 2.

d = N-4;

a = chi2inv(1-PFalseAlarm,d);

rd = sqrt(a*(SigmaS)^2/(N-4));

%Calculate ARP
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

SLOPE_Max = max(SLOPE);

maxslopsatellite = find(SLOPE == SLOPE_Max) ;%i_slope indice to the max satellite

ARP = SLOPE_Max*rd;

%set alarm limit above ARP:
%Alarm_Limit(i) = ARP(i); %1.7 rule of thumb
%Alarm_Limit = 339; %339 metres

%Check for bad geometry

if ARP > Alarm_Limit

    BadGeometry = 100;
else
    BadGeometry = 0;

end

SSE = ResVec'*ResVec;

%Form test statistic
r = sqrt((SSE/(N-4)));

if SLOPE_Max*r > Alarm_Limit
    RAIM_ALARM = 100;
else
    RAIM_ALARM = 0;

end






