function D_IONO = ionomodel(GPSTime, Xu, Su, alpha, beta)
% GPS ICD200 Single Frequency Ionospheric Model Correction
% Written by Duncan Greer for the GARDSim
%
% This function implements the ionospheric model for single frequency GPS
% pseudorange measurements based on the broadcast model parameters.  See
% ICD200 Figure 20-4 Page 126.
% 
%   Function returns the ionospheric delay (D_IONO) in metres.  
%
%   D_IONO = ionomodel(GPSTime, Xu, Su, alpha, beta)
%
%   - GPSTime is the receiver calculated GPS Time
%   - Xu is the user position in ECEF
%   - Su is the satellite measurement position in ECEF
%   - alpha and beta are the coefficients of the amplitude and period of the
%   delay model given by the GPS nav message. They are a 4 x 1 vector
%   - 
global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;

speedoflight = 2.99792458e8;  % speed of light m/s

% get user geodetic latitude and longitude
[phi_u,lamda_u,h_u] = ECEF2LLH(Xu);

% calculate satellite line of sight vector
SV_LOS_ecef = (Su - Xu)';

% generate rotation matrix
TMatrix_ECEF2ENU = T_ECEF2ENU(lamda_u,phi_u);

% convert to LTP coordinates
SV_LOS_enu = TMatrix_ECEF2ENU * SV_LOS_ecef;

% calculate azimuth and elevation
A = atan2(SV_LOS_enu(1),SV_LOS_enu(2)) / pi; % semi-circles
E = atan2(SV_LOS_enu(3),sqrt(SV_LOS_enu(1)^2 + SV_LOS_enu(2)^2)) / pi; % semi-circles


% if the elevation is less than 0, return a delay of 0
if E < 0
    D_IONO = 0;
    return;
end

% calculate the Earth-centred angle, psi
psi = 0.0137 / (E + 0.11) - 0.022; % semi-circles

% calculate the subionospheric latitude, phi_i
phi_i = phi_u/pi + psi * cos(A*pi);

% limit phi_i
if phi_i > 0.416
    phi_i = 0.416;
else if phi_i < -0.416
        phi_i = -0.416;
    end
end

% geodetic longitude of the earth projection of hte ionopsheric
% intersection point (subionospheric longitude)
lamda_i = lamda_u/pi + psi * sin(A*pi)/cos(phi_i*pi); % semi-circles

% geomagnetic latitude of the earth projection of the ionophseric intersection point
phi_m = phi_i + 0.064 * cos(lamda_i - 1.617); % semi-circles

% calculate local time
t_local = 4.32e4 * lamda_i + GPSTime; % seconds

% limit local time to +- 86400 seconds (24 hours)
if t_local < -86400
    %t_local = t_local + 86400;
    t_local = mod(t_local,-86400);
else if t_local >= 86400
    %    t_local = t_local - 86400;
    t_local = mod(t_local,86400);
    end
end



% calculate the slant factor, F
F = 1 + 16 * (0.53 - E) ^3;



% calculate the amplitude and phase of the delay
AMP = alpha(1) + alpha(2) * phi_m + alpha(3) * phi_m^2 + alpha(4) * phi_m^3;

if AMP < 0
    AMP = 0;
end

PER = beta(1) + beta(2)*phi_m + beta(3) * phi_m^2 + beta(4) * phi_m^3;

if PER < 72000
    PER = 72000;
end

% calculate x - what is x?!??
x = 2 * pi * (t_local - 50400) / PER;

% calculate the delay in seconds using x
if abs(x) < 1.57
    T_IONO = F * (5e-9 + AMP*(1 - (x^2)/2 + (x^4)/24));
else 
    T_IONO = F * 5e-9;
end

% convert seconds to metres by multiplying by speed of light.
D_IONO = T_IONO * speedoflight;

