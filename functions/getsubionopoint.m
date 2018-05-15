function [lamda_i, phi_i] = getsubionopoint(lamda_u, phi_u, E, A);

global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;
% function calculates the sub ionospheric longitude and latitude for a
% visible satellite from the user position (lamda_u, phi_u) given satellite
% elevation (E) and azimuth (A).

% Note:  This function is based on the ICD200 and has not been validated.

%   lamda_u - user longitude in radians
%   phi_u - user latitude in radians
%   E - SV Elevation in radians (LTP)
%   A - SV Azimuth in radians (LTP)

% written by Duncan Greer 15 Aug 2005

% convert longitude and latitude to semi-circles
lamda_u = lamda_u / GPS_PI;
phi_u = phi_u / GPS_PI;
E = E / GPS_PI;
A = A / GPS_PI;


% calculate the Earth-centred angle, psi
psi = 0.0137 / (E + 0.11) - 0.022; % semi-circles


% calculate the subionospheric latitude, phi_i
phi_i = phi_u + psi * cos(A);

% limit phi_i
if phi_i > 0.416
    phi_i = 0.416;
else if phi_i < -0.416
        phi_i = -0.416;
    end
end

% geodetic longitude of the earth projection of hte ionopsheric
% intersection point (subionospheric longitude)
lamda_i = lamda_u + psi * sin(A)/cos(phi_i); % semi-circles

% convert output longitude and latitude to radians
lamda_i = lamda_i * GPS_PI;
phi_i = phi_i * GPS_PI;
