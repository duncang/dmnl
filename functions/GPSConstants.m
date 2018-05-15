%GPSConstants for GARDSim

global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;


GPS_PI = 3.1415926535898; %ICD value of PI
%OMEGAedot = 7.2921151467e-5; % Earth rotation rate in radians per second

OMEGAedot = 7.292115e-5; % this is the value aerosim uses, changed 11.1.08.. Earth rotation rate in radians per second

mu = 3.986005e14;  % WGS-84 valeu of the earths universal gravitational parameter in m^3/s^2
Earthradius = 6378136;  %m
Speedoflight = 2.99792458e8;  % speed of light m/s
c = 2.99792458e8; % speed of light m/s

F = -4.442807633e-10; % a random number from the ICD page 88


L1_f = 1575.42e6; %Hz
L2_f = 1227.6e6; %Hz

gamma = (L1_f/L2_f)^2;  % unitless

L1_Wavelength = Speedoflight/L1_f; %Metres



