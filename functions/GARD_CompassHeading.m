function Heading = GARD_CompassHeading(mag, Xsf,Ysf,Xoff,Yoff,phi,theta)
%
% mag = 3x1 vector of raw magnetometer measuremnts
% Xsf, Ysf = X and Y magnetometer scale factors
% Xoff,Yoff = X and Y magnetometer offset factors
% phi, theta = roll and pitch angles for tilt compenstaion

% Calibration Correction

X = mag(1)*Xsf + Xoff;
Y = mag(2)*Ysf + Yoff;
Z = mag(3);


% Tilt Correction
Xh = X * cos(phi) + Y * sin(theta) * sin(phi) - Z * cos(theta) * sin(phi) ;
Yh = Y * cos(theta) + Z * sin(theta) ;


% Heading Calculation

Heading = atan2(Yh,Xh);


