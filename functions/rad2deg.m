function D=rad2deg(R)
%global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;
%Version 1.00
%
%DEG2RAD Converts angles from degrees to radians
%
%  rad = DEG2RAD(deg) converts angles from degrees to radians.
%
%  See also RAD2DEG, DEG2DMS, ANGLEDIM, ANGL2STR
%
%  Copyright 1996-2003 The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1883 $    $Date: 2008-07-15 15:53:55 +1000 (Tue, 15 Jul 2008) $


if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(R)
     warning('Imaginary parts of complex ANGLE argument ignored')
     R = real(R);
end

D = R*180/pi;
