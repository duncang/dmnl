function R=deg2rad(D)
global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;
%Version 1.00

%DEG2RAD Converts angles from degrees to radians
%
%  rad = DEG2RAD(deg) converts angles from degrees to radians.
%
%  See also RAD2DEG, DEG2DMS, ANGLEDIM, ANGL2STR

%  Copyright 1996-2003 The MathWorks, Inc.
%  Written by:  E. Byrns, E. Brown
%   $Revision: 1874 $    $Date: 2008-07-15 14:42:16 +1000 (Tue, 15 Jul 2008) $


if nargin==0
	error('Incorrect number of arguments')
elseif ~isreal(D)
     warning('Imaginary parts of complex ANGLE argument ignored')
     D = real(D);
end

R = D*GPS_PI/180;
