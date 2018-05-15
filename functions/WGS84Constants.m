% WGS-84 Earth Ellipsoid Constants
% Added by Duncan Greer 22 January 2006
%
% $Id: WGS84Constants.m 1944 2008-07-28 07:51:03Z greerd $
%

% declare as global
global a f e2 e nm2m m2nm;

a = 6378137.0;   % semi-major axis (metres)
f = 1/298.2572; % flattening
e2 = f * (2-f); % eccentricity squared
e = sqrt(e2);   % first eccentricity

nm2m = 1852;  % nautical miles to metres conversion eg. 1 nm in metres = 1nm * nm2m = 1852m;
m2nm = 1/nm2m; % metres to nautical miles eg.  1nm = 1m * m2nm;