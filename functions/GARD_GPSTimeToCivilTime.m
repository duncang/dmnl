function [YY MM DD HH mm ss] = GARD_GPSTimeToCivilTime(GPS_WEEK, GPSsecs)
%
%

%GPSsecs = rem(dayssince/7,1)*7*24*60*60;
dayssince = (GPS_WEEK)*7 + GPSsecs / (24*60*60);

%days since GPS std epoch
%dayssince = JD - 2444244.5;
JD = dayssince + 2444244.5;

datenum_m = JD -  1721058.5;

[YY MM DD HH mm ss] = datevec(datenum_m);

