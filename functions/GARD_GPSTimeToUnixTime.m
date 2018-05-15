function unixtime = GARD_GPSTimeToUnixTime(GPS_WEEK, GPSsecs)
%
%

%GPSsecs = rem(dayssince/7,1)*7*24*60*60;
dayssince = (GPS_WEEK)*7 + GPSsecs / (24*60*60);

%days since GPS std epoch
%dayssince = JD - 2444244.5;
JD = dayssince + 2444244.5;




unixtime = (JD - 2440587.5) * 86400;

