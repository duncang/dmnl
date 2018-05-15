function [GPS_WEEK, GPSsecs] = ftime(civil)
% [GPS_WEEK, GPSsecs] = ftime(civil)
%
% input - civil date in [year month day hour minute second] format.
% only works for data year 2000 or later
% $Id: ftime.m 1874 2008-07-15 04:42:16Z n2523710 $
%


Y = 2000 + civil(1);
M = civil(2);
D = civil(3);
UT = civil(4) + civil(5)/60 + civil(6)/3600;



if M <= 2
   y = Y-1;
   m = M+12;
end

if M > 2
  y = Y;
  m = M;
end;


JD = fix(365.25*y) + fix(30.6001*(m+1)) + D + UT/24 + 1720981.5;

%JD is working properly


hh = civil(4);
mm = civil(5);
ss = civil(6);

%alternative way to calculate.

%JD = (367.0*y - floor((7.0*(y + floor((m + 9.0) / 12.0)))/4.0)...
%				+ floor((275.0 * m)/9.0) + D + 1721013.5...
%				+ ((((ss/60.0) + mm)/60.0 + hh)/24.0));


%days since GPS std epoch
dayssince = JD - 2444244.5;
GPS_WEEK = fix(dayssince/7);
GPSsecs = rem(dayssince/7,1)*7*24*60*60;

