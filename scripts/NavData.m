% this script compiles aa list of waypoints and navaids based on X-plane
% data.  The x-plane data is based on teh NIMA DAFIF navaid data and should
% be representative of true data.  Any errros in data should be reported to
% Robin Peel at x-plane.org.
% written by Duncan Greer 1 Aug 2005
%
% $Id: NavData.m 1884 2008-07-15 05:54:33Z n2523710 $
%

NavigData = 'data/airportdata/nav.dat';
FixData = 'data/airportdata/fix.dat';
AirpData = 'data/airportdata/apt.dat';

if(~exist('AustAirports'))
    AustAirports = ReadAirports(AirpData,-45,-10,110,155);
end

figure();
if ~exist('austcoast')
    load 'data/austcoast.dat';
end

% load the GRAS reference station locations
GRAS_GRS;

figure();
hold on;
plot(austcoast(:,1),austcoast(:,2));
plot(AustAirports(:,6)*180/pi,AustAirports(:,5)*180/pi,'r.');
plot(GRS_LLH(:,2),GRS_LLH(:,1),'kd');
grid on;
hold off;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

