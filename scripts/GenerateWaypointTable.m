% generate waypoint list
% this script generates a table of waypoints (in LLH) based on the config
% file.  Essentially, it looks through the list of australian airports and
% matches the ICAO codes, extracts the LLH coordinates and outputs to a
% table.
% Written by Duncan Greer 25 January 2006
%
% $Id: GenerateWaypointTable.m 1879 2008-07-15 05:20:21Z n2523710 $
%

WptFile = 'data/long_flight_route.txt';
WptIndex = 0;

% load airports
load 'data/AustAirports.mat';


% load waypoint list
fid = fopen(WptFile,'r');

while(~feof(fid))
    WptIndex = WptIndex + 1;
    
    
    % read airport code
    Wptcode = fgetl(fid);
    if(length(Wptcode) ~= 4)
        continue;
    end
    
    % search for airport code in airports list
    for(AirpIndex = 1:length(AustAirports))
       Airpcode = char(AustAirports(AirpIndex,1:4));
       if(strcmp(Wptcode,Airpcode) == 1)
           WPTTable(WptIndex,:) = AustAirports(AirpIndex,5:7);
           break;
       end
    end
end
fclose(fid);

% save waypoints
save data/WPTTable WPTTable;