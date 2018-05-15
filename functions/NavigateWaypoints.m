function Output = NavigateWaypoints(Arguments)
% function Output = NavigateWaypoints(Arguments)
%
% performs waypoint navigation.
%
% Written by Duncan Greer 25 January 2006
% 
% CurrentWaypoint  = Arguments(1);
% CurrentLatitude = Arguments(2);
% CurrentLongitude = Arguments(3);
% CurrentHeight = Arguments(4);
%
% $Id: NavigateWaypoints.m 1883 2008-07-15 05:53:55Z n2523710 $
%

% load waypoint list
if(~exist('WPTTable'))
    load 'data\long_flight\WPTTable.mat';
end

CurrentWaypoint  = Arguments(1);
CurrentLatitude = Arguments(2);
CurrentLongitude = Arguments(3);
CurrentHeight = Arguments(4);

% ensure it stays in scope
%global WPTTable;

if(CurrentWaypoint == 0)
    CurrentWaypoint = 1;
end

StartApproach = 0;

% get current waypoint coords
WptLat = WPTTable(CurrentWaypoint,1);
WptLong = WPTTable(CurrentWaypoint,2);
WptHeight = WPTTable(CurrentWaypoint,3);

[Bearing,Distance] = CalculateGC(CurrentLatitude,CurrentLongitude,WptLat,WptLong);
Waypoint = CurrentWaypoint;

% if we are within 500m of this waypoint, declare it captured 
if(Distance < 500)
    % check if we are on our last waypoint
    if(CurrentWaypoint == length(WPTTable))
        Bearing = 000;
        Distance = 000;
        Waypoint = 1;
        StartApproach = 1;
    else
        WptLat = WPTTable(CurrentWaypoint+1,1);
        WptLong = WPTTable(CurrentWaypoint+1,2);
        WptHeight = WPTTable(CurrentWaypoint+1,3);
        [Bearing,Distance] = CalculateGC(CurrentLatitude,CurrentLongitude,WptLat,WptLong);
        Waypoint = CurrentWaypoint + 1;
        StartApproach = 0;
    end
end

Output = [Bearing,Distance,Waypoint,StartApproach];
