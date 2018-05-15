function Output = GARD_FlyApproach(Arguments)
% function Output = GARD_FlyApproach(Arguments)
%
% Written by Duncan Greer and Troy Bruggemann 
% $Id: GARD_FlyApproach.m 1850 2008-07-14 04:52:47Z greerd $
%

% The approach sequence is 
% 1 - approach.IAF
% 2 - approach.IF
% 3 - approach.FAF
% 4 - approach.MAPWP
% 5 - approach.MAPHP

% if we are on WP1, then we are navigating to the IAF from our present
% position.  should aim to hit the IAF at the correct track (i.e. the track
% from IAF to IF
% 
%

% CurrentWaypoint is an index into teh approach waypoint table (WPTTable)

CurrentWaypoint  = Arguments(1);
CurrentLatitude = Arguments(2);
CurrentLongitude = Arguments(3);
CurrentHeight = Arguments(4);

if ~exist('approach','var');
    load data/approach.mat
end

% find the waypoint index of the approach fixes
for i=1:length(approach.WPTName)
    if strcmp(approach.WPTName(i,:),approach.IAF)
        IAF_WPIndex = i;
    end
    
    if strcmp(approach.WPTName(i,:),approach.IF)
        IF_WPIndex = i;
    end
    
    if strcmp(approach.WPTName(i,:),approach.FAF)
        FAF_WPIndex = i;
    end
    
    if strcmp(approach.WPTName(i,:),approach.MAPWP)
        MAPWP_WPIndex = i;
    end
    
    if strcmp(approach.WPTName(i,:),approach.MAPHP)
        MAPHP_WPIndex = i;
    end
end
    
% initialisation
if(CurrentWaypoint == 0)
    CurrentWaypoint = IAF_WPIndex;
end



CurrentWaypointName = approach.WPTName(CurrentWaypoint,:);



% get current waypoint coords
WptLat = approach.WPTTable(CurrentWaypoint,1);
WptLong = approach.WPTTable(CurrentWaypoint,2);
WptHeight = approach.WPTTable(CurrentWaypoint,3);

Holding = 0;

if CurrentWaypoint == IAF_WPIndex
    StartWPLat = CurrentLatitude;
    StartWPLong = CurrentLongitude;
    NextWPIndex = IF_WPIndex;
else
    
    switch CurrentWaypoint
        case IF_WPIndex
            LastWPIndex = IAF_WPIndex;
            NextWPIndex = FAF_WPIndex;
        case FAF_WPIndex
            LastWPIndex = IF_WPIndex;
            NextWPIndex = MAPWP_WPIndex;
        case MAPWP_WPIndex
            LastWPIndex = FAF_WPIndex;
            NextWPIndex = MAPHP_WPIndex;
        case MAPHP_WPIndex
            LastWPIndex = MAPWP_WPIndex;
            NextWPIndex = IAF_WPIndex;
            Holding = 1;
    end
    
    if ~exist('LastWPIndex','var')
        LastWPIndex = IAF_WPIndex;
    end
    
    StartWPLat = approach.WPTTable(LastWPIndex,1);
    StartWPLong = approach.WPTTable(LastWPIndex,2);
end

[RequiredTrack,TrackDistance] =  CalculateGC(StartWPLat,StartWPLong,WptLat,WptLong);
[Bearing_to_WP,Distance_to_WP] = CalculateGC(CurrentLatitude,CurrentLongitude,WptLat,WptLong);
Waypoint = CurrentWaypoint;



% determine if we need to sequence to next WP
if Distance_to_WP < (0.5*1852) && CurrentWaypoint ~= MAPWP_WPIndex
    Waypoint = NextWPIndex;
end

if Distance_to_WP < (0.05*1852) && CurrentWaypoint == MAPWP_WPIndex
    Waypoint = NextWPIndex;
end



Track_Error = RequiredTrack - Bearing_to_WP;
XTE = Distance_to_WP * sin(Track_Error);

%Distance_to_MAPWP = 0;
Glide_path_height = 0;

%MAPWptLat = approach.WPTTable(MAPWP_WPIndex,1);
%MAPWptLong = approach.WPTTable(MAPWP_WPIndex,2);
%[Bearing_to_MAPWP,Distance_to_MAPWP] = CalculateGC(CurrentLatitude,CurrentLongitude,MAPWptLat,MAPWptLong);

[Bearing_to_GPI,Distance_to_GPI] = CalculateGC(CurrentLatitude,CurrentLongitude,approach.GPI(1),approach.GPI(2));


Glide_path_height = (Distance_to_GPI*tan(approach.GPA)+approach.AD_ELEV);

if CurrentWaypoint == MAPHP_WPIndex
    % if we are tracking the Missed approach holding point, climb.
    Glide_path_height = approach.WPTTable(MAPHP_WPIndex,3);
end

if CurrentWaypoint == IAF_WPIndex
    % if we are tracking the Missed approach holding point, climb.
    Glide_path_height = approach.WPTTable(IAF_WPIndex,3);
end


Output = [Bearing_to_WP,Distance_to_WP,Waypoint,XTE,Holding,Track_Error,Glide_path_height,Distance_to_GPI];

