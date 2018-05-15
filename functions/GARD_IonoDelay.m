function [IonoDelay] = GARD_IonoDelay(Elevation, Azimuth, UserLatitude, UserLongitude, IonoTECMap, IonoLatitude, IonoLongitude)
% function [IonoDelay] = GARD_IonoDelay(Elevation, Azimuth, UserLatitude, UserLongitude, IonoTECMap)
%  
% Written by Duncan Greer 28 Oct 2005
%
% This function returns the 'true' ionospheric delay in metres for a given satellite
% slant vector (defined by the Azimuth and Elevation) calculated from the 
% post-processed ionospheric total-electron-count TEC maps provided by the IGS.
%
% * Elevation and Azimuth are the LTP line-of-sight values in radians.
% * UserLatitude and Longitude are self-explanatory
% * IonoTECMap - TEC map for the desired time epoch.  Note that ReadIonex
%      returns a vector of TEC Maps, you need to provide the correct map 
%      for the time.
%
% >>>>> COMPLETION STATUS - UNCOMPLETED <<<<<<
% >>>>> VALIDATION STATUS - UNVALIDATED <<<<<<
%
% =========================================================================



% calculate sub-iono point
[lamda_i, phi_i] = getsubionopoint(UserLongitude, UserLatitude, Elevation, Azimuth);

% lookup VTEC value from TEC Map
 long_index = find(IonoLongitude > lamda_i*180/pi, 1, 'first'); % return the first iono tec map bin that matches the iono longitude
 lat_index = find(IonoLatitude > phi_i*180/pi, 1, 'last'); % return the first inoo tec map bin that matches the iono latitude
 VTEC = IonoTECMap(lat_index,long_index);

% convert to Slant TEC
z = pi/2 - Elevation;
STEC = vtec2stec(VTEC,z,6371.8,6371.8,450);

%  calculate iono delay using 16.2cm per TECU
IonoDelay = STEC * 0.162;
