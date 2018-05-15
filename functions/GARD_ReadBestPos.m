function data = GARD_ReadBestPos(filename)
% function data = GARD_ReadBestPos(filename)
% Reads the data stored in filename to a data structure
% Written by Duncan Greer, 1 March 2007
% 
% Last Update: $Id: GARD_ReadBestPos.m 1850 2008-07-14 04:52:47Z greerd $
%
% Column Data:
% 1 - Packet Type (#BESTPOS)
% 2 - Port
% 3 - 
% 4 - 
% 5 - 
% 6 - GPS Week
% 7 - GPS Second of Week
% 8 - 
% 9 - RX Status?
% 10 - 
% 11 - Solution Status
% 12 - Solution Type
% 13 - Latitude (deg)
% 14 - Longitude (deg)
% 15 - Height (m)
% 16 - Undulation (m)
% 17 - Datum
% 18 - Latitude Variance/DOP (m)
% 19 - Longitude Variance/DOP (m)
% 20 - Height Variance/DOP (m) 
% 21 - Base Station ID
% 22 - Differential Age (sec)
% 23 - Differential Latency (sec)
% 24 - Satellites Tracked
% 25 - Satellites Used
% 26 - RTK Satellites Tracked ??
% 27 - RTK Satellites Used ??
% 28 - RTK .. ???
% 29 - RTK .. ???
% 30 - ???


%
% #BESTPOSA,COM1,0,79.5,FINESTEERING,1416,91867.000,00000000,4ca6,
% 2580;SOL_COMPUTED,SINGLE,-27.47737340335,153.02715632184,51.4889,
% 40.7703,WGS84,1.9235,1.3927,3.2740,"",0.000,0.000,10,10,0,0,0,0,0,0*e07ae1a2

fid = fopen(filename);
data = textscan(fid,'%s%s%d%f%s%d%f%d%d%s%s%s%f%f%f%f%s%f%f%f%q%f%f%d%d%d%d%d%d%d%*[^\n]','delimiter',',');
fclose(fid);