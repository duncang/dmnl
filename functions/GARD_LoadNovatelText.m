function [data] = GARD_LoadNovatelText(filename)
% function [data] = GARD_LoadNovatelText(filename)
% Returns observations from novatel receiver in a data structure.
% Written by Duncan Greer 7 May 2007
%
% $Id: GARD_LoadNovatelText.m 1850 2008-07-14 04:52:47Z greerd $
%

% the header and data part are separated by a ';' character
% the data fields are separated by a ',' character
% the packet type is given by a '#PACKET' e.g. '#RANGEA'


%% example range packet
% #RANGEA,COM1,0,57.5,COARSESTEERING,1416,345668.000,00480000,5103,2580;
%  6,
%  1,0,25717905.036,0.108,-135148502.140034,0.010,2508.078,46.1,36.790,08009c04,
%  23,0,22304758.188,0.120,-117212294.257946,0.010,5296.723,47.1,23.620,18009c44,
%  4,0,26368940.202,0.113,-138569714.998066,0.013,4307.789,45.8,36.180,18009c64,
%  25,0,23477788.123,0.120,-123376608.427005,0.010,2201.750,46.6,26.640,08009c84,
%  16,0,23491342.725,0.108,-123447838.307516,0.010,2327.582,46.7,32.400,08009d04,
%  20,0,22205662.225,0.104,-116691541.888192,0.008,2991.227,47.1,31.360,18009da4
%  *b632e08d

% open the file
inputfd = fopen(filename,'r');

data = 0;

linenumber = 0;
rangepackets = 0;

while ~feof(inputfd)
    line = fgetl(inputfd);
    linenumber = linenumber+1;
    % get the header and determine string type
    [token, rest] = strtok(line,',');
    
    if strcmp(token,'#RANGEA')
       % decode 
       
        rangepackets = rangepackets+1;
        [rangedata,NumberMeasurements,GPSTime_Week,GPSTime_Sec] = DecodeRange(line);
        
        
        disp(sprintf('Found Range Packet with %d Measurements',NumberMeasurements));
        
        if NumberMeasurements > 0
            MaxPRN = size(rangedata,1);
            outputrangedata(rangepackets,1:MaxPRN,1:10) = rangedata;
        end
        outputrangetime(rangepackets,:) = [GPSTime_Week, GPSTime_Sec];
        clear rangedata;
    end
    
    
    
    
end



% close the file
fclose(inputfd);

disp(sprintf('Processed %d Lines',linenumber));
disp(sprintf('Found %d Range Packets',rangepackets));

data.RANGE.GPSTime_Week = outputrangetime(:,1);
data.RANGE.GPSTime_Sec = outputrangetime(:,2);
data.RANGE.LockTime = outputrangedata(:,:,9);
data.RANGE.C1 = outputrangedata(:,:,3);
data.RANGE.C1_sigma = outputrangedata(:,:,4);
data.RANGE.L1 = outputrangedata(:,:,5);
data.RANGE.L1_sigma = outputrangedata(:,:,6);
data.RANGE.D1 = outputrangedata(:,:,7);
data.RANGE.S1 = outputrangedata(:,:,8);


%% sub function to decode #RANGEA data
function [rangedata,N,GPSTime_Week, GPSTime_Sec] = DecodeRange(string)

% separate header from data
[header,packetdata] = strtok(string,';');

% decode the header data to get GPSTime of measurement
[GPSTime_Week, GPSTime_Sec] = DecodeHeader(header);


% separate checksum from data
[packetdata,checksum] = strtok(packetdata,'*');

% the range packet has 10 fields for each measurement.  the first
% field contains the number of measurements
[nummeas,remain] = strtok(packetdata,',');
NumberMeasurements = str2num(nummeas);

rangedata = zeros(1,10);

% loop through measurements
for MeasurementNumber = 1:NumberMeasurements

    for MeasurementType = 1:10
        [measurement,remain] = strtok(remain,',');

        if MeasurementType == 1
            PRN = str2num(measurement);
            rangedata(PRN,MeasurementType) = PRN;
        elseif MeasurementType == 10
            rangedata(PRN,MeasurementType) = hex2dec(measurement);
        else
            rangedata(PRN,MeasurementType) = str2num(measurement);
        end
    end
end



N = NumberMeasurements;
        
% #RANGEA,COM1,0,57.5,COARSESTEERING,1416,345668.000,00480000,5103,2580
%% decode header data
function [GPSTime_Week, GPSTime_Sec] = DecodeHeader(string)
%[a1 a2 a3 a4 a5 a6 a7 a8 a9 a10] = strread(string,'%s %s %d %f %s %d %f %d %d %d','delimiter',',')

remain = string;

for field = 1:10
       
    [data,remain] = strtok(remain,',');
    
    if field == 6
        GPSTime_Week = str2num(data);      
    end
    
    if field == 7
        GPSTime_Sec = str2num(data);
    end
end



