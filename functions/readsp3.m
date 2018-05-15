function [NumberSVs, VehicleIDs, EpochNumber, GPSTime, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data] = readsp3(SP3FILE);
% [NumberSVs, VehicleIDs, EpochNumber, GPSTime, SV_X_Data, SV_Y_Data, SV_Z_Data, SV_T_Data] = readsp3(SP3FILE) 
% Written by Duncan Greer 23 May 2005 
% $Id: readsp3.m 1883 2008-07-15 05:53:55Z n2523710 $
% function to read the contents of an SP3 GPS data file and return as an
% array
%

% open file for reading
DataFile = fopen(SP3FILE, 'r');
DataLineNumber = 0;
EpochNumber = 0;

% initialise output data with 31 satellites

SV_X_Data(1,1:31) = 0;
SV_Y_Data(1,1:31) = 0;
SV_Z_Data(1,1:31) = 0;
SV_T_Data(1,1:31) = 0;
    
% loop until we reach the end of the file
while ~feof(DataFile)
    
    % read a line of data
    DataLine = fgetl(DataFile);
    DataLineNumber = DataLineNumber + 1;
    
    % read the number of satellites
    if DataLineNumber == 3
        NumberSVs = str2num(DataLine(5:6));
        LogString = sprintf('There is data for %d SVs in this file.', NumberSVs);
        disp(LogString);
    end
    
    % check if we have reached the start of the SV position data
%     if DataLineNumber == 23
%         EpochNumber = EpochNumber + 1;
%         
%     end
    


    if DataLineNumber > 22 & DataLine(1) == '*'
        % we have reached a new record
        EpochNumber = EpochNumber + 1;
        SVNumber = 0;
        
        % get the time and date for this record
        % Format: '*  YYYY MM DD HH MM SS.ssssssss' eg. '*  2005  1  5 23 45   .00000000'
        
        year = str2num(DataLine(4:7)) - 2000;
        month = str2num(DataLine(9:10));
        day = str2num(DataLine(12:13));
        hour = str2num(DataLine(15:16));
        minute = str2num(DataLine(18:19));
        second = str2num(DataLine(21:31));
        
        [GPSWeek(EpochNumber),GPSSecs(EpochNumber)] = ftime([year month day hour minute second]);
        GPSTime(EpochNumber,:) = [GPSWeek(EpochNumber),GPSSecs(EpochNumber)];
        
        % read data for number of satellites
        while SVNumber < NumberSVs
            % read a line of data
            DataLine = fgetl(DataFile);
            DataLineNumber = DataLineNumber + 1;
            SVNumber = SVNumber + 1;
            %VehicleID = str2num(DataLine(2:4));
            VehicleID = str2num(DataLine(3:4));   %note just put 3:4 so it could read it in properly- troy 18.3.08 
            
            SV_X_Data(EpochNumber,VehicleID) = str2num(DataLine(5:18)) * 1000;
            SV_Y_Data(EpochNumber,VehicleID) = str2num(DataLine(19:32)) * 1000;
            SV_Z_Data(EpochNumber,VehicleID) = str2num(DataLine(33:46)) * 1000;
            SV_T_Data(EpochNumber,VehicleID) = str2num(DataLine(47:60)) * 1e-6;
            VehicleIDs(VehicleID) = 1;
            
        end
        
        %LogLine = sprintf('Processed record number %d',EpochNumber);
        %disp(LogLine);
        
    end
        
    
end



