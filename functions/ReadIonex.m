function [IonoTECMap, GPSTime_Obs, Longitude_out, Latitude_out] = ReadIonex(Filename);
% Fucntion reads the IONEX file (Version 1) given by Filename and returns the IonoTEC
% maps as a grid
% Written by Duncan Greer 2 August 2005 for GARDSim
%
% See http://igscb.jpl.nasa.gov/components/compindex.html for more info.

% open the file for reading
IonoData = fopen(Filename,'r');

% TODO: file error handling

% get a line
Newline = fgetl(IonoData);

Epoch = 0;

while Newline ~= -1
    
    % check if this is the start of a TEC map record - NOte, there are also
    % RMS Maps which are not yet supported.
    if Newline(61:76) == 'START OF TEC MAP'
        
        Epoch = Epoch + 1;
        
        % get the next line
        Newline = fgetl(IonoData);
        % get the time of this record
        yyyy = str2num(Newline(3:6));
        mm = str2num(Newline(11:12));
        dd = str2num(Newline(17:18));
        hh = str2num(Newline(23:24));
        mm = str2num(Newline(29:30));
        ss = str2num(Newline(35:36));

        GPSTime_Obs(Epoch) = ftime([yyyy mm dd hh mm ss]);

        % get the next line
        Newline = fgetl(IonoData);
            
        LatitudeIndex = 0;
        % loop until we get to the end of this record
        while ~strcmp(Newline(61:74),'END OF TEC MAP')
            % record index
            LatitudeIndex = LatitudeIndex + 1;
            

            % get the latitude block
            Latitude(LatitudeIndex) = str2num(Newline(4:8));
            
            % get the longitude range and increment
            LongitudeStart(LatitudeIndex) = str2num(Newline(9:14));
            LongitudeEnd(LatitudeIndex) = str2num(Newline(15:20));
            LongitudeInc(LatitudeIndex) = str2num(Newline(22:26));
            % get the height
            Height(LatitudeIndex) = str2num(Newline(28:32));

            % calculate the number of measurements
            NumberMeasurements(LatitudeIndex) = ((LongitudeEnd - LongitudeStart) / LongitudeInc) + 1;

            % calculate hte number of lines (16 measurements per line)
            NumberLines = ceil(NumberMeasurements / 16);


            % read data in
            for Line = 1:NumberLines
                % get the next line
                Newline = fgetl(IonoData);

                for MeasurementSlotNumber = 1:16
                    MeasurementNumber = MeasurementSlotNumber + (Line - 1) * 16;
                    if (Line == NumberLines(LatitudeIndex)) && (MeasurementSlotNumber >= mod(NumberMeasurements(LatitudeIndex),16))
                        % this line wont be full of measurements
                        break;
                    end

                    Measurement(LatitudeIndex,MeasurementNumber) = str2num(Newline(MeasurementSlotNumber*5 - 4:MeasurementSlotNumber*5)) / 10;

                end
            end

            % start next latitude block
            % get the next line
            Newline = fgetl(IonoData);
            
        end % while not at end of epoch
    end % if at start of new epoch
    
    % get a line
    Newline = fgetl(IonoData);
    
    IonoTECMap(:,:,Epoch) = Measurement;
    Longitude_out(Epoch,:) = [LongitudeStart(1):LongitudeInc(1):LongitudeEnd(1)];
    Latitude_out(Epoch,:) = Latitude;
end % while not EOF




