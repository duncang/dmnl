function [Airports] = ReadAirports(Filename,LatMin,LatMax,LongMin,LongMax)
% function [Airports] = ReadAirports(Filename,LatMin,LatMax,LongMin,LongMax)
% Reads in the apt.dat file given by Filename and returns a vector of
% airports in LLH (rad,rad,m-AMSL) which are within the geographical limits
% provided by LatMin,LatMax,LongMin,LongMax
%
% this program reads the first runway centre point provided in the file
%
% Written by Duncan Greer January 2006
%
% $Id: ReadAirports.m 1883 2008-07-15 05:53:55Z n2523710 $
%

% open file for reading
fid = fopen(Filename,'r');

NumAirports = 0;
d2r = pi/180;
r2d = 180/pi;
f2m = 1/3.2;

% cycle through line by line
% Line 1 - New airport description
% Line 10 - runway descriptoin
while(~feof(fid))
    newline = fgetl(fid);
    if(length(newline) < 2)
        continue;
    end
    if(str2num(newline(1:2)) == 1)
        % start of a new airport line
        

        Elevation  = str2num(newline(3:8));
        ICAO = newline(14:17);
        
        % get the first runway location
        newline = fgetl(fid);
        while(str2num(newline(1:2)) ~= 10)
            newline = fgetl(fid);
        end
        
        Lat = str2num(newline(4:13));
        Long = str2num(newline(15:25));
        
        if(Lat > LatMin && Lat < LatMax && Long > LongMin && Long < LongMax)

            NumAirports = NumAirports+1;

            % construct airport vector
            NewAirport = [double(ICAO),Lat*d2r,Long*d2r,Elevation*f2m];
        
            Airports(NumAirports,:) = NewAirport;
        end
    end
    

    
    
    
end
fclose(fid);
