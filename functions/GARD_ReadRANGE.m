function data = GARD_ReadRANGE(filename,type)
% function data = GARD_ReadBESTXYZ(filename,type)
% Written by Duncan Greer 25 May 2007
%
% Inputs
% ======
% filename - filename to read
% type - type of file format - 1 = binary, 2 = ascii, 3 = text
% 
% Outputs
% =======
% data - array containing data
%
% $Id: GARD_ReadBESTXYZ.m 1941 2008-07-28 06:16:46Z greerd $
%


%% Arguments
switch type
    case 1
        disp('Binary Not Supported');
        data = 0;
    case 2
        disp('Decoding ASCII format');
        data = Decode_RANGE_Ascii(filename);

    otherwise
        disp(sprintf('Type %d Not Supported',type));
        data = 0;
end


%% Ascii Data Files (commaa delimited Type 2)


function ascii_data = Decode_RANGE_Ascii(filename)

%   #RANGEA,SPECIAL,0,67.5,FINESTEERING,1591,430249.000,00000000,5103,5683;
%   27,
%   9,0,24554678.069,0.075,-139547873.976903,0.019,-2577.672,43.7,321.950,18109c04,
%   9,0,24554675.952,0.194,-108767911.289318,0.024,-2008.582,35.3,304.660,11309c0b,
%   25,0,20334591.275,0.029,-106858996.123086,0.004,-574.465,51.8,145.300,88109c24,
%   25,0,20334590.190,0.041,-83266746.279904,0.004,-447.641,49.0,139.700,81309c2b,
%   30,0,21099336.293,0.032,-121389930.085841,0.006,1548.605,50.9,329.500,08109c44,
%   30,0,21099332.261,0.057,-94618863.175043,0.006,1206.703,46.0,306.220,01309c4b,
%   27,0,26020835.331,0.159,-136740403.070072,0.022,-2702.863,37.1,216.870,18109c64,
%   27,0,26020834.425,0.147,-106550959.990348,0.029,-2106.137,37.8,216.880,11309c6b,
%   5,0,20918318.859,0.035,-120438675.749125,0.006,498.039,50.3,338.550,08109c84,
%   5,0,20918315.021,0.049,-93877626.659493,0.006,388.082,47.3,304.680,01309c8b,
%   31,0,25182340.718,0.111,-132334091.545232,0.015,900.258,40.2,280.974,18109ca4,
%   31,0,25182336.794,0.131,-103117462.565257,0.021,701.496,38.8,273.620,11309cab,
%   21,0,24337769.935,0.060,-127895844.081394,0.012,3015.094,45.6,283.130,08109cc4,
%   21,0,24337766.612,0.143,-99659091.233764,0.014,2349.422,38.0,270.600,01309ccb,
%   4,0,25581274.519,0.141,-134430506.199297,0.023,-2773.914,38.1,100.640,18109ce4,
%   4,0,25581273.361,0.366,-104751039.674362,0.046,-2161.492,29.8,100.640,11309ceb,
%   15,0,25790020.605,0.079,-135527467.015187,0.011,1746.695,43.9,112.880,08109d04,
%   15,0,25790019.762,0.143,-105605824.505846,0.013,1361.066,38.9,106.660,01309d0b,
%   10,0,22474705.743,0.063,-128617552.225893,0.016,-1513.340,45.1,336.060,18109d24,
%   10,0,22474702.594,0.093,-100250773.881796,0.017,-1179.227,41.7,301.720,11309d2b,
%   29,0,22581517.408,0.048,-129178846.608901,0.012,2069.031,47.4,329.500,08109d44,
%   29,0,22581514.147,0.090,-100688146.203068,0.012,1612.227,42.0,306.220,01309d4b,
%   12,0,20257692.923,0.031,-116967062.998755,0.006,-902.285,51.4,336.946,08109d64,
%   12,0,20257688.127,0.039,-91172470.829964,0.005,-703.082,49.4,306.220,01309d6b,
%   2,0,22366950.581,0.045,-117539117.211095,0.008,-2697.043,48.0,157.874,18109d84,
%   2,0,22366944.632,0.085,-91588900.668320,0.009,-2101.594,42.5,153.200,11309d8b,
%   25,0,20334606.955,0.010,-79797365.755143,0.005,-429.004,48.7,105.710,89c03dc4
%   *4c76da4f

inputfile = fopen(filename);


if inputfile < 0
    return ;
end

record = 0;





while ~feof(inputfile)
   
   line = fgetl(inputfile);
   
   
   
   if strfind(line,'#RANGEA')
       record = record + 1;
       

        ascii_data(record) = DecodeRANGERecord(line);
       

       
       if mod(record,1000) == 0
           disp(sprintf('Decoded record #%d',record));
       end
       
       
       
   else
       
      continue; 
   end
   

   
end


fclose(inputfile);



function rangerecord = DecodeRANGERecord(line)


% check if we actually have a range measurement
if strcmp(line(1:7),'#RANGEA') ~= 1
    rangerecord = 0;
    return;
end

% i is column index in line
lastdelim = 0;
field = 0;
for i=1:length(line)
   
    if line(i) == ','
        % reached end of field
        field = field+1;
        fielddata{field} = line(lastdelim+1:i-1);
        lastdelim = i;
    end

    if line(i) == ';'
        % reached end of header
        field = field+1;
        fielddata{field} = line(lastdelim+1:i-1);
        lastdelim = i;
        
        % record number of fields in header
        headerfields = field;
        
    end
    
    if line(i) == '*'
       % reached end of field
        % next field is the checksum
        field = field+1;
        fielddata{field} = line(lastdelim+1:i-1);
        lastdelim = i;
        
        % get the checksum
        checksum = line(lastdelim+1:length(line));
    end
end

% sanity check on header length
if headerfields ~= 10
    disp('Warning - header had unexpected number of fields');
end

% start processing the data
rangerecord.GPSWeek = str2num(fielddata{6});
rangerecord.GPSSec = str2num(fielddata{7});
rangerecord.NumberMeasurements = str2num(fielddata{11});

% check we have the expected number of fields
if size(fielddata,2) ~= 11 + rangerecord.NumberMeasurements*10
   disp('Warning - do not have expected number of fields'); 
end

for measurement = 1:rangerecord.NumberMeasurements

    % get a copy of the measurement data fields
    for j = 1:10
       measurementfield{j} = fielddata{11 + (measurement-1)*10 + j};
    end
    
    
    rangerecord.Measurements(measurement).usPRN = str2num(measurementfield{1});
    rangerecord.Measurements(measurement).glonass_freq = str2num(measurementfield{2});
    rangerecord.Measurements(measurement).dPseudorange = str2num(measurementfield{3});
    rangerecord.Measurements(measurement).dPseudorange_sigma = str2num(measurementfield{4});
    rangerecord.Measurements(measurement).carrier_phase = str2num(measurementfield{5});
    rangerecord.Measurements(measurement).carrier_phase_sigma = str2num(measurementfield{6});
    rangerecord.Measurements(measurement).fDoppler = str2num(measurementfield{7});
    rangerecord.Measurements(measurement).CNo = str2num(measurementfield{8});
    rangerecord.Measurements(measurement).locktime = str2num(measurementfield{9});
    rangerecord.Measurements(measurement).ch_tr_status = hex2dec(measurementfield{10});
    
    % determine the signal type - 0=L1, 9=L2P codeless, 14=L5C -see Novatel OM-20000094 table 72, p.400
    rangerecord.Measurements(measurement).signal_type = bitshift(bitand(uint64(rangerecord.Measurements(measurement).ch_tr_status),uint64(hex2dec('3E00000'))),-21);
    
    
end
