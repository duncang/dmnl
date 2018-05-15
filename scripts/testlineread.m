
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
    
    
    rangerecord.Measurements(measurement).prn = str2num(measurementfield{1});
    rangerecord.Measurements(measurement).glonass_freq = str2num(measurementfield{2});
    rangerecord.Measurements(measurement).pseudorange = str2num(measurementfield{3});
    rangerecord.Measurements(measurement).pseudorange_sigma = str2num(measurementfield{4});
    rangerecord.Measurements(measurement).carrier_phase = str2num(measurementfield{5});
    rangerecord.Measurements(measurement).carrier_phase_sigma = str2num(measurementfield{6});
    rangerecord.Measurements(measurement).doppler = str2num(measurementfield{7});
    rangerecord.Measurements(measurement).CNo = str2num(measurementfield{8});
    rangerecord.Measurements(measurement).locktime = str2num(measurementfield{9});
    rangerecord.Measurements(measurement).ch_tr_status = hex2dec(measurementfield{10});
    
    % determine the signal type - 0=L1, 9=L2P codeless, 14=L5C -see Novatel OM-20000094 table 72, p.400
    rangerecord.Measurements(measurement).signal_type = bitshift(bitand(uint64(ch_tr_status),uint64(hex2dec('3E00000'))),-21);
    
    
end


