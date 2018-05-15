function data = GARD_ReadBESTXYZ(filename,type)
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

% example of type 3 format data
%[COM1]<BESTXYZ COM1 0 55.5 FINESTEERING 1428 269421.900 00000000 f798 1810
%<     SOL_COMPUTED SINGLE -5046958.9079 2567791.9695 -2925495.1877 6.3757 3.4164 9.4475 SOL_COMPUTED DOPPLER_VELOCITY 0.2202 -0.1010 0.1847 0.7620 0.4083 1.1292 "" 0.150 0.000 0.000 5 4 0 0 0 0 0 0

%% Arguments
switch type
    case 1
        disp('Binary Not Supported');
        data = 0;
    case 2
        disp('Decoding ASCII format');
        data = Decode_BESTXYZ_Ascii(filename);
    case 3
        % run blah
        disp('Decoding Text Format');
        data = Decode_BESTXYZ_Text(filename);
    case 4
        % Read the one with the system timestamp
        disp('Decoding Text Format with Timestamp');
        data = Decode_BESTXYZ_TextWithSystemTimestamp(filename);
    otherwise
        disp(sprintf('Type %d Not Supported',type));
        data = 0;
end


%% Ascii Data Files (commaa delimited Type 2)


function ascii_data = Decode_BESTXYZ_Ascii(filename)

% #BESTXYZA,COM2,0,65.5,FINESTEERING,1489,538924.000,00000000,f798,2580;SOL_COMPUTED,SINGLE,
% -5029388.2086,2551115.1555,-2969796.3019,3.3845,1.7086,3.2199,SOL_COMPUTED,DOPPLER_VELOCITY,
% 0.0153,-0.0058,-0.0406,0.4045,0.2042,0.3849,"",0.150,0.000,0.000,9,9,0,0,0,0,0,0*2be300cb


inputfile = fopen(filename);


if inputfile < 0
    return ;
end

record = 0;



% preallocate matrices for speed
allocated_records_chunk = 5000;
allocated_records = allocated_records_chunk;
ascii_data.PosECEF = zeros(allocated_records,3);
ascii_data.PosECEF_Sigma = zeros(allocated_records,3);
ascii_data.VelECEF = zeros(allocated_records,3);
ascii_data.VelECEF_Sigma = zeros(allocated_records,3);
ascii_data.VelocityLatency = zeros(allocated_records,1);

ascii_data.DifferentialAge = zeros(allocated_records,1);
ascii_data.SolutionAge = zeros(allocated_records,1); 
ascii_data.SatellitesTracked  = zeros(allocated_records,1); 
ascii_data.SatellitesUsed_L1 = zeros(allocated_records,1); 
ascii_data.RTK_SatellitesUsed_L1 = zeros(allocated_records,1);
ascii_data.RTK_SatellitesUsed_L2 = zeros(allocated_records,1);

ascii_data.GPSTime_Week = zeros(allocated_records,1);
ascii_data.GPSTime_Sec = zeros(allocated_records,1);

while ~feof(inputfile)
   %% get 2 lines at a time
   
   line1 = fgetl(inputfile);
   
   
   
   if strfind(line1,'#BESTXYZA')
       record = record + 1;
       
       if record > allocated_records
           
            % allocate another chunk
            disp('Allocating another chunk of records');
            allocated_records = allocated_records+allocated_records_chunk;
            ascii_data.PosECEF(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,3);
            ascii_data.PosECEF_Sigma(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,3);
            ascii_data.VelECEF(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,3);
            ascii_data.VelECEF_Sigma(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,3);
            ascii_data.VelocityLatency(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
          
            ascii_data.DifferentialAge(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            ascii_data.SolutionAge(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            ascii_data.SatellitesTracked(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            ascii_data.SatellitesUsed_L1(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            ascii_data.RTK_SatellitesUsed_L1(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            ascii_data.RTK_SatellitesUsed_L2(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            
            ascii_data.GPSTime_Week(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            ascii_data.GPSTime_Sec(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            
       end
       
       if length(line1) < 68
           break; 
       end
       
        % #BESTXYZA,COM2,0,65.5,FINESTEERING,1489,538924.000,00000000,f798,2580;SOL_COMPUTED,SINGLE,
        % -5029388.2086,2551115.1555,-2969796.3019,3.3845,1.7086,3.2199,SOL_COMPUTED,DOPPLER_VELOCITY,
        % 0.0153,-0.0058,-0.0406,0.4045,0.2042,0.3849,"",0.150,0.000,0.000,9,9,0,0,
        % 0,0,0,0*2be300cb

       
       linedata = textscan(line1,'%s%s%d%f%s%d%f%d%s%s%s%f%f%f%f%f%f%s%s%f%f%f%f%f%f%s%f%f%f%d%d%d%d%d%d%d%s','delimiter',',');
   
   
       ascii_data.GPSTime_Week(record) = linedata{6};
       ascii_data.GPSTime_Sec(record) = linedata{7};
       

       
       try 
       ascii_data.PosECEF(record,1) = linedata{12};
       ascii_data.PosECEF(record,2) = linedata{13};
       ascii_data.PosECEF(record,3) = linedata{14};
       
       ascii_data.PosECEF_Sigma(record,1) = linedata{15};
       ascii_data.PosECEF_Sigma(record,2) = linedata{16};
       ascii_data.PosECEF_Sigma(record,3) = linedata{17};
       
       
       ascii_data.VelECEF(record,1) = linedata{20};
       ascii_data.VelECEF(record,2) = linedata{21};
       ascii_data.VelECEF(record,3) = linedata{22};
       
       ascii_data.VelECEF_Sigma(record,1) = linedata{23};
       ascii_data.VelECEF_Sigma(record,2) = linedata{24};
       ascii_data.VelECEF_Sigma(record,3) = linedata{25};
       
       ascii_data.VelocityLatency(record,1) = linedata{27};
       ascii_data.DifferentialAge(record,1) = linedata{28};
       ascii_data.SolutionAge(record,1) = linedata{29};
       
       ascii_data.SatellitesTracked(record,1) = linedata{30};
       ascii_data.SatellitesUsed_L1(record,1) = linedata{31};
       ascii_data.RTK_SatellitesUsed_L1(record,1) = linedata{32};
       ascii_data.RTK_SatellitesUsed_L2(record,1) = linedata{33};
       
       
       catch
          disp('Warning: failed to allcoate line (may be incomplete)');
          
       end
       
       if mod(record,1000) == 0
           disp(sprintf('Decoded record #%d',record));
       end
       
       
       
   else
       
      continue; 
   end
   

   
end

% trim results
ascii_data.PosECEF = ascii_data.PosECEF(1:record,:);
ascii_data.PosECEF_Sigma = ascii_data.PosECEF_Sigma(1:record,:);
ascii_data.VelECEF = ascii_data.VelECEF(1:record,:);
ascii_data.VelECEF_Sigma = ascii_data.VelECEF_Sigma(1:record,:);
ascii_data.VelocityLatency = ascii_data.VelocityLatency(1:record);
ascii_data.DifferentialAge = ascii_data.DifferentialAge(1:record);
ascii_data.SolutionAge = ascii_data.SolutionAge(1:record);
ascii_data.SatellitesTracked = ascii_data.SatellitesTracked(1:record);
ascii_data.SatellitesUsed_L1 = ascii_data.SatellitesUsed_L1(1:record);
ascii_data.RTK_SatellitesUsed_L1 = ascii_data.RTK_SatellitesUsed_L1(1:record);
ascii_data.RTK_SatellitesUsed_L2 = ascii_data.RTK_SatellitesUsed_L2(1:record);

ascii_data.GPSTime_Week = ascii_data.GPSTime_Week(1:record);
ascii_data.GPSTime_Sec = ascii_data.GPSTime_Sec(1:record);

fclose(inputfile);



%% Text Data Files (Type 3)
function text_data = Decode_BESTXYZ_Text(filename)



inputfile = fopen(filename);
record = 0;


% preallocate matrices for speed
allocated_records_chunk = 5000;
allocated_records = allocated_records_chunk;
text_data.PosECEF = zeros(allocated_records,3);
text_data.PosECEF_Sigma = zeros(allocated_records,3);
text_data.VelECEF = zeros(allocated_records,3);
text_data.VelECEF_Sigma = zeros(allocated_records,3);
text_data.VelocityLatency = zeros(allocated_records,1);

text_data.DifferentialAge = zeros(allocated_records,1);
text_data.SolutionAge = zeros(allocated_records,1); 
text_data.SatellitesTracked  = zeros(allocated_records,1); 
text_data.SatellitesUsed_L1 = zeros(allocated_records,1); 
text_data.RTK_SatellitesUsed_L1 = zeros(allocated_records,1);
text_data.RTK_SatellitesUsed_L2 = zeros(allocated_records,1);



while ~feof(inputfile)
   %% get 2 lines at a time
   
   line1 = fgetl(inputfile);
   
   
   
   if strfind(line1,'BESTXYZ')
       record = record + 1;
       
       if record > allocated_records
           
            % allocate another chunk
            disp('Allocating another chunk of records');
            allocated_records = allocated_records+allocated_records_chunk;
            text_data.PosECEF(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,3);
            text_data.PosECEF_Sigma(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,3);
            text_data.VelECEF(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,3);
            text_data.VelECEF_Sigma(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,3);
            text_data.VelocityLatency(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
          
            text_data.DifferentialAge(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            text_data.SolutionAge(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            text_data.SatellitesTracked(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            text_data.SatellitesUsed_L1(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            text_data.RTK_SatellitesUsed_L1(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
            text_data.RTK_SatellitesUsed_L2(record:record+allocated_records_chunk-1,:) = zeros(allocated_records_chunk,1);
       end
       
       if length(line1) < 68
           break; 
       end
       
       % check for EOF before reding second line...may not be there at all!
       if ~feof(inputfile)
           line2 = fgetl(inputfile);
       else
           break;
       end
       
       linedata = textscan(line1,'%s%s%d%f%s%d%f%d%s%s','delimiter',' ');
   
   
       text_data.GPSTime_Week(record) = linedata{6};
       text_data.GPSTime_Sec(record) = linedata{7};
       
% <     SOL_COMPUTED SINGLE -5046958.9224 2567791.9820 -2925495.2079 6.3756 3.4163 9.4476 
% SOL_COMPUTED DOPPLER_VELOCITY 0.0044 -0.0088 0.0278 0.7620 0.4083 1.1292 "" 0.150 0.000 0.000 5 4 0 0 0 0 0 0

       linedata = textscan(line2,'%s%s%s%s%s%s%s%f%f%f%f%f%f%s%s%f%f%f%f%f%f%s%f%f%f%d%d%d%d%d%d%d%d  ','delimiter',' ');
          
       
       try 
       text_data.PosECEF(record,1) = linedata{8};
       text_data.PosECEF(record,2) = linedata{9};
       text_data.PosECEF(record,3) = linedata{10};
       
       text_data.PosECEF_Sigma(record,1) = linedata{11};
       text_data.PosECEF_Sigma(record,2) = linedata{12};
       text_data.PosECEF_Sigma(record,3) = linedata{13};
       
       
       text_data.VelECEF(record,1) = linedata{16};
       text_data.VelECEF(record,2) = linedata{17};
       text_data.VelECEF(record,3) = linedata{18};
       
       text_data.VelECEF_Sigma(record,1) = linedata{19};
       text_data.VelECEF_Sigma(record,2) = linedata{20};
       text_data.VelECEF_Sigma(record,3) = linedata{21};
       
       text_data.VelocityLatency(record,1) = linedata{23};
       text_data.DifferentialAge(record,1) = linedata{24};
       text_data.SolutionAge(record,1) = linedata{25};
       
       text_data.SatellitesTracked(record,1) = linedata{26};
       text_data.SatellitesUsed_L1(record,1) = linedata{27};
       text_data.RTK_SatellitesUsed_L1(record,1) = linedata{28};
       text_data.RTK_SatellitesUsed_L2(record,1) = linedata{29};
       
       
       catch
          disp('Warning: failed to allcoate line (may be incomplete)');
          
       end
       
       if mod(record,1000) == 0
           disp(sprintf('Decoded record #%d',record));
       end
       
       
       
   else
       
      continue; 
   end
   

   
end

% trim results
text_data.PosECEF = text_data.PosECEF(1:record,:);
text_data.PosECEF_Sigma = text_data.PosECEF_Sigma(1:record,:);
text_data.VelECEF = text_data.VelECEF(1:record,:);
text_data.VelECEF_Sigma = text_data.VelECEF_Sigma(1:record,:);
text_data.VelocityLatency = text_data.VelocityLatency(1:record);

fclose(inputfile);

%% Wrapper function to read the ones with system timestamps and place it in the data format expected by the function
function text_data = Decode_BESTXYZ_TextWithSystemTimestamp(filename)

% Read the log file using Damien's log file reader 
text_data = NMEA_ReadLogFile(filename, 'BESTXYZA', 'new_reader');

% % Preallocate the data array to return
% allocated_records = length(SentenceData);
% text_data.PosECEF = zeros(allocated_records,3);
% text_data.PosECEF_Sigma = zeros(allocated_records,3);
% text_data.VelECEF = zeros(allocated_records,3);
% text_data.VelECEF_Sigma = zeros(allocated_records,3);
% text_data.VelocityLatency = zeros(allocated_records,1);
% 
% text_data.DifferentialAge = zeros(allocated_records,1);
% text_data.SolutionAge = zeros(allocated_records,1); 
% text_data.SatellitesTracked  = zeros(allocated_records,1); 
% text_data.SatellitesUsed_L1 = zeros(allocated_records,1); 
% text_data.RTK_SatellitesUsed_L1 = zeros(allocated_records,1);
% text_data.RTK_SatellitesUsed_L2 = zeros(allocated_records,1);
% 
% % Translate each of the fields in the structure
% for i = 1:allocated_records
% 
%     % SentenceData = NMEA_ReadOEM4Header;
%     % SentenceData.Sentence = '#BESTXYZA';
%     % SentenceData.PositionSolutionStatus = 'POSITIONSOLUTIONSTATUS';
%     % SentenceData.PositionType = 'POSITIONTYPE';
%     text_data.PosECEF(i,1) = SentenceData(i).XCoordinate;
%     text_data.PosECEF(i,2) = SentenceData(i).YCoordinate;
%     text_data.PosECEF(i,3) = SentenceData(i).ZCoordinate;
%     text_data.PosECEF_Sigma(i,1) = SentenceData(i).XPositionStd;
%     text_data.PosECEF_Sigma(i,2) = SentenceData(i).YPositionStd;
%     text_data.PosECEF_Sigma(i,3) = SentenceData(i).ZPositionStd;
%     % SentenceData(i).VelocitySolutionStatus;
%     % SentenceData(i).VelocityType;
%     text_data.VelECEF(i,1) = SentenceData(i).XVelocity;
%     text_data.VelECEF(i,2) = SentenceData(i).YVelocity;
%     text_data.VelECEF(i,3) = SentenceData(i).ZVelocity;
%     text_data.VelECEF_Sigma(i,1) = SentenceData(i).XVelocityStd;
%     text_data.VelECEF_Sigma(i,2) = SentenceData(i).YVelocityStd;
%     text_data.VelECEF_Sigma(i,3) = SentenceData(i).ZVelocityStd;
%     % SentenceData(i).StationID;
%     text_data.VelocityLatency(i) = SentenceData(i).VelocityLatency;
%     text_data.DifferentialAge(i) = SentenceData(i).DifferentialAge;
%     text_data.SolutionAge(i) = SentenceData(i).SolutionAge;
%     text_data.SatellitesTracked(i) = SentenceData(i).NumberObservations;
%     text_data.SatellitesUsed_L1(i) = SentenceData(i).NumberGPSL1RangesUsed;
%     text_data.RTK_SatellitesUsed_L1(i) = SentenceData(i).NumberGPSL1;
%     text_data.RTK_SatellitesUsed_L2(i) = SentenceData(i).NumberGPSL2;
%     % SentenceData(i).Reserved1;
%     % SentenceData(i).Reserved2;
%     % SentenceData(i).Reserved3;
%     % SentenceData(i).Reserved4;
%     % SentenceData(i).Checksum = -1;
%     % SentenceData(i).ChecksumValid = 0;
%     % SentenceData(i).Preamble = 'PREAMBLE';
%     % SentenceData(i).Remainder = 'REMAINDER';
% end
% 



