function [GPSTime_Week, GPSTime_Sec,NumberRinexObsTypes,ApproxPos, DATA_STRUCT] = ReadRinexGRS(Filename)
%
%  [GPSTime_Week, GPSTime_Sec,NumberRinexObsTypes,ApproxPos, DATA_STRUCT] = ReadRinexGRS(Filename)
%
% Written by Troy Bruggemann (c) CRCSS 2005
% Modified by Duncan Greer to suit Novatel OEM4/V Receiver
% Modified again by TB
% $Id: GPSandModelEKFwithAeroModel.m,v 1 2007/02/05 06:13:24 n2523710 Exp
%

%read the DATA_STRUCT by
%DATASTRUCT.C1(PRN,:); etc.

feature accel on



filename_obsAMicro = fopen(Filename,'r');

%% first read the header
for i=1:24
    obslineMicro = fgetl(filename_obsAMicro);
    %% look for observation list
    if(strfind(obslineMicro,'TYPES OF OBSERV'))
        NumberRinexObsTypes = str2num(obslineMicro(5:6));
        disp(sprintf('This observation file contains %d observation types',NumberRinexObsTypes));

        for p = 1:NumberRinexObsTypes
            ObsType(p,:) = obslineMicro(5+p*6:6+p*6);
        end
    end

    %% read the approximate position
    if(strfind(obslineMicro,'APPROX POSITION XYZ'))
        ApproxPos(1) = str2num(obslineMicro(1:13));
        ApproxPos(2) = str2num(obslineMicro(15:26));
        ApproxPos(3) = str2num(obslineMicro(29:41));
        [ApproxPos_LLH(1),ApproxPos_LLH(2),ApproxPos_LLH(3)] = ECEF2LLH(ApproxPos);
        disp(sprintf('Approx Starting Position: %fN %fE %fm',ApproxPos_LLH(1)*180/pi,ApproxPos_LLH(2)*180/pi,ApproxPos_LLH(3)));
    end

    if(strfind(obslineMicro,'END OF HEADER'))
        %disp(sprintf('Found EOH at line %d',i));
        break;
    end
end


%
%% ----------------------------------------------------
%read in all observation data from observation file
%----------------------------------------------------

i = 1;  % i is the line count, record count for each epoch

while 1,

    %read in start of data
    obslineMicro = fgetl(filename_obsAMicro);

    if obslineMicro == -1       %EOF reached
        %numberofobs = i;      %POSSIBLE ERROR SOURCE
        break;
    end;


     sscanf(obslineMicro,'%f');

    yearAMicro(i) = str2num([obslineMicro(2) obslineMicro(3)]);
    monthAMicro(i) = str2num([obslineMicro(5) obslineMicro(6)]);
    dayAMicro(i) = str2num([obslineMicro(8) obslineMicro(9)]);
    hourAMicro(i) = str2num([obslineMicro(11) obslineMicro(12)]);
    minuteAMicro(i) = str2num([obslineMicro(14) obslineMicro(15)]);
    secondAMicro(i) = str2num([obslineMicro(17:26)]);
    epochflagAMicro(i) = str2num(obslineMicro(29));
    numsatsAMicro(i) = str2num([obslineMicro(31) obslineMicro(32)]);
    civiliandate_vecAMicro(i,1:6) = [yearAMicro(i) monthAMicro(i) dayAMicro(i) hourAMicro(i) minuteAMicro(i) secondAMicro(i)];
    [GPSTime_Week(i) GPSTime_Sec(i)] = ftime(civiliandate_vecAMicro(i,1:6));
    
    % Put PRN's in array - This is if the rinex file has G's separating
    % the PRN numbers
    %----------------------------------------------
    position = findstr('G',obslineMicro);

    for countMicro = 1:length(position)
        PRNvecAMicro(i,countMicro) = str2num(obslineMicro((position(countMicro)+1):(position(countMicro)+2)));
    end
    
    if length(position) < numsatsAMicro(i)
        % there is another line of obs header
        obslineMicro = fgetl(filename_obsAMicro);
        
        position = findstr('G',obslineMicro);

        for countMicro2 = 1:length(position)
            PRNvecAMicro(i,countMicro+countMicro2) = str2num(obslineMicro((position(countMicro2)+1):(position(countMicro2)+2)));
        end
    end
    
    %----------------------------------------------

    %Use this if the rinex file does NOT have G's separating the PRN
    %numbers
%     for countMicro = 1:numsatsAMicro(i)
%         PRNvecAMicro(i,countMicro)=  str2num(obslineMicro((32+3*(countMicro-1)+1):(32+3*(countMicro-1)+3)));
%     end


    %Get observation Data


    for j = 1:numsatsAMicro(i)   
        if NumberRinexObsTypes >= 5
            numobs = 5;
        else
            numobs = NumberRinexObsTypes;
        end
        
        nchar = 84;  %number of characters to read
        obslineMicro = fgets(filename_obsAMicro,nchar);
        
        %search for end of line (end of line might terminate early for
        %a line with missing records
        
        endofline = length(obslineMicro);
        
        %this endofline position tells the next piece of code how far to
        %read into obslineMicro
        
              
        %if ~ischar(obslineMicro)
        for p = 1:numobs
            
            readstart = (p-1)*16+1;
            readend = (p-1)*16+14;
            
            %this stops it trying to read past the end of obslineMicro
            if readend > endofline
                readend = endofline;
            end
            
            OBS_VEC(j,p,i,:) = str2double(obslineMicro(readstart:readend));
            
             if OBS_VEC(j,p,i,:) == NaN
                    OBS_VEC(j,p,i,:) = 0;
                end
            
        end
        %if second line exists, read it in
        if NumberRinexObsTypes >5
            nchar = 84;  %number of characters to read
            %get the next line
            obslineMicro = fgets(filename_obsAMicro,nchar);
            endofline = length(obslineMicro);
            
            for k = 1:NumberRinexObsTypes-5
                
            readstart = (k-1)*16+1;
            readend = (k-1)*16+14;
            
            %this stops it trying to read past the end of obslineMicro
            if readend > endofline
                readend = endofline;
            end
                
                OBS_VEC(j,k+5,i,:) = str2double(obslineMicro(readstart:readend));
                
                if OBS_VEC(j,k+5,i,:) == NaN
                    OBS_VEC(j,k+5,i,:) = 0;
                end
                
            end
        end
    end
        
         
    
   %order according to PRN                     

   j = 1;  
   for countMicro = 1:numsatsAMicro(i)
              
       OBS_VEC_PRN(PRNvecAMicro(i,countMicro),:,i,:) = OBS_VEC(j,:,i,:);   
       j = j+1;
       
   end%     
    
    
    i = i+1;



    disp(sprintf('Processed %d observations',i));
    
end


 for q = 1:NumberRinexObsTypes
     
     expression = ObsType(q,1:2);
 
 DATA_STRUCT.(expression) = [(OBS_VEC_PRN(:,:,:,q))];

 end
 

fclose(filename_obsAMicro);

disp(sprintf('Finished reading %s.  Found %d Observations',Filename,i));

