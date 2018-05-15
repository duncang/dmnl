function [GPSTime_Week, GPSTime_Sec,NumberRinexObsTypes,ValidDataRinexObs,ApproxPos, Novatel_C1, Novatel_L1, Novatel_D1, Novatel_S1, Novatel_P2, Novatel_L2, Novatel_D2, Novatel_S2] = ReadRinexNovatel(Filename)
% Written by Troy Bruggemann (c) CRCSS 2005
% Modified by Duncan Greer to suit Novatel OEM4/V Receiver
% $Id: ReadRinexNovatel.m 3547 2010-06-17 13:27:24Z greerd $
%
%  C1    L1    D1    S1    P2    L2    D2    S2
%
global GPS_PI OMEGAedot mu Earthradius Speedoflight c F L1_f L2_f gamma L1_Wavelength;


% GMearth = 3986005e8; %m^3/s^2
% Earthradius = 6378136;  %m
%    %Earth's rotation rate
% omegadotE = 7.2921151467e-5; %rad/sec
%    %speed of light
% Speedoflight = 2.99792458e8; %m/s
% 
% c = 2.99792458e8;
% L1 = 1575.42e6;
%  
% L1_wavelength = c/L1;


%Filename = 'data/Flight_Data/00282031/00282031.09O';




filename_obsAMicro = fopen(Filename,'r');


%filename_obsAMicro = fopen('data\Ground_Test_Data\2Feb2007\TEST02020330.07O','r');


%% first read the header
for i=1:25


    obslineMicro = fgetl(filename_obsAMicro);
    
    %% look for observation list
    if(strfind(obslineMicro,'TYPES OF OBSERV'))
       NumberRinexObsTypes = str2num(obslineMicro(5:6));
       disp(sprintf('This observation file contains %d observation types',NumberRinexObsTypes));
       
       
       
      for p = 1:NumberRinexObsTypes
       
       ObsType(p,:) = obslineMicro(5+p*6:6+p*6);
   
      end
          
          
       %C1    L1    D1    S1    P2    L2    D2    S2      
       
       
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

%% ----------------------------------------------------
%read in all observation data from observation file
%----------------------------------------------------


i = 1;  % i is the line count, record count for each epoch


%bsline = 1;
%while i < 20,
while 1,   
%while i < 500,
%while i < 2000,
%a(i) = fscanf(filename,'%s',inf);
%[c,d] = fscanf(filename,'%f',inf);

   %read in start of data
   obslineMicro = fgetl(filename_obsAMicro);
   
   if strcmp(obslineMicro,'')
      continue; 
   end
   
   
   if obslineMicro == -1       %EOF reached
      %numberofobs = i;      %POSSIBLE ERROR SOURCE
      break;
   end;

   
   a = sscanf(obslineMicro,'%f');
   
 
     
           
     
     yearAMicro(i) = str2num([obslineMicro(2) obslineMicro(3)]);
     
     monthAMicro(i) = str2num([obslineMicro(5) obslineMicro(6)]);
     
     dayAMicro(i) = str2num([obslineMicro(8) obslineMicro(9)]);
     
     hourAMicro(i) = str2num([obslineMicro(11) obslineMicro(12)]);
     
     minuteAMicro(i) = str2num([obslineMicro(14) obslineMicro(15)]);
     
     secondAMicro(i) = str2num([obslineMicro(17:26)]);
                 
     epochflagAMicro(i) = str2num([obslineMicro(29)]);
          
     numsatsAMicro(i) = str2num([obslineMicro(31) obslineMicro(32)]);
     
     
     civiliandate_vecAMicro(i,1:6) = [yearAMicro(i) monthAMicro(i) dayAMicro(i) hourAMicro(i) minuteAMicro(i) secondAMicro(i)] ;
     
     [GPSTime_Week(i) GPSTime_Sec(i)] = ftime(civiliandate_vecAMicro(i,1:6));
     
     
     
     
     
      
      % Put PRN's in array
      
     % position = findstr(' ',obsline);
       
     % for count = 1:length(position)
      
     % PRNvecA(i,count) = str2num(obsline((position(count)+1):(position(count)+2))); 
      
     % end
     
     
     %position = findstr(' ',obsline);
     
     
     % Put PRN's in array - This is if the rinex file has G's separating
     % the PRN numbers
      %----------------------------------------------
      position = findstr('G',obslineMicro);
      
      Glonass = findstr('R',obslineMicro);
      numberGlonass(i) = length(Glonass);
      
      for countMicro = 1:length(position)
      
      PRNvecAMicro(i,countMicro) = str2num(obslineMicro((position(countMicro)+1):(position(countMicro)+2)));
      
      end
      
      firstcountMicro = length(position);
      
      for countGlonass = 1:numberGlonass(i)
         GlonassVec(i,countGlonass) = str2num(obslineMicro((Glonass(countGlonass)+1):(Glonass(countGlonass)+2)));
      end
      
      firstcountGlonass = numberGlonass(i);
      
      
      %% check if there is another line of satellite PRNs
      if firstcountMicro + firstcountGlonass ~= numsatsAMicro(i)
          obslineMicro = fgetl(filename_obsAMicro);
        
          position = findstr('G',obslineMicro);

          Glonass = findstr('R',obslineMicro);
          numberGlonass(i) = numberGlonass(i) + length(Glonass);

          for countMicro = 1:length(position)

          PRNvecAMicro(i,firstcountMicro+countMicro) = str2num(obslineMicro((position(countMicro)+1):(position(countMicro)+2)));

          end

          for countGlonass = 1:length(Glonass)
             GlonassVec(i,firstcountGlonass+countGlonass) = str2num(obslineMicro((Glonass(countGlonass)+1):(Glonass(countGlonass)+2)));
          end
      
        
        
      end
       %----------------------------------------------
     
    %Use this if the rinex file does NOT have G's separating the PRN
    %numbers
      
%         for countMicro = 1:numsatsAMicro(i)
%         
%                 
%         PRNvecAMicro(i,countMicro)=  str2num(obslineMicro((32+3*(countMicro-1)+1):(32+3*(countMicro-1)+3))); 
%         
%         end
%         
     
        
    
      
          
      
      %rxclockoffset(i) = a(9+numsats(i));
    
      %end         
      
      %Get observation Data
    
      
    for j = 1:numsatsAMicro(i)
       obslineMicro = fgetl(filename_obsAMicro);
       
       
       %extend obslineMicro to full vector length with white space
       
      VecLength =  size(obslineMicro);
           
      
      if VecLength(2) < 79
                 
         BlankVec = blanks(79-VecLength(2));
         obslineMicro = [obslineMicro, BlankVec];
       
      end
                               
%        L1AMicro(j,i) = str2double(obslineMicro(1:14));
%        LLindL1AMicro(j,i) = str2double(obslineMicro(15));        %loss of lock indicator for L1
%        %signalstrengthL1(i,j) = str2double(obsline(16));
% 
% 	    
%        L2AMicro(j,i) = str2double(obslineMicro(17:30));
%        LLindL2AMicro(j,i) = str2double(obslineMicro(31));        %loss of lock indicator for L1
%        
%        
%        
%        C1AMicro(j,i) = str2double(obslineMicro(33:46));
%        C1indL2AMicro(j,i) = str2double(obslineMicro(47));        %loss of lock indicator for L1
%        %signalstrengthL2(i,j) = str2double(obsline(32));
%        
%        P1AMicro(j,i) = str2double(obslineMicro(49:61));
%        C1indL2AMicro(j,i) = str2double(obslineMicro(62));        %loss of lock indicator for L1
%        
%        P2AMicro(j,i) = str2double(obslineMicro(65:77));
%        C1indL2AMicro(j,i) = str2double(obslineMicro(78));        %loss of lock indicator for L1
%        %L1(j) = L1tmp(1:15);
%        %L2(j) = a(2);
%        
%         obslineMicro = fgetl(filename_obsAMicro); % GET THE DOPPLER DATA WHICH IS ON THE NEXT LINE
%        %D1A(j,i) = str2double(obsline(33:46));
%        %D1(i,j) = str2double(obsline(47:62));

        %%C1    L1    D1    S1    P2    L2    D2    S2
%          07 02 02 05 14 10.5000000  0  5G13G28G 8G19G27                                 
%           21373676.73749 117487306.59249     -2523.30949        47.903    21373676.51246
%           87521633.30346     -1966.21946        43.297  
%           20763681.37649 114281762.48249      1031.07449        50.007    20763679.70348
%           89065127.95548       803.43048        46.038  
%           21349011.42149 117357689.08349       767.30149        49.726    21349010.61147
%           91461963.82447       597.89547        44.161  
%           23064934.23347 126374929.36947     -2110.11347        45.491    23064932.43745
%           98488378.11245     -1644.24245        39.778  
%           20592734.89649 113383431.41549      -442.95749        49.399    20592733.64747
%           88365138.03047      -345.16447        44.685  
%       

        C1(j,i) = str2double(obslineMicro(3:16));
        L1(j,i) = str2double(obslineMicro(18:32));
        D1(j,i) = str2double(obslineMicro(34:48));
        S1(j,i) = str2double(obslineMicro(50:62));

        % if were only expecting 4 observersions (L1 only), go around here
        if NumberRinexObsTypes == 8
            % check if we have the other obs
            if strcmp(obslineMicro(64:79),'                ')
                P2(j,i) = 0;
                L2(j,i) = 0;
                D2(j,i) = 0;
                S2(j,i) = 0;
                obslineMicro = fgetl(filename_obsAMicro);
            else
             P2(j,i) = str2double(obslineMicro(64:79));
             obslineMicro = fgetl(filename_obsAMicro);
             if strcmp(obslineMicro,'')
                    L2(j,i) = 0;
                    D2(j,i) = 0;
                    S2(j,i) = 0;

             else
                 L2(j,i) = str2double(obslineMicro(3:16));
                 D2(j,i) = str2double(obslineMicro(18:32));
                 S2(j,i) = str2double(obslineMicro(34:48));
             end
            end
        else
            P2(j,i) = 0;
            L2(j,i) = 0;
            D2(j,i) = 0;
            S2(j,i) = 0;
        end
        

     end 
     
     
     %Order according to PRN
     
      

   

            for prnMicro = 1:32
                 for k = 1:numsatsAMicro(i) - numberGlonass(i)    % first element is just gps time
                           
                
                     if (PRNvecAMicro(i,k) == prnMicro)  
                        
                        %SV_PRN(prn,i) = SVtoplot(i,j);
%                         L1_PRNAMicro(prnMicro,i) = L1AMicro(k,i);   
%                   
%                         L2_PRNAMicro(prnMicro,i) = L2AMicro(k,i);  
%                         
%                         C1_PRNAMicro(prnMicro,i) = C1AMicro(k,i); 
%                                                   
%                   
%                         P1_PRNAMicro(prnMicro,i) = P1AMicro(k,i); 
%                   
%                         P2_PRNAMicro(prnMicro,i) = P2AMicro(k,i); 
                        %D1_PRNA(prn,i) = D1A(k,i); 
                        % C1    L1    D1    S1    P2    L2    D2    S2
                         Novatel_C1(prnMicro,i) = C1(k,i);
                         Novatel_L1(prnMicro,i) = L1(k,i);
                         Novatel_D1(prnMicro,i) = D1(k,i);
                         Novatel_S1(prnMicro,i) = S1(k,i);
                         Novatel_P2(prnMicro,i) = P2(k,i);
                         Novatel_L2(prnMicro,i) = L2(k,i);
                         Novatel_D2(prnMicro,i) = D2(k,i);
                         Novatel_S2(prnMicro,i) = S2(k,i);
%                          Novatel_P2(prnMicro,i) = 0;
%                          Novatel_L2(prnMicro,i) = 0;
%                          Novatel_D2(prnMicro,i) = 0;
%                          Novatel_S2(prnMicro,i) = 0;
                         
                        
                        
                        ValidDataRinexObs(prnMicro,i) = 1;
                                                                            
                    
                     %else
                        %ValidDataRinexObs(prnMicro,i) = 0;
                     end
                 end
            end                           
            
                       
            
    i = i+1;
end   


%extend to 32 elements if necessary (easier to deal with fixed array size
%rather than variable

sizevec = size(Novatel_C1);
shift = 32 - sizevec(1);
catB = zeros(shift,sizevec(2));

Novatel_C1 = cat(1, Novatel_C1, catB);
Novatel_L1 = cat(1, Novatel_L1, catB);
Novatel_D1 = cat(1, Novatel_D1, catB);
Novatel_S1 = cat(1, Novatel_S1, catB);
Novatel_P2 = cat(1, Novatel_P2, catB);
Novatel_L2 = cat(1, Novatel_L2, catB);
Novatel_D2 = cat(1, Novatel_D2, catB);
Novatel_S2 = cat(1, Novatel_S2, catB);

ValidDataRinexObs = cat(1, ValidDataRinexObs, catB);



fclose(filename_obsAMicro);

disp(sprintf('Finished reading %s.  Found %d Observations',Filename,i));

