function [GPStimeAMicroWeek, GPStimeAMicro,ValidDataRinexObs, L1_PRNAMicro, L2_PRNAMicro, C1_PRNAMicro, P1_PRNAMicro, P2_PRNAMicro] = ReadRinexSblockRoof(Filename)


%THIS IS READER FOR ASHTECH ROOF GPS RECEIVER OKAY! ASHTECH MICRO-Z
%RECEIVER in S BLOCK ROOF GROUND STATION QUT

% written by Troy Bruggemann (c) CRCSS 2005
% last modified 25 May 2005

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


filename_obsAMicro = fopen(Filename,'r');


%filename_obsAMicro = fopen('C:\Documents and Settings\n2523710\Desktop\ThesisTests\things\29.1124HrRoof\compareMatlab\EESE3343.04O','r');


%----------------------------------------------------
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
     
     [GPStimeAMicroWeek(i) GPStimeAMicro(i)] = ftime(civiliandate_vecAMicro(i,1:6));
     
     
     
     
     
      
      % Put PRN's in array
      
     % position = findstr(' ',obsline);
       
     % for count = 1:length(position)
      
     % PRNvecA(i,count) = str2num(obsline((position(count)+1):(position(count)+2))); 
      
     % end
     
     
     %position = findstr(' ',obsline);
     
     
     % Put PRN's in array - This is if the rinex file has G's separating
     % the PRN numbers
      %----------------------------------------------
%       position = findstr('G',obslineMicro);
%        
%       for countMicro = 1:length(position)
%       
%       PRNvecAMicro(i,countMicro) = str2num(obslineMicro((position(countMicro)+1):(position(countMicro)+2))); 
%       
%       end
       %----------------------------------------------
     
    %Use this if the rinex file does NOT have G's separating the PRN
    %numbers
      
        for countMicro = 1:numsatsAMicro(i)
        
                
        PRNvecAMicro(i,countMicro)=  str2num(obslineMicro((32+3*(countMicro-1)+1):(32+3*(countMicro-1)+3))); 
        
        end
        
     
        
    
      
          
      
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
                               
       L1AMicro(j,i) = str2double(obslineMicro(1:14));
       LLindL1AMicro(j,i) = str2double(obslineMicro(15));        %loss of lock indicator for L1
       %signalstrengthL1(i,j) = str2double(obsline(16));

	    
       L2AMicro(j,i) = str2double(obslineMicro(17:30));
       LLindL2AMicro(j,i) = str2double(obslineMicro(31));        %loss of lock indicator for L1
       
       
       
       C1AMicro(j,i) = str2double(obslineMicro(33:46));
       C1indL2AMicro(j,i) = str2double(obslineMicro(47));        %loss of lock indicator for L1
       %signalstrengthL2(i,j) = str2double(obsline(32));
       
       P1AMicro(j,i) = str2double(obslineMicro(49:61));
       C1indL2AMicro(j,i) = str2double(obslineMicro(62));        %loss of lock indicator for L1
       
       P2AMicro(j,i) = str2double(obslineMicro(65:77));
       C1indL2AMicro(j,i) = str2double(obslineMicro(78));        %loss of lock indicator for L1
       %L1(j) = L1tmp(1:15);
       %L2(j) = a(2);
       
        obslineMicro = fgetl(filename_obsAMicro); % GET THE DOPPLER DATA WHICH IS ON THE NEXT LINE
       %D1A(j,i) = str2double(obsline(33:46));
       %D1(i,j) = str2double(obsline(47:62));
     end 
     
     
     %Order according to PRN
     
      

   

            for prnMicro = 1:32
                 for k = 1:numsatsAMicro(i)    % first element is just gps time
                           
                
                     if (PRNvecAMicro(i,k) == prnMicro)  
                        
                        %SV_PRN(prn,i) = SVtoplot(i,j);
                        L1_PRNAMicro(prnMicro,i) = L1AMicro(k,i);   
                  
                        L2_PRNAMicro(prnMicro,i) = L2AMicro(k,i);  
                        
                        C1_PRNAMicro(prnMicro,i) = C1AMicro(k,i); 
                                                  
                  
                        P1_PRNAMicro(prnMicro,i) = P1AMicro(k,i); 
                  
                        P2_PRNAMicro(prnMicro,i) = P2AMicro(k,i); 
                        %D1_PRNA(prn,i) = D1A(k,i); 
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

sizevec = size(L1_PRNAMicro);
shift = 32 - sizevec(1);
catB = zeros(shift,sizevec(2));

L1_PRNAMicro = cat(1, L1_PRNAMicro, catB);
L2_PRNAMicro = cat(1, L2_PRNAMicro, catB);
C1_PRNAMicro = cat(1, C1_PRNAMicro, catB);
P1_PRNAMicro = cat(1, P1_PRNAMicro, catB);
P2_PRNAMicro = cat(1, P2_PRNAMicro, catB);
ValidDataRinexObs = cat(1, ValidDataRinexObs, catB);



fclose(filename_obsAMicro);