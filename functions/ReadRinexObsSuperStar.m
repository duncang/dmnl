function [GPStimeAMicroWeek, GPStimeAMicro,ValidDataRinexObs, L1_PRNAsh, C1_PRNAsh] = ReadRinexObsSuperStar(Filename)


% By Troy Bruggemann  2 June 2004
%This Reads in observation data from a Rinex 2.0 observation file, NOT
%including the header, so delete the header section off the rinex
%observation file to use this.
 

%clear all





%filename_obs = fopen('C:\Documents and Settings\n2523710\Desktop\Roc6ppsZFinalVersionForCP\25.5compareWithAshtekTests\O.txt','r');

%filename_obs = fopen('C:\Documents and Settings\n2523710\Desktop\Documentation\CarrierPhase\TestData\Roc6ppsZFinalVersionForCP\25.5compareWithAshtekTests\2.6CompareReceivers\O.txt','r');

% filename_obsAsh = fopen('D:\Documents and Settings\n2523710\Desktop\SuperStar\formatlabclock_not_compensated.05o','r');


filename_obsAsh = fopen(Filename,'r');

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
   obsline = fgetl(filename_obsAsh);
   
   
   
   if obsline == -1       %EOF reached
      %numberofobs = i;      %POSSIBLE ERROR SOURCE
      break;
   end;

   
   a = sscanf(obsline,'%f');
   
 
     
           
     
     year(i) = str2num([obsline(2) obsline(3)]);
     
     month(i) = str2num([obsline(5) obsline(6)]);
     
     day(i) = str2num([obsline(8) obsline(9)]);
     
     hour(i) = str2num([obsline(11) obsline(12)]);
     
     minute(i) = str2num([obsline(14) obsline(15)]);
     
     second(i) = str2num([obsline(17:26)]);
                 
     epochflag(i) = str2num([obsline(29)]);
     
     numsats(i) = str2num([obsline(31) obsline(32)]);
         
       
     civiliandate_vecAMicro(i,1:6) = [year(i) month(i) day(i) hour(i) minute(i) second(i)] ;
     
               
     [GPStimeAMicroWeek(i) GPStimeAMicro(i)] = ftime(civiliandate_vecAMicro(i,1:6));
     
      
      % Put PRN's in array
      
      position = findstr('G',obsline);
       
      for count = 1:length(position)
      
      PRNvec(i,count) = str2num(obsline((position(count)+1):(position(count)+2))); 
      
      end
      
          
      
      %rxclockoffset(i) = a(9+numsats(i));
    
      %end         
      
      %Get observation Data
    
    for j = 1:numsats(i)
       obsline = fgetl(filename_obsAsh);
                               
       C1Ash(j,i) = str2double(obsline(1:14));
       C1indL1(j,i) = str2double(obsline(15));        %loss of lock indicator for L1
       %signalstrengthL1(i,j) = str2double(obsline(16));

		 
       L1Ash(j,i) = str2double(obsline(17:30));
       L1indL2(j,i) = str2double(obsline(31));        %loss of lock indicator for L1
       %signalstrengthL2(i,j) = str2double(obsline(32));
       
       %L1(j) = L1tmp(1:15);
       %L2(j) = a(2);
       %D1Ash(j,i) = str2double(obsline(33:46));
       %D1(i,j) = str2double(obsline(47:62));
     end 
     
     
     %Order according to PRN


            for prn = 1:32
                 for k = 1:numsats(i)    % first element is just gps time
                           
                
                     if (PRNvec(i,k) == prn)  
                        
                        %SV_PRN(prn,i) = SVtoplot(i,j);
                        C1_PRNAsh(prn,i) = C1Ash(k,i);   %REVERSED THESE FOR REPEATER TEST, ASHTECH PUTS L1 FIRST NOT C1
                  
                        L1_PRNAsh(prn,i) = L1Ash(k,i); 
                  
                       % D1_PRNAsh(prn,i) = D1Ash(k,i); 
                        ValidDataRinexObs(prn,i) = 1;
                                                            
                    end    
                        
                 end
            end
       
    i = i+1;
end   

%extend to 32 elements if necessary (easier to deal with fixed array size
%rather than variable

sizevec = size(L1_PRNAsh);
shift = 32 - sizevec(1);
catB = zeros(shift,sizevec(2));

L1_PRNAsh = cat(1, L1_PRNAsh, catB);
C1_PRNAsh = cat(1, C1_PRNAsh, catB);
ValidDataRinexObs = cat(1, ValidDataRinexObs, catB);

fclose(filename_obsAsh);


