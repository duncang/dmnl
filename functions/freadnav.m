function [nav_vec, ALPHA, BETA] = freadnav(Filename)
% [nav_vec, ALPHA, BETA] = freadnav(Filename)
%By Troy Bruggemann
% $Id: freadnav.m 1874 2008-07-15 04:42:16Z n2523710 $
%--------------------------------------------
%read in navigation data from navigation file
%--------------------------------------------
%     1: PRN_nav_vec(i) 
%     2: toc_nav_sec 
%     3: svclkbias_nav 
%     4: svclkdrift_nav 
%     5: svclkdriftrate_nav 
%     6: IODE_nav
%     7: crs_nav 
%     8: deltan_nav 
%     9: M0_nav 
%     10: cuc_nav 
%     11: e_nav 
%     12: cus_nav 
%     13: sqrta_nav 
%     14: toe_nav 
%     15: cic_nav 
%     16: OMEGA0_nav
%     17: cis_nav 
%     18: i0_nav 
%     19: crc_nav 
%     20: omega_nav 
%     21: OMEGA_DOT_nav 
%     22: idot_nav 
%     23: CodeL2_nav 
%     24: GPSwknum_nav 
%     25: L2Pdataflag_nav 
%     26: SVaccuracy_nav 
%     27: SVhealth_nav 
%     28: tgd_nav 
%     29: IODC_nav
%     30: txtimeofmessage_nav 
%     31: fitinterval_nav
%     32:  toc_nav_week
%
%
%

disp('[freadnav] Did you remember to replace Ds with Es?');

filename_nav = fopen(Filename,'r');

%% read nav header
for i=1:8
    obsline_nav = fgetl(filename_nav);

    
    %% find the iono parameters
    if(strfind(obsline_nav,'ION ALPHA'))
        ALPHA(1) = str2double(obsline_nav(5:14));
        ALPHA(2) = str2double(obsline_nav(16:26));
        ALPHA(3) = str2double(obsline_nav(28:38));
        ALPHA(4) = str2double(obsline_nav(40:50));
    end
    

    if(strfind(obsline_nav,'ION BETA'))
        BETA(1) = str2double(obsline_nav(5:14));
        BETA(2) = str2double(obsline_nav(16:26));
        BETA(3) = str2double(obsline_nav(28:38));
        BETA(4) = str2double(obsline_nav(40:50));
    end
    
    if(strfind(obsline_nav,'END OF HEADER'))
        disp(sprintf('Found EOH at line %d',i));
        break;
    end
end

if ~exist('ALPHA','var')
    ALPHA = zeros(1,4);
    disp('No Ion ALPHA data');
end

if ~exist('BETA','var');
    BETA = zeros(1,4);
    disp('No Ion BETA data');
end

%% read remainder of file


i = 1;
k = 0;
j = 1;
%PRN_nav_vec = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %initialise

%PRN_nav_vec = zeros(100,1); %initialise  
while 1,
%while i < numberofobs,
%for j = 1:numsats,
   repeatedflag = 0;
   %read in start of data
   obsline_nav = fgetl(filename_nav);
   if obsline_nav == -1 
      %numberofobs = i;      %POSSIBLE ERROR SOURCE
      break;
   end;

   
   a_nav = sscanf(obsline_nav,'%f');
   
   %PRN_nav_vec(j) = a_nav(1);
   PRN_nav = a_nav(1);
   
   %for count = 1:i,
      %if PRN_nav == PRN_nav_vec(count)   %check for repeated satellite data
       %  repeatedflag = 1;     %satellite already read in from nav file
      %   break;
      %end  
   %end     
          
   if repeatedflag ~= 1
      PRN_nav_vec(i) = PRN_nav;  %vector for holding satellite list
   end         

   toc_nav = [a_nav(2),a_nav(3),a_nav(4),a_nav(5),a_nav(6),a_nav(7)];
   svclkbias_nav = a_nav(8);
   svclkdrift_nav = a_nav(9);
   svclkdriftrate_nav = a_nav(10);
   
   %next line 
   obsline_nav = fgetl(filename_nav);
   a_nav = sscanf(obsline_nav,'%f');
   
   
   %data
   IODE_nav = a_nav(1);
   crs_nav = a_nav(2);
   deltan_nav = a_nav(3);
   M0_nav = a_nav(4);
   
   %next line 
   obsline_nav = fgetl(filename_nav);
   a_nav = sscanf(obsline_nav,'%f');
   
   
   cuc_nav = a_nav(1);
   e_nav = a_nav(2);
   cus_nav = a_nav(3);
   sqrta_nav = a_nav(4);
   %next line 
   obsline_nav = fgetl(filename_nav);
   a_nav = sscanf(obsline_nav,'%f');
   
   
   toe_nav = a_nav(1);   %sec of GPS week
   cic_nav = a_nav(2);
   OMEGA0_nav = a_nav(3);
   cis_nav = a_nav(4);
   %next line 
   obsline_nav = fgetl(filename_nav);
   a_nav = sscanf(obsline_nav,'%f');
   
   
   i0_nav = a_nav(1);
   crc_nav = a_nav(2);
   omega_nav = a_nav(3);
   OMEGA_DOT_nav = a_nav(4);
   %next line 
   obsline_nav = fgetl(filename_nav);
   a_nav = sscanf(obsline_nav,'%f');
   
   
   idot_nav = a_nav(1);
   CodeL2_nav = a_nav(2);
   GPSwknum_nav = a_nav(3);
   L2Pdataflag_nav = a_nav(4);
   %next line 
   obsline_nav = fgetl(filename_nav);
   a_nav = sscanf(obsline_nav,'%f');
   
   
   SVaccuracy_nav = a_nav(1);
   SVhealth_nav = a_nav(2);
   tgd_nav = a_nav(3);
   IODC_nav = a_nav(4);
   
   %next line 
   obsline_nav = fgetl(filename_nav);
   a_nav = sscanf(obsline_nav,'%f');
   
   txtimeofmessage_nav = a_nav(1);
   fitinterval_nav = a_nav(2);
   
   
         
 
 	[toc_nav_week toc_nav_sec] = ftime(toc_nav); %toc in GPS week and GPS seconds, input in epoch vector format
 
	
	if repeatedflag ~= 1

		nav_vec(i,1:32) =  [PRN_nav_vec(i) toc_nav_sec svclkbias_nav svclkdrift_nav svclkdriftrate_nav IODE_nav ...
      						crs_nav deltan_nav M0_nav cuc_nav e_nav cus_nav sqrta_nav toe_nav cic_nav OMEGA0_nav ...
                  		cis_nav i0_nav crc_nav omega_nav OMEGA_DOT_nav idot_nav CodeL2_nav GPSwknum_nav ...
                  		L2Pdataflag_nav SVaccuracy_nav SVhealth_nav tgd_nav IODC_nav ...
      						txtimeofmessage_nav fitinterval_nav toc_nav_week];  %note added toc_nav_week 6/9/05 TB.
      
      
    %if no repeats then continue   
    	i = i + 1; 
	end
   
     
      
   end
   
   fclose(filename_nav);
   
if ~exist('ALPHA','var')
   ALPHA = zeros(1,4); 
   disp('Warning: ION ALPHA values not found!');
end
   
if ~exist('BETA','var')
   BETA = zeros(1,4); 
   disp('Warning: ION BETA values not found!');
end
 
   
   
