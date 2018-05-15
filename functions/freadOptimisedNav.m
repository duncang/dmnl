function [nav_vec, TruePosVec] = freadOptimisedNav(Filename)

%generates optimised 24 GPS satellite constellation data from a file which
%is Table B-2 of Appendix B in DO229C WAAS MOPS.

filename_nav = fopen(Filename,'r');

   


i = 1;
k = 0;
j = 1;
readheaderflag = 1; %to make sure read header only once. header contains semi major, eccentricity, and orbital plane inclination. 

%PRN_nav_vec = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %initialise

PRN_nav_vec = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; %initialise  % for 31 satellites
while 1,
%while i < numberofobs,
%for j = 1:numsats,



if readheaderflag ==1 
   %read in start of data
   obsline_nav = fgetl(filename_nav);
   if obsline_nav == -1 
      %numberofobs = i;      
      break;
   end;
   
   
   readheaderflag = 0; 
   
   a_nav = sscanf(obsline_nav,'%f');  

   
   toc_nav = [a_nav(1),a_nav(2),a_nav(3),a_nav(4),a_nav(5),a_nav(6)];
   toc_nav_week = a_nav(7);
   toc_nav_sec = a_nav(8);   

   %these will be constants over time in the optimised thingy	
   semimajor = a_nav(9);
   eccentricity = a_nav(10);
   ii =  a_nav(11);
   inclination = ii*pi/180;   
   sqrt_a = sqrt(semimajor);
   
   
else  %read the rest of the data
    
 obsline_nav = fgetl(filename_nav);
   if obsline_nav == -1 
      %numberofobs = i;      %POSSIBLE ERROR SOURCE
      break;
   end;

   
     a_nav = sscanf(obsline_nav,'%f');      
     


   PRN_nav_vec(i) = a_nav(1);  %vector for holding satellite list


 
   svclkbias_nav = 0;
   svclkdrift_nav = 0;
   svclkdriftrate_nav = 0;   
   
   
   %data
   IODE_nav = 0;
   crs_nav = 0;
   deltan_nav = 0;
   monavtemp = a_nav(2);
   M0_nav = monavtemp*pi/180; %done
   
    
   cuc_nav = 0;
   e_nav = eccentricity;
   cus_nav = 0;
   sqrta_nav = sqrt_a;
     
   
   toe_nav = toc_nav_sec;   %sec of GPS week
   cic_nav =0;
   onavtemp = a_nav(3);
   OMEGA0_nav = onavtemp*pi/180;
   cis_nav = 0;
    
   
   i0_nav = inclination;
   crc_nav = 0;
   omega_nav = 0;
   OMEGA_DOT_nav = 0;
 
   
   idot_nav = 0;
   CodeL2_nav = 0;
   GPSwknum_nav = toc_nav_week; %made the toe the same as toc
   L2Pdataflag_nav = 0;
 
   
   
   SVaccuracy_nav = 0;
   SVhealth_nav = 0;
   tgd_nav = 0;
   IODC_nav = 0;
  
   
   txtimeofmessage_nav = 0;
   fitinterval_nav = 0;
   
   
   
         
 
 	%[toc_nav_week toc_nav_sec] = ftime(toc_nav); %toc in GPS week and GPS seconds, input in epoch vector format
 
   
	

		nav_vec(i,1:32) =  [PRN_nav_vec(i) toc_nav_sec svclkbias_nav svclkdrift_nav svclkdriftrate_nav IODE_nav ...
      						crs_nav deltan_nav M0_nav cuc_nav e_nav cus_nav sqrta_nav toe_nav cic_nav OMEGA0_nav ...
                  		cis_nav i0_nav crc_nav omega_nav OMEGA_DOT_nav idot_nav CodeL2_nav GPSwknum_nav ...
                  		L2Pdataflag_nav SVaccuracy_nav SVhealth_nav tgd_nav IODC_nav ...
      						txtimeofmessage_nav fitinterval_nav toc_nav_week];  %note added toc_nav_week 6/9/05 TB.
                        
                            
                            
        %true x , y , z position
   
   TrueX(i) = a_nav(4);
   TrueY(i) = a_nav(5);
   TrueZ(i) = a_nav(6);
                                 
   TruePosVec(i,1:3) = [TrueX(i) TrueY(i) TrueZ(i)];                          
      
      
    
    i = i + 1; 
     
      
end
   
end %if readheaderflag ==1 
   
   fclose(filename_nav);






























