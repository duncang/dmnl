%This is for seeing the effect of different geometries
%24 November 2008



warning off; 



%USE THIS ONE FOR PLOTTING MONTECARLO RESULTS

%For one type of ramp fault, 2.5 m/s , it calculates with and without MMF. 

feature accel on


%need this for HPC


%matlabfiles\functions

%add paths
% DataPath = 'matlabfiles/data/ThesisFlight/';
% DataPath = 'matlabfiles/functions/';
% DataPath = 'matlabfiles/scripts/';

%cd('/home/n2523710/matlabfiles/')
%cd('L:home/n2523710/matlabfiles/data/ThesisFlight')
%cd('/home/n2523710/matlabfiles/functions/')


%path(path,'/home/n2523710/matlabfiles/functions/');


% 
% 
% path(path,'/home/n2523710/matlabfiles/');
% path(path,'/home/n2523710/matlabfiles/scripts/');
% path(path,'/home/n2523710/matlabfiles/scripts/data/');


%function which calls the main function with different statistics 

tstart = cputime;

%start
%number of montecarlo runs

Nmc = 6;   %83 runs is the maximum without getting an out of memory error


%change INS statistics
%use seed (rand(n))

%change GPS statistics

%change ADM statistics
MonteCarlo = 1;

%end


%do geometry simulation
montecarlo_geometry = 1;



runload = 1;





%now preallocate these with zeros


%amount to allocate



%these need to be the same as in the other file ie julier 3.m
startepoch = 3;
endepoch = 137;




 timeallocate = endepoch;
%save relevant parameters

HPL_INS_NFmc  = zeros(Nmc,timeallocate);
HPL_INS_Cmc = zeros(Nmc,timeallocate);
VPL_INS_NFmc = zeros(Nmc,timeallocate);
VPL_INS_Cmc = zeros(Nmc,timeallocate);



HPL_MODEL_NFmc = zeros(Nmc,timeallocate);
HPL_MODEL_Cmc = zeros(Nmc,timeallocate);
VPL_MODEL_NFmc = zeros(Nmc,timeallocate);
VPL_MODEL_Cmc = zeros(Nmc,timeallocate);



HPL_GPS_mc = zeros(Nmc,timeallocate);
VPL_GPS_mc = zeros(Nmc,timeallocate);




%save HPL and VPL H0

HPL_H0_INS_NFmc = zeros(Nmc,timeallocate);
HPL_H0_INS_Cmc = zeros(Nmc,timeallocate);
VPL_H0_INS_NFmc = zeros(Nmc,timeallocate);
VPL_H0_INS_Cmc = zeros(Nmc,timeallocate);


HPL_H0_MODEL_NFmc = zeros(Nmc,timeallocate);
HPL_H0_MODEL_Cmc = zeros(Nmc,timeallocate);
VPL_H0_MODEL_NFmc = zeros(Nmc,timeallocate);
VPL_H0_MODEL_Cmc = zeros(Nmc,timeallocate);


HPL_GPS_H0mc = zeros(Nmc,timeallocate);
VPL_GPS_H0mc = zeros(Nmc,timeallocate);


%save test statistics and thresholds




aTDINS_Hmc = zeros(Nmc,timeallocate);
aTDINS_Vmc = zeros(Nmc,timeallocate);





lambda_ss_outINS_H_NFmc = zeros(Nmc, timeallocate,32);
lambda_ss_outINS_Hmc   = zeros(Nmc, timeallocate,32);


lambda_ss_outINS_V_NFmc  = zeros(Nmc, timeallocate,32);
lambda_ss_outINS_Vmc  = zeros(Nmc, timeallocate,32);


lambda_ss_outMODEL_H_NFmc  = zeros(Nmc, timeallocate,32);
lambda_ss_outMODEL_Hmc  = zeros(Nmc, timeallocate,32);

lambda_ss_outMODEL_V_NFmc = zeros(Nmc, timeallocate,32);
lambda_ss_outMODEL_Vmc = zeros(Nmc, timeallocate,32);

lambda_ss_outGPS_H_NFmc  = zeros(Nmc, timeallocate,32);
lambda_ss_outGPS_Hmc  = zeros(Nmc, timeallocate,32);

lambda_ss_outGPS_V_NFmc  = zeros(Nmc, timeallocate,32);
lambda_ss_outGPS_Vmc  = zeros(Nmc, timeallocate,32);

                       
             



HPE_INS_NFmc = zeros(Nmc,timeallocate);
HPE_INS_Cmc = zeros(Nmc,timeallocate);
VPE_INS_NFmc = zeros(Nmc,timeallocate); 
VPE_INS_Cmc = zeros(Nmc,timeallocate);



HPE_MODEL_NFmc = zeros(Nmc,timeallocate);
HPE_MODEL_Cmc = zeros(Nmc,timeallocate);
VPE_MODEL_NFmc = zeros(Nmc,timeallocate);
VPE_MODEL_Cmc = zeros(Nmc,timeallocate);






 HPE_GPS_mc = zeros(Nmc,timeallocate);
 VPE_GPS_mc = zeros(Nmc,timeallocate);




%for INS
 sigma_vCmc_INS = zeros(Nmc,timeallocate);
        sigma_v_maxCmc_INS  = zeros(Nmc,timeallocate);
        dmajorCmc_INS = zeros(Nmc,timeallocate);
        dmajor_maxCmc_INS = zeros(Nmc,timeallocate);
    

      %for ADM
        sigma_v_Cmc_M = zeros(Nmc,timeallocate);
        sigma_v_max_Cmc_M  = zeros(Nmc,timeallocate);    
        dmajor_Cmc_M = zeros(Nmc,timeallocate);
        dmajor_max_Cmc_M = zeros(Nmc,timeallocate);



%for INS unfused
 sigma_vNFmc_INS = zeros(Nmc,timeallocate);
        sigma_v_maxNFmc_INS  = zeros(Nmc,timeallocate);
        dmajorNFmc_INS = zeros(Nmc,timeallocate);
        dmajor_maxNFmc_INS = zeros(Nmc,timeallocate);
    

        
 %for ADM unfused     
        sigma_v_NFmc_M=  zeros(Nmc,timeallocate);
        sigma_v_max_NFmc_M = zeros(Nmc,timeallocate);  
        dmajor_NFmc_M = zeros(Nmc,timeallocate);
        dmajor_max_NFmc_M = zeros(Nmc,timeallocate);
        
        





x_state_total_INSBIAS_NFmc =  zeros(Nmc,6, timeallocate);
x_state_total_INSBIAS_Cmc =  zeros(Nmc,6, timeallocate);


x_state_total_MODELBIAS_NFmc =  zeros(Nmc,6, timeallocate);
x_state_total_MODELBIAS_Cmc =  zeros(Nmc,6, timeallocate);







%positions
N_PosErrorsINS_NFmc = zeros(Nmc,timeallocate);
E_PosErrorsINS_NFmc= zeros(Nmc,timeallocate);
D_PosErrorsINS_NFmc = zeros(Nmc,timeallocate);

N_PosErrorsMODEL_NFmc = zeros(Nmc,timeallocate);
E_PosErrorsMODEL_NFmc = zeros(Nmc,timeallocate);
D_PosErrorsMODEL_NFmc = zeros(Nmc,timeallocate);

N_PosErrorsINS_Cmc = zeros(Nmc,timeallocate);
E_PosErrorsINS_Cmc = zeros(Nmc,timeallocate);
D_PosErrorsINS_Cmc = zeros(Nmc,timeallocate);

N_PosErrorsMODEL_Cmc = zeros(Nmc,timeallocate);
E_PosErrorsMODEL_Cmc = zeros(Nmc,timeallocate);
D_PosErrorsMODEL_Cmc = zeros(Nmc,timeallocate);


%velocities

N_VelErrorsINS_NFmc = zeros(Nmc,timeallocate);
E_VelErrorsINS_NFmc = zeros(Nmc,timeallocate);
D_VelErrorsINS_NFmc = zeros(Nmc,timeallocate);


N_VelErrorsMODEL_NFmc = zeros(Nmc,timeallocate);
E_VelErrorsMODEL_NFmc = zeros(Nmc,timeallocate);
D_VelErrorsMODEL_NFmc = zeros(Nmc,timeallocate);

N_VelErrorsINS_Cmc = zeros(Nmc,timeallocate);
E_VelErrorsINS_Cmc = zeros(Nmc,timeallocate);
D_VelErrorsINS_Cmc = zeros(Nmc,timeallocate);

N_VelErrorsMODEL_Cmc = zeros(Nmc,timeallocate);
E_VelErrorsMODEL_Cmc = zeros(Nmc,timeallocate);
D_VelErrorsMODEL_Cmc = zeros(Nmc,timeallocate);

%commented out the correlation stuff, not necessary for thesis results 
%
% corrcoeffAtt_NFmc = zeros(Nmc,3,timeallocate);
% corrcoeffVel_NFmc = zeros(Nmc,3,timeallocate);
% corrcoeffPos_NFmc = zeros(Nmc,3,timeallocate);
% 
% corrcoeffCLKbias_NFmc = zeros(Nmc,1,timeallocate);
% corrcoeffCLKdrift_NFmc = zeros(Nmc,1,timeallocate);
% 
% corrcoeffAccbias_NFmc = zeros(Nmc,3,timeallocate);
% corrcoeffGyrobias_NFmc = zeros(Nmc,3,timeallocate);
% 
% 
% corrcoeffAtt_Cmc = zeros(Nmc,3,timeallocate);
% corrcoeffVel_Cmc = zeros(Nmc,3,timeallocate);
% corrcoeffPos_Cmc = zeros(Nmc,3,timeallocate);
% 
% corrcoeffCLKbias_Cmc = zeros(Nmc,1,timeallocate);
% corrcoeffCLKdrift_Cmc = zeros(Nmc,1,timeallocate);
% 
% corrcoeffAccbias_Cmc = zeros(Nmc,3,timeallocate);
% corrcoeffGyrobias_Cmc = zeros(Nmc,3,timeallocate);




  RollErrorINS_NFmc = zeros(Nmc,timeallocate);
   PitchErrorINS_NFmc =  zeros(Nmc,timeallocate); 
   YawErrorINS_NFmc = zeros(Nmc,timeallocate);
   
   RollErrorMODEL_NFmc = zeros(Nmc,timeallocate);
   PitchErrorMODEL_NFmc = zeros(Nmc,timeallocate);
   YawErrorMODEL_NFmc = zeros(Nmc,timeallocate);


  RollErrorINS_Cmc = zeros(Nmc,timeallocate);
   PitchErrorINS_Cmc= zeros(Nmc,timeallocate);
   YawErrorINS_Cmc = zeros(Nmc,timeallocate);
   
   RollErrorMODEL_Cmc = zeros(Nmc,timeallocate);
   PitchErrorMODEL_Cmc = zeros(Nmc,timeallocate);
   YawErrorMODEL_Cmc = zeros(Nmc,timeallocate);



 Pminus_save_NFmc =  zeros(Nmc,34,34,timeallocate);

 P_save_NFmc  =  zeros(Nmc,34,34,timeallocate);


x_hat_save_NFmc =  zeros(Nmc,34,timeallocate);
V_C_save_NFmc =  zeros(Nmc,4*12+3,4*12+3,timeallocate);  %assume max number of sats is 12
z_save_NFmc = zeros(Nmc,timeallocate,4*12+3);


 
 
%with fusion 
 Pminus_save_Cmc =  zeros(Nmc,34,34,timeallocate);

 P_save_Cmc  =  zeros(Nmc,34,34,timeallocate);
 

x_hat_save_Cmc = zeros(Nmc,34,timeallocate);
V_C_save_Cmc =  zeros(Nmc,4*12+3,4*12+3,timeallocate); 
z_save_Cmc = zeros(Nmc,timeallocate,4*12+3);
N_save_Cmc  = zeros(Nmc,timeallocate);



%I use these primarily for estimating Qk1k2, also for verifying P minus
x_state_total_minusINSsaveCmc =  zeros(Nmc,10,timeallocate);
x_state_total_minusMODELsaveCmc =  zeros(Nmc,10,timeallocate);


x_state_total_minusINSsaveNFmc =  zeros(Nmc,10,timeallocate);
x_state_total_minusMODELsaveNFmc =  zeros(Nmc,10,timeallocate);










for iM = 1:Nmc
    
    
       
    
    iM
    
    %this is to make code only run once in the newjulier3 code
%     if iM == 1
%         
%         runOnce = 0;
%         
%         
%     else
%         runOnce = 1;
%         
%     end
%     
    
  
    
    
    
    
%clear variables need to do this to free up memory
 close all %close any open matlab figures    
    
%do this to try and free up memory
    
clear HPL_INS_NF HPL_INS_C VPL_INS_NF clearVPL_INS_C HPL_MODEL_NF HPL_MODEL_C VPL_MODEL_NF VPL_MODEL_C HPE_INS_NF HPE_INS_C VPE_INS_NF VPE_INS_C HPE_MODEL_NF HPE_MODEL_C VPE_MODEL_NF
clear  VPE_MODEL_C N_PosErrorsINS_NF E_PosErrorsINS_NF D_PosErrorsINS_NF N_PosErrorsMODEL_NF E_PosErrorsMODEL_NF D_PosErrorsMODEL_NF N_PosErrorsINS_C E_PosErrorsINS_C D_PosErrorsINS_C
clear  N_PosErrorsMODEL_C E_PosErrorsMODEL_C D_PosErrorsMODEL_C N_VelErrorsINS_NF E_VelErrorsINS_NF D_VelErrorsINS_NF N_VelErrorsMODEL_NF E_VelErrorsMODEL_NF D_VelErrorsMODEL_NF
clear  N_VelErrorsINS_C E_VelErrorsINS_C D_VelErrorsINS_C N_VelErrorsMODEL_C E_VelErrorsMODEL_C D_VelErrorsMODEL_C corrcoeffAtt_NF corrcoeffVel_NF corrcoeffPos_NF
clear corrcoeffCLKbias_NF corrcoeffCLKdrift_NF corrcoeffAccbias_NF corrcoeffGyrobias_NF corrcoeffAtt corrcoeffVel corrcoeffPos corrcoeffCLKbias corrcoeffCLKdrift corrcoeffAccbias
clear  corrcoeffGyrobias RollErrorINS_NF PitchErrorINS_NF YawErrorINS_NF RollErrorMODEL_NF PitchErrorMODEL_NF YawErrorMODEL_NF RollErrorINS_C PitchErrorINS_C YawErrorINS_C
clear RollErrorMODEL_C RollErrorMODEL_C YawErrorMODEL_C Pminus_save_NF P_minus_C100Hz_save_NF P_save_NF x_hat_save_NF V_C_save_NF z_save_NF N_save_NF
clear  P_minus_save P_minus_C100Hz_save P_save x_hat_save V_C_save z_save N_save    
    
    %not enough physical memory RAM
    
    
%call the script

%for yy = 1:2  %with this it was running it 4 times instead of twice...
    
    
Fused = 0;  

runSensors = 1;  %use a new random number generator seed in airdata sensors



%GPSINSMODEL_ONLY_EKF100HzNewJulier3;
%GPSINSMODEL_ONLY_EKF100HzNewJulier2withAdaptiveFilter;
%GPSINSMODEL_ONLY_EKF100HzNewJulier2;

GPSINSMODEL_ONLY_EKF100HzNewJulier3NoBias;

%GPSTWOINS_EKF100HzNewJulier;
runload = 0;  

Fused = 1;
runSensors = 0;



%GPSINSMODEL_ONLY_EKF100HzNewJulier3;

GPSINSMODEL_ONLY_EKF100HzNewJulier3NoBias;
%GPSINSMODEL_ONLY_EKF100HzNewJulier2withAdaptiveFilter;
%GPSINSMODEL_ONLY_EKF100HzNewJulier2;
%GPSTWOINS_EKF100HzNewJulier;
%end

%save relevant parameters

HPL_INS_NFmc(iM,:) = HPL_INS_NF;
HPL_INS_Cmc(iM,:) = HPL_INS_C;
VPL_INS_NFmc(iM,:) =  VPL_INS_NF;
VPL_INS_Cmc(iM,:) = VPL_INS_C;



HPL_MODEL_NFmc(iM,:) = HPL_MODEL_NF;
HPL_MODEL_Cmc(iM,:) = HPL_MODEL_C;
VPL_MODEL_NFmc(iM,:) = VPL_MODEL_NF;
VPL_MODEL_Cmc(iM,:) = VPL_MODEL_C;



HPL_GPS_mc(iM,:) = HPL_GPS;
VPL_GPS_mc(iM,:) = VPL_GPS;




%save HPL and VPL H0

HPL_H0_INS_NFmc(iM,:) = HPL_INS_H0_NF;
HPL_H0_INS_Cmc(iM,:) = HPL_H0;
VPL_H0_INS_NFmc(iM,:) = VPL_INS_H0_NF;
VPL_H0_INS_Cmc(iM,:) = VPL_H0;


HPL_H0_MODEL_NFmc(iM,:) = HPL_MODEL_H0_NF;
HPL_H0_MODEL_Cmc(iM,:) = HPL_H0_M;
VPL_H0_MODEL_NFmc(iM,:) = VPL_MODEL_H0_NF;
VPL_H0_MODEL_Cmc(iM,:) = VPL_H0_M;


HPL_GPS_H0mc(iM,:) = HPL_H0_GPS;
VPL_GPS_H0mc(iM,:) = VPL_H0_GPS;


%save test statistics and thresholds



aTDINS_Hmc(iM,:) = aTDINS_H;
aTDINS_Vmc(iM,:) = aTDINS_V;



lambda_ss_outINS_H_NFmc(iM,:,:) = lambda_ss_outINS_H_NF;
lambda_ss_outINS_Hmc(iM,:,:) = lambda_ss_outINS_H;


lambda_ss_outINS_V_NFmc(iM,:,:) = lambda_ss_outINS_V_NF;
lambda_ss_outINS_Vmc(iM,:,:) = lambda_ss_outINS_V;


lambda_ss_outMODEL_H_NFmc(iM,:,:) = lambda_ss_outMODEL_H_NF;
lambda_ss_outMODEL_Hmc(iM,:,:) = lambda_ss_outMODEL_H;

lambda_ss_outMODEL_V_NFmc(iM,:,:) = lambda_ss_outMODEL_V_NF;
lambda_ss_outMODEL_Vmc(iM,:,:) = lambda_ss_outMODEL_V;

lambda_ss_outGPS_H_NFmc(iM,:,:) = lambda_ss_outGPS_H_NF;
lambda_ss_outGPS_Hmc(iM,:,:) = lambda_ss_outGPS_H;

lambda_ss_outGPS_V_NFmc(iM,:,:) = lambda_ss_outGPS_V_NF;
lambda_ss_outGPS_Vmc(iM,:,:) = lambda_ss_outGPS_V;




MaxSlopeSatelliteToDo_mc(iM,1) = MaxSlopeSatelliteToDo;



% 
% 
%         lambda_ss_outINS_H_NF = lambda_ss_outINS_H;
%          lambda_ss_outINS_V_NF = lambda_ss_outINS_V;
%         
%         
%              lambda_ss_outMODEL_H_NF = lambda_ss_outMODEL_H;
%              lambda_ss_outMODEL_V_NF = lambda_ss_outMODEL_V;
%              
%              
%         
%              lambda_ss_outGPS_H_NF = lambda_ss_outGPS_H;
%              lambda_ss_outGPS_V_NF = lambda_ss_outGPS_V;
%         

             
             
             
             
             



HPE_INS_NFmc(iM,:) = HPE_INS_NF;
HPE_INS_Cmc(iM,:) = HPE_INS_C;
VPE_INS_NFmc(iM,:) = VPE_INS_NF;
VPE_INS_Cmc(iM,:) = VPE_INS_C;



HPE_MODEL_NFmc(iM,:) = HPE_MODEL_NF;
HPE_MODEL_Cmc(iM,:) = HPE_MODEL_C;
VPE_MODEL_NFmc(iM,:) = VPE_MODEL_NF;
VPE_MODEL_Cmc(iM,:) = VPE_MODEL_C;








 HPE_GPS_mc(iM,:) = HPE_GPS_C;    %there's no difference between fused and not fused for the GPS
 VPE_GPS_mc(iM,:) = VPE_GPS_C;




%for INS
 sigma_vCmc_INS(iM,:) = sigma_v_INS;
        sigma_v_maxCmc_INS(iM,:) = sigma_v_max_INS;        
        dmajorCmc_INS(iM,:) = dmajor_INS;
        dmajor_maxCmc_INS(iM,:) = dmajor_max_INS;
    

      %for ADM
        sigma_v_Cmc_M(iM,:)  = sigma_v_M;
        sigma_v_max_Cmc_M(iM,:) = sigma_v_max_M;        
        dmajor_Cmc_M(iM,:) = dmajor_M;
        dmajor_max_Cmc_M(iM,:) = dmajor_max_M;



%for INS unfused
 sigma_vNFmc_INS(iM,:) = sigma_v_INSNF;
        sigma_v_maxNFmc_INS(iM,:) = sigma_v_max_INSNF;        
        dmajorNFmc_INS(iM,:) = dmajor_INSNF;
        dmajor_maxNFmc_INS(iM,:) = dmajor_max_INSNF;
    

        
 %for ADM unfused     
        sigma_v_NFmc_M(iM,:)  = sigma_v_MNF;
        sigma_v_max_NFmc_M(iM,:) = sigma_v_max_MNF;        
        dmajor_NFmc_M(iM,:) = dmajor_MNF;
        dmajor_max_NFmc_M(iM,:) = dmajor_max_MNF;
        
        





x_state_total_INSBIAS_NFmc(iM,:,:) =  x_state_total_INSBIAS_NF(:,1:endepoch);
x_state_total_INSBIAS_Cmc(iM,:,:) =  x_state_total_INSBIAS(:,1:endepoch);


x_state_total_MODELBIAS_NFmc(iM,:,:) =  x_state_total_MODELBIAS_NF(:,1:endepoch);
x_state_total_MODELBIAS_Cmc(iM,:,:) =  x_state_total_MODELBIAS(:,1:endepoch);



%need these to calculate ensemble variance

%positions
N_PosErrorsINS_NFmc(iM,:) = N_PosErrorsINS_NF;
E_PosErrorsINS_NFmc(iM,:) = E_PosErrorsINS_NF;
D_PosErrorsINS_NFmc(iM,:) = D_PosErrorsINS_NF;

N_PosErrorsMODEL_NFmc(iM,:) = N_PosErrorsMODEL_NF;
E_PosErrorsMODEL_NFmc(iM,:) = E_PosErrorsMODEL_NF;
D_PosErrorsMODEL_NFmc(iM,:) = D_PosErrorsMODEL_NF;


N_PosErrorsINS_Cmc(iM,:) = N_PosErrorsINS_C;
E_PosErrorsINS_Cmc(iM,:) = E_PosErrorsINS_C;
D_PosErrorsINS_Cmc(iM,:) = D_PosErrorsINS_C;

N_PosErrorsMODEL_Cmc(iM,:) = N_PosErrorsMODEL_C;
E_PosErrorsMODEL_Cmc(iM,:) = E_PosErrorsMODEL_C;
D_PosErrorsMODEL_Cmc(iM,:) = D_PosErrorsMODEL_C;


%velocities

N_VelErrorsINS_NFmc(iM,:) = N_VelErrorsINS_NF;
E_VelErrorsINS_NFmc(iM,:) = E_VelErrorsINS_NF;
D_VelErrorsINS_NFmc(iM,:) = D_VelErrorsINS_NF;


N_VelErrorsMODEL_NFmc(iM,:) = N_VelErrorsMODEL_NF;
E_VelErrorsMODEL_NFmc(iM,:) = E_VelErrorsMODEL_NF;
D_VelErrorsMODEL_NFmc(iM,:) = D_VelErrorsMODEL_NF;

N_VelErrorsINS_Cmc(iM,:) = N_VelErrorsINS_C;
E_VelErrorsINS_Cmc(iM,:) = E_VelErrorsINS_C;
D_VelErrorsINS_Cmc(iM,:) = D_VelErrorsINS_C;

N_VelErrorsMODEL_Cmc(iM,:) = N_VelErrorsMODEL_C;
E_VelErrorsMODEL_Cmc(iM,:)= E_VelErrorsMODEL_C;
D_VelErrorsMODEL_Cmc(iM,:) = D_VelErrorsMODEL_C;


% 
% corrcoeffAtt_NFmc(iM,:,:) = corrcoeffAtt_NF;
% corrcoeffVel_NFmc(iM,:,:) = corrcoeffVel_NF;
% corrcoeffPos_NFmc(iM,:,:) = corrcoeffPos_NF;
% 
% corrcoeffCLKbias_NFmc(iM,:,:) = corrcoeffCLKbias_NF;
% corrcoeffCLKdrift_NFmc(iM,:,:) = corrcoeffCLKdrift_NF;
% 
% corrcoeffAccbias_NFmc(iM,:,:) = corrcoeffAccbias_NF;
% corrcoeffGyrobias_NFmc(iM,:,:) = corrcoeffGyrobias_NF;
% 
% 
% corrcoeffAtt_Cmc(iM,:,:) = corrcoeffAtt;
% corrcoeffVel_Cmc(iM,:,:) = corrcoeffVel;
% corrcoeffPos_Cmc(iM,:,:) = corrcoeffPos;
% 
% corrcoeffCLKbias_Cmc(iM,:,:) = corrcoeffCLKbias;
% corrcoeffCLKdrift_Cmc(iM,:,:) = corrcoeffCLKdrift;
% 
% corrcoeffAccbias_Cmc(iM,:,:) = corrcoeffAccbias;
% corrcoeffGyrobias_Cmc(iM,:,:) = corrcoeffGyrobias;




  RollErrorINS_NFmc(iM,:) = RollErrorINS_NF;
   PitchErrorINS_NFmc(iM,:) = PitchErrorINS_NF;
   YawErrorINS_NFmc(iM,:) = YawErrorINS_NF;  
   
   RollErrorMODEL_NFmc(iM,:) = RollErrorMODEL_NF;
   PitchErrorMODEL_NFmc(iM,:) = PitchErrorMODEL_NF;
   YawErrorMODEL_NFmc(iM,:) = YawErrorMODEL_NF; 


  RollErrorINS_Cmc(iM,:) = RollErrorINS_C;
   PitchErrorINS_Cmc(iM,:) = PitchErrorINS_C;
   YawErrorINS_Cmc(iM,:) = YawErrorINS_C;  
   
   RollErrorMODEL_Cmc(iM,:) = RollErrorMODEL_C;
   PitchErrorMODEL_Cmc(iM,:) = PitchErrorMODEL_C;
   YawErrorMODEL_Cmc(iM,:) = YawErrorMODEL_C; 



 Pminus_save_NFmc(iM,:,:,:) = Pminus_save_NF;
 %P_minus_C100Hz_save_NFmc(iM,:,:,:) = P_minus_C100Hz_save_NF;
 P_save_NFmc(iM,:,:,:) = P_save_NF;
%  Q_save100Hz_NF = Q_save100Hz; 
%  W_INS_Save_NF = W_INS_Save; 
%  W_MODEL_Save_NF = W_MODEL_Save;  

x_hat_save_NFmc(iM,:,:,:) = x_hat_save_NF;
V_C_save_NFmc(iM,:,:,:) = V_C_save_NF;
z_save_NFmc(iM,:,:,:) = z_save_NF;


 
 
%with fusion 
 Pminus_save_Cmc(iM,:,:,:) = P_minus_save;
 %P_minus_C100Hz_save_Cmc(iM,:,:,:)  = P_minus_C100Hz_save;
 P_save_Cmc(iM,:,:,:)  = P_save; 
 

x_hat_save_Cmc(iM,:,:,:)  = x_hat_save;
V_C_save_Cmc(iM,:,:,:)  = V_C_save;
z_save_Cmc(iM,:,:,:)  = z_save;
N_save_Cmc(iM,:,:,:)  = N_save;



%I use these primarily for estimating Qk1k2, also for verifying P minus
x_state_total_minusINSsaveCmc(iM,:,:) = x_state_total_minusINSsave;
x_state_total_minusMODELsaveCmc(iM,:,:) = x_state_total_minusMODELsave;


x_state_total_minusINSsaveNFmc(iM,:,:) = x_state_total_minusINSsaveNF;
x_state_total_minusMODELsaveNFmc(iM,:,:) = x_state_total_minusMODELsaveNF;









%save('MonteCarloData.mat');   %save all variables in workspace in current directory

end


save('MonteCarloDataTest.mat');   %at end, save all variables in workspace in current directory


%quit


%i just commented it all so matlab doesnt have to run


% 
% processdata = 0; %don't process the data... (if running on HPC etc)
% 
% if processdata == 1
% 
% 
% booboo = aart;  %just to stop code from going on by making it give an error
% 
% %===================
% %PROCESS DATA 
% %===================
% 
% startprocessepoch = startepoch + 20;  %start it after the filter has reached steady state
% endprocessepoch = endepoch-3;
% 
% %Calculate mean values  state errors
% 
% 
% %Position
% 
% 
% for ko = startprocessepoch:endprocessepoch
% 
%     
%     
%    %I think this is right, should result in an ix6 vector
% % meanx_state_total_INSBIAS_NFmc(ko,1)   = mean(x_state_total_INSBIAS_NFmc(:,1,ko));
% % meanx_state_total_INSBIAS_NFmc(ko,2)   = mean(x_state_total_INSBIAS_NFmc(:,2,ko));
% % meanx_state_total_INSBIAS_NFmc(ko,3)   = mean(x_state_total_INSBIAS_NFmc(:,3,ko));
% % meanx_state_total_INSBIAS_NFmc(ko,4)   = mean(x_state_total_INSBIAS_NFmc(:,4,ko));
% % meanx_state_total_INSBIAS_NFmc(ko,5)   = mean(x_state_total_INSBIAS_NFmc(:,5,ko));
% % meanx_state_total_INSBIAS_NFmc(ko,6)   = mean(x_state_total_INSBIAS_NFmc(:,6,ko));
% % 
% % 
% % 
% % meanx_state_total_INSBIAS_Cmc(ko,1) = mean(x_state_total_INSBIAS_Cmc(:,1,ko));
% % meanx_state_total_INSBIAS_Cmc(ko,2) = mean(x_state_total_INSBIAS_Cmc(:,2,ko));
% % meanx_state_total_INSBIAS_Cmc(ko,3) = mean(x_state_total_INSBIAS_Cmc(:,3,ko));
% % meanx_state_total_INSBIAS_Cmc(ko,4) = mean(x_state_total_INSBIAS_Cmc(:,4,ko));
% % meanx_state_total_INSBIAS_Cmc(ko,5) = mean(x_state_total_INSBIAS_Cmc(:,5,ko));
% % meanx_state_total_INSBIAS_Cmc(ko,6) = mean(x_state_total_INSBIAS_Cmc(:,6,ko));
% % 
% % 
% % 
% % meanx_state_total_MODELBIAS_NFmc(ko,1)  = mean(x_state_total_MODELBIAS_NFmc(:,1,ko)); 
% % meanx_state_total_MODELBIAS_NFmc(ko,2)  = mean(x_state_total_MODELBIAS_NFmc(:,2,ko)); 
% % meanx_state_total_MODELBIAS_NFmc(ko,3)  = mean(x_state_total_MODELBIAS_NFmc(:,3,ko)); 
% % meanx_state_total_MODELBIAS_NFmc(ko,4)  = mean(x_state_total_MODELBIAS_NFmc(:,4,ko)); 
% % meanx_state_total_MODELBIAS_NFmc(ko,5)  = mean(x_state_total_MODELBIAS_NFmc(:,5,ko)); 
% % meanx_state_total_MODELBIAS_NFmc(ko,6)  = mean(x_state_total_MODELBIAS_NFmc(:,6,ko)); 
% % 
% % 
% % 
% % meanx_state_total_MODELBIAS_Cmc(ko,1)  = mean(x_state_total_MODELBIAS_Cmc(:,1,ko)) ;
% % meanx_state_total_MODELBIAS_Cmc(ko,2)  = mean(x_state_total_MODELBIAS_Cmc(:,2,ko)) ;
% % meanx_state_total_MODELBIAS_Cmc(ko,3)  = mean(x_state_total_MODELBIAS_Cmc(:,3,ko)) ;
% % meanx_state_total_MODELBIAS_Cmc(ko,4)  = mean(x_state_total_MODELBIAS_Cmc(:,4,ko)) ;
% % meanx_state_total_MODELBIAS_Cmc(ko,5)  = mean(x_state_total_MODELBIAS_Cmc(:,5,ko)) ;
% % meanx_state_total_MODELBIAS_Cmc(ko,6)  = mean(x_state_total_MODELBIAS_Cmc(:,6,ko)) ;
% % 
% 
% 
%     
%     
% 
% meanHPL_INS_NFmc(ko) = mean(HPL_INS_NFmc(:,ko));
% meanHPL_INS_Cmc(ko) =  mean(HPL_INS_Cmc(:,ko));
% meanVPL_INS_NFmc(ko) =  mean(VPL_INS_NFmc(:,ko));
% meanVPL_INS_Cmc(ko) =  mean(VPL_INS_Cmc(:,ko));
% 
% meanHPL_MODEL_NFmc(ko) =  mean(HPL_MODEL_NFmc(:,ko));
% meanHPL_MODEL_Cmc(ko) =  mean(HPL_MODEL_Cmc(:,ko));
% meanVPL_MODEL_NFmc(ko) =  mean(VPL_MODEL_NFmc(:,ko));
% meanVPL_MODEL_Cmc(ko) = mean( VPL_MODEL_Cmc(:,ko));
% 
% 
% 
% meanHPE_INS_NFmc(ko) = mean(HPE_INS_NFmc(:,ko));
% meanHPE_INS_Cmc(ko) = mean(HPE_INS_Cmc(:,ko));
% meanVPE_INS_NFmc(ko) = mean(VPE_INS_NFmc(:,ko));
% meanVPE_INS_Cmc(ko) = mean(VPE_INS_Cmc(:,ko));
% 
% meanHPE_MODEL_NFmc(ko) = mean(HPE_MODEL_NFmc(:,ko));
% meanHPE_MODEL_Cmc(ko) = mean(HPE_MODEL_Cmc(:,ko));
% meanVPE_MODEL_NFmc(ko) = mean(VPE_MODEL_NFmc(:,ko));
% meanVPE_MODEL_Cmc(ko) = mean(VPE_MODEL_Cmc(:,ko));
% 
% % 
% % meanHPE_GPS_mc(ko) = mean(HPE_GPS_mc(:,ko));
% % meanVPE_GPS_mc(ko) = mean(VPE_GPS_mc(:,ko));
% 
% 
% 
% meanN_PosErrorsINS_NFmc(ko) = mean(N_PosErrorsINS_NFmc(:,ko));
% meanE_PosErrorsINS_NFmc(ko) = mean(E_PosErrorsINS_NFmc(:,ko));
% meanD_PosErrorsINS_NFmc(ko) = mean(D_PosErrorsINS_NFmc(:,ko));
% 
% 
% meanN_PosErrorsMODEL_NFmc(ko) = mean(N_PosErrorsMODEL_NFmc(:,ko));
% meanE_PosErrorsMODEL_NFmc(ko) = mean(E_PosErrorsMODEL_NFmc(:,ko));
% meanD_PosErrorsMODEL_NFmc(ko) = mean(D_PosErrorsMODEL_NFmc(:,ko));
% 
% 
% meanN_PosErrorsINS_Cmc(ko) = mean(N_PosErrorsINS_Cmc(:,ko));
% meanE_PosErrorsINS_Cmc(ko) = mean(E_PosErrorsINS_Cmc(:,ko));
% meanD_PosErrorsINS_Cmc(ko) = mean(D_PosErrorsINS_Cmc(:,ko));
% 
% 
% meanN_PosErrorsMODEL_Cmc(ko) = mean(N_PosErrorsMODEL_Cmc(:,ko));
% meanE_PosErrorsMODEL_Cmc(ko) = mean(E_PosErrorsMODEL_Cmc(:,ko));
% meanD_PosErrorsMODEL_Cmc(ko) = mean(D_PosErrorsMODEL_Cmc(:,ko));
% 
% meanN_VelErrorsINS_NFmc(ko) = mean(N_VelErrorsINS_NFmc(:,ko));
% meanE_VelErrorsINS_NFmc(ko) = mean(E_VelErrorsINS_NFmc(:,ko));
% meanD_VelErrorsINS_NFmc(ko) = mean(D_VelErrorsINS_NFmc(:,ko));
% 
% meanN_VelErrorsMODEL_NFmc(ko) = mean(N_VelErrorsMODEL_NFmc(:,ko));
% meanE_VelErrorsMODEL_NFmc(ko) = mean(E_VelErrorsMODEL_NFmc(:,ko));
% meanD_VelErrorsMODEL_NFmc(ko) = mean(D_VelErrorsMODEL_NFmc(:,ko));
% 
% meanN_VelErrorsINS_Cmc(ko) = mean(N_VelErrorsINS_Cmc(:,ko));
% meanE_VelErrorsINS_Cmc(ko) = mean(E_VelErrorsINS_Cmc(:,ko));
% meanD_VelErrorsINS_Cmc(ko) = mean(D_VelErrorsINS_Cmc(:,ko));
% 
% 
% meanN_VelErrorsMODEL_Cmc(ko) = mean(N_VelErrorsMODEL_Cmc(:,ko));
% meanE_VelErrorsMODEL_Cmc(ko) = mean(E_VelErrorsMODEL_Cmc(:,ko));
% meanD_VelErrorsMODEL_Cmc(ko) = mean(D_VelErrorsMODEL_Cmc(:,ko));
% 
% 
% meanRollErrorINS_NFmc(ko) = mean(RollErrorINS_NFmc(:,ko));
% meanPitchErrorINS_NFmc(ko) = mean(PitchErrorINS_NFmc(:,ko));
% meanYawErrorINS_NFmc(ko) = mean(YawErrorINS_NFmc(:,ko));
% 
% meanRollErrorMODEL_NFmc(ko) = mean(RollErrorMODEL_NFmc(:,ko));
% meanPitchErrorMODEL_NFmc(ko) = mean(PitchErrorMODEL_NFmc(:,ko));
% meanYawErrorMODEL_NFmc(ko) = mean(YawErrorMODEL_NFmc(:,ko));
% 
% 
% meanRollErrorINS_Cmc(ko) = mean(RollErrorINS_Cmc(:,ko));
% meanPitchErrorINS_Cmc(ko) = mean(PitchErrorINS_Cmc(:,ko));
% meanYawErrorINS_Cmc(ko) = mean(YawErrorINS_Cmc(:,ko));
% 
% 
% meanRollErrorMODEL_Cmc(ko) = mean(RollErrorMODEL_Cmc(:,ko));
% meanPitchErrorMODEL_Cmc(ko) = mean(PitchErrorMODEL_Cmc(:,ko));
% meanYawErrorMODEL_Cmc(ko) = mean(YawErrorMODEL_Cmc(:,ko));
% 
% 
% 
% 
% 
%  meanPminus_save_NFmc(:,:,ko) = mean(Pminus_save_NFmc(:,:,:,ko));
%  
%  
%  %P_minus_C100Hz_save_NFmc(iM,:,:,:) = P_minus_C100Hz_save_NF;
% meanP_save_NFmc(:,:,ko) = mean(P_save_NFmc(:,:,:,ko));
% %  Q_save100Hz_NF = Q_save100Hz; 
% %  W_INS_Save_NF = W_INS_Save; 
% %  W_MODEL_Save_NF = W_MODEL_Save;  
% 
% meanx_hat_save_NFmc(:,ko) = mean(x_hat_save_NFmc(:,:,ko));
% meanV_C_save_NFmc(:,:,ko) = mean(V_C_save_NFmc(:,:,:,ko));
% meanz_save_NFmc(ko,:)  = mean(z_save_NFmc(:,ko,:));
% %meanN_save_NFmc(iM,ko)  =  mean(N_save_NFmc(iM,ko));
% 
%   
% %with fusion 
%  meanPminus_save_Cmc(:,:,ko) = mean(Pminus_save_Cmc(:,:,:,ko));
%  %P_minus_C100Hz_save_Cmc(iM,:,:,:)  = P_minus_C100Hz_save;
%  meanP_save_Cmc(:,:,ko)  = mean(P_save_Cmc(:,:,:,ko));
%  
% 
% meanx_hat_save_Cmc(:,ko) = mean(x_hat_save_Cmc(:,:,ko));
% meanV_C_save_Cmc(:,:,ko) = mean(V_C_save_Cmc(:,:,:,ko));
% meanz_save_Cmc(ko,:) = mean(z_save_Cmc(:,ko,:));
% %meanN_save_Cmc(iM,ko) = mean(N_save_Cmc(iM,ko));
% 
% 
% 
% for iM = 1:6
% 
% lambda_ss_outINS_H_NFmcProcessed(iM,ko) = max(lambda_ss_outINS_H_NFmc(iM,ko,:));
% lambda_ss_outINS_HmcProcessed(iM,ko) =  max(lambda_ss_outINS_Hmc(iM,ko,:));
% 
% 
% lambda_ss_outINS_V_NFmcProcessed(iM,ko) =  max(lambda_ss_outINS_V_NFmc(iM,ko,:));
% lambda_ss_outINS_VmcProcessed(iM,ko)=  max(lambda_ss_outINS_Vmc(iM,ko,:));
% 
% 
% lambda_ss_outMODEL_H_NFmcProcessed(iM,ko) =  max(lambda_ss_outMODEL_H_NFmc(iM,ko,:));
% lambda_ss_outMODEL_HmcProcessed(iM,ko)=  max(lambda_ss_outMODEL_Hmc(iM,ko,:));
% 
% lambda_ss_outMODEL_V_NFmcProcessed(iM,ko) =  max(lambda_ss_outMODEL_V_NFmc(iM,ko,:));
% lambda_ss_outMODEL_VmcProcessed(iM,ko) =  max(lambda_ss_outMODEL_Vmc(iM,ko,:));
% 
% lambda_ss_outGPS_H_NFmcProcessed(iM,ko) =  max(lambda_ss_outGPS_H_NFmc(iM,ko,:));
% lambda_ss_outGPS_HmcProcessed(iM,ko) =  max(lambda_ss_outGPS_Hmc(iM,ko,:));
% 
% lambda_ss_outGPS_V_NFmcProcessed(iM,ko) =  max(lambda_ss_outGPS_V_NFmc(iM,ko,:));
% lambda_ss_outGPS_VmcProcessed(iM,ko) =  max(lambda_ss_outGPS_Vmc(iM,ko,:));
% 
% 
% 
% end
% 
% 
% 
% 
% 
% end
% 
% %calculate variances
% 
% for ko = startprocessepoch:endprocessepoch
%     
% 
% %     for no = 1:Nmc
% %         
% %         N_PosErrorsINS_NFmc(no,km)
% %         
% %     end
% %     
% % end
%         
% %     
% %     for i =1:17
% %         bell(i) = (vect(i) - mean(vect))^2;
% %     end
% %     
% %     vart = (1/17)*sum(bell);
%     
%     
% %Note, need to calculate it this way instead of using matlabs var function
% %, similarly for calculating std deviation. Because the expected value of
% %these errors is zero. 
% 
% varN_PosErrorsINS_NFmc(ko) = (1/Nmc)*sum(N_PosErrorsINS_NFmc(:,ko).^2);
%  varE_PosErrorsINS_NFmc(ko) = (1/Nmc)*sum(E_PosErrorsINS_NFmc(:,ko).^2);   
%    varD_PosErrorsINS_NFmc(ko) = (1/Nmc)*sum(D_PosErrorsINS_NFmc(:,ko).^2); 
% 
% 
% 
%     
% %varN_PosErrorsINS_Cmc(ko) = var(N_PosErrorsINS_Cmc(:,ko));
% 
% 
% 
% varN_PosErrorsINS_Cmc(ko) = (1/Nmc)*sum(N_PosErrorsINS_Cmc(:,ko).^2);
% varE_PosErrorsINS_Cmc(ko) = (1/Nmc)*sum(E_PosErrorsINS_Cmc(:,ko).^2);
% varD_PosErrorsINS_Cmc(ko) = (1/Nmc)*sum(D_PosErrorsINS_Cmc(:,ko).^2);
% 
% 
% 
% 
% 
% varN_VelErrorsINS_NFmc(ko) = (1/Nmc)*sum(N_VelErrorsINS_NFmc(:,ko).^2);
%  varE_VelErrorsINS_NFmc(ko) = (1/Nmc)*sum(E_VelErrorsINS_NFmc(:,ko).^2);   
%    varD_VelErrorsINS_NFmc(ko) = (1/Nmc)*sum(D_VelErrorsINS_NFmc(:,ko).^2); 
%  
%    
%    
% 
% varN_VelErrorsINS_Cmc(ko) = (1/Nmc)*sum(N_VelErrorsINS_Cmc(:,ko).^2);
% varE_VelErrorsINS_Cmc(ko) = (1/Nmc)*sum(E_VelErrorsINS_Cmc(:,ko).^2);
% varD_VelErrorsINS_Cmc(ko) = (1/Nmc)*sum(D_VelErrorsINS_Cmc(:,ko).^2);
%     
% 
% 
% 
% 
% varRollErrorINS_NFmc(ko) = (1/Nmc)*sum(RollErrorINS_NFmc(:,ko).^2);
%  varPitchErrorINS_NFmc(ko) = (1/Nmc)*sum(PitchErrorINS_NFmc(:,ko).^2);   
%    varYawErrorINS_NFmc(ko) = (1/Nmc)*sum(YawErrorINS_NFmc(:,ko).^2); 
%  
%    
%    
% 
% varRollErrorINS_Cmc(ko) = (1/Nmc)*sum(RollErrorINS_Cmc(:,ko).^2);
% varPitchErrorINS_Cmc(ko) = (1/Nmc)*sum(PitchErrorINS_Cmc(:,ko).^2);
% varYawErrorINS_Cmc(ko) = (1/Nmc)*sum(YawErrorINS_Cmc(:,ko).^2);
%     
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %for ADM
% 
% 
% 
% 
% varN_PosErrorsMODEL_NFmc(ko) = (1/Nmc)*sum(N_PosErrorsMODEL_NFmc(:,ko).^2);
%  varE_PosErrorsMODEL_NFmc(ko) = (1/Nmc)*sum(E_PosErrorsMODEL_NFmc(:,ko).^2);   
%    varD_PosErrorsMODEL_NFmc(ko) = (1/Nmc)*sum(D_PosErrorsMODEL_NFmc(:,ko).^2); 
% 
% 
% 
%     
% %varN_PosErrorsMODEL_Cmc(ko) = var(N_PosErrorsMODEL_Cmc(:,ko));
% 
% 
% 
% varN_PosErrorsMODEL_Cmc(ko) = (1/Nmc)*sum(N_PosErrorsMODEL_Cmc(:,ko).^2);
% varE_PosErrorsMODEL_Cmc(ko) = (1/Nmc)*sum(E_PosErrorsMODEL_Cmc(:,ko).^2);
% varD_PosErrorsMODEL_Cmc(ko) = (1/Nmc)*sum(D_PosErrorsMODEL_Cmc(:,ko).^2);
% 
% 
% 
% 
% 
% varN_VelErrorsMODEL_NFmc(ko) = (1/Nmc)*sum(N_VelErrorsMODEL_NFmc(:,ko).^2);
%  varE_VelErrorsMODEL_NFmc(ko) = (1/Nmc)*sum(E_VelErrorsMODEL_NFmc(:,ko).^2);   
%    varD_VelErrorsMODEL_NFmc(ko) = (1/Nmc)*sum(D_VelErrorsMODEL_NFmc(:,ko).^2); 
%  
%    
%    
% 
% varN_VelErrorsMODEL_Cmc(ko) = (1/Nmc)*sum(N_VelErrorsMODEL_Cmc(:,ko).^2);
% varE_VelErrorsMODEL_Cmc(ko) = (1/Nmc)*sum(E_VelErrorsMODEL_Cmc(:,ko).^2);
% varD_VelErrorsMODEL_Cmc(ko) = (1/Nmc)*sum(D_VelErrorsMODEL_Cmc(:,ko).^2);
%     
% 
% 
% 
% 
% varRollErrorMODEL_NFmc(ko) = (1/Nmc)*sum(RollErrorMODEL_NFmc(:,ko).^2);
%  varPitchErrorMODEL_NFmc(ko) = (1/Nmc)*sum(PitchErrorMODEL_NFmc(:,ko).^2);   
%    varYawErrorMODEL_NFmc(ko) = (1/Nmc)*sum(YawErrorMODEL_NFmc(:,ko).^2); 
%  
%    
%    
% 
% varRollErrorMODEL_Cmc(ko) = (1/Nmc)*sum(RollErrorMODEL_Cmc(:,ko).^2);
% varPitchErrorMODEL_Cmc(ko) = (1/Nmc)*sum(PitchErrorMODEL_Cmc(:,ko).^2);
% varYawErrorMODEL_Cmc(ko) = (1/Nmc)*sum(YawErrorMODEL_Cmc(:,ko).^2);
%     
% 
% 
% 
% 
% 
% 
% %for GPS
% 
% 
% 
% 
% 
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% %get variances out for comparing with ensemble variance
% 
% for ko = startprocessepoch:endprocessepoch
%     
%     
%     [Rn, Re] = WGS84_calcRnRe(Lat_truth1Hz(ko));     
%     Rnh = Rn + Hgt_truth1Hz(ko);
%     Reh = Re + Hgt_truth1Hz(ko);
% 
%      
%     
% 
%     %need to convert to metres
% PnNF(ko) = (sqrt(meanP_save_NFmc(7,7,ko))*Rnh)^2;
% 
%     %need to convert to metres
% PeNF(ko) = (sqrt(meanP_save_NFmc(8,8,ko))*Reh)^2;
%     %need to convert to metres
% PdNF(ko) = (sqrt(meanP_save_NFmc(9,9,ko)))^2;
% 
% 
% %need to convert to metres
% PnC(ko) = (sqrt(meanP_save_Cmc(7,7,ko))*Rn)^2;
% 
%     %need to convert to metres
% PeC(ko) = (sqrt(meanP_save_Cmc(8,8,ko))*Re)^2;
%     %need to convert to metres
% PdC(ko) = (sqrt(meanP_save_Cmc(9,9,ko)))^2;
% 
% 
% 
% 
% Pn_minusNF(ko) = (sqrt(meanPminus_save_NFmc(7,7,ko))*Rn)^2;
% 
%     %need to convert to metres
% Pe_minusNF(ko) = (sqrt(meanPminus_save_NFmc(8,8,ko))*Re)^2;
%     %need to convert to metres
% Pd_minusNF(ko) = (sqrt(meanPminus_save_NFmc(9,9,ko)))^2;
% 
% 
% Pn_minusC(ko) = (sqrt(meanPminus_save_Cmc(7,7,ko))*Rn)^2;
% 
%     %need to convert to metres
% Pe_minusC(ko) = (sqrt(meanPminus_save_Cmc(8,8,ko))*Re)^2;
%     %need to convert to metres
% Pd_minusC(ko) = (sqrt(meanPminus_save_Cmc(9,9,ko)))^2;
% 
% 
% 
% 
% 
% 
% %FOR THE ADM STATE
% 
%     
% 
%     %need to convert to metres
% PnNF_MODEL(ko) = (sqrt(meanP_save_NFmc(24,24,ko))*Rnh)^2;
% 
%     %need to convert to metres
% PeNF_MODEL(ko) = (sqrt(meanP_save_NFmc(25,25,ko))*Reh)^2;
%     %need to convert to metres
% PdNF_MODEL(ko) = (sqrt(meanP_save_NFmc(25,26,ko)))^2;
% 
% 
% %need to convert to metres
% PnC_MODEL(ko) = (sqrt(meanP_save_Cmc(24,24,ko))*Rn)^2;
% 
%     %need to convert to metres
% PeC_MODEL(ko) = (sqrt(meanP_save_Cmc(25,25,ko))*Re)^2;
%     %need to convert to metres
% PdC_MODEL(ko) = (sqrt(meanP_save_Cmc(26,26,ko)))^2;
% 
% 
% 
% 
% Pn_minusNF_MODEL(ko) = (sqrt(meanPminus_save_NFmc(24,24,ko))*Rn)^2;
% 
%     %need to convert to metres
% Pe_minusNF_MODEL(ko) = (sqrt(meanPminus_save_NFmc(25,25,ko))*Re)^2;
%     %need to convert to metres
% Pd_minusNF_MODEL(ko) = (sqrt(meanPminus_save_NFmc(25,26,ko)))^2;
% 
% 
% Pn_minusC_MODEL(ko) = (sqrt(meanPminus_save_Cmc(24,24,ko))*Rn)^2;
% 
%     %need to convert to metres
% Pe_minusC_MODEL(ko) = (sqrt(meanPminus_save_Cmc(25,25,ko))*Re)^2;
%     %need to convert to metres
% Pd_minusC_MODEL(ko) = (sqrt(meanPminus_save_Cmc(26,26,ko)))^2;
% 
% 
% 
% end
% 
% 
% 
% 
% 
% %plot ensemble mean and variance
% figure;
% hold;
% 
% title 'Ensemble mean and variance unfused INS'
% 
% plot(N_PosErrorsINS_NFmc(:,startprocessepoch:endprocessepoch)','g')
% plot(2*sqrt(varN_PosErrorsINS_NFmc(startprocessepoch:endprocessepoch)),'r')
% plot(-2*sqrt(varN_PosErrorsINS_NFmc(startprocessepoch:endprocessepoch)),'r')
% plot(meanN_PosErrorsINS_NFmc(startprocessepoch:endprocessepoch))
% plot(2*sqrt(PnNF(startprocessepoch:endprocessepoch)),'b')
% plot(-2*sqrt(PnNF(startprocessepoch:endprocessepoch)),'b')
% plot(2*sqrt(Pn_minusNF(startprocessepoch:endprocessepoch)),'k')
% plot(-2*sqrt(Pn_minusNF(startprocessepoch:endprocessepoch)),'k')
% legend('green Pos Errors', 'red var pos errors', 'blue mean pos errors', 'blue posterior var', 'black apriori var');
% 
% 
% 
% figure;
% hold;
% title 'Ensemble mean and variance fused INS'
% plot(N_PosErrorsINS_Cmc(:,startprocessepoch:endprocessepoch)','g')
% plot(2*sqrt(varN_PosErrorsINS_Cmc(startprocessepoch:endprocessepoch)),'r')
% plot(-2*sqrt(varN_PosErrorsINS_Cmc(startprocessepoch:endprocessepoch)),'r')
% plot(meanN_PosErrorsINS_Cmc(startprocessepoch:endprocessepoch))
% plot(2*sqrt(PnC(startprocessepoch:endprocessepoch)),'b')
% plot(-2*sqrt(PnC(startprocessepoch:endprocessepoch)),'b')
% plot(2*sqrt(Pn_minusC(startprocessepoch:endprocessepoch)),'k')
% plot(-2*sqrt(Pn_minusC(startprocessepoch:endprocessepoch)),'k')
% 
% legend('green Pos Errors', 'red var pos errors', 'blue mean pos errors', 'blue posterior var', 'black apriori var');
% 
% 
% 
% plot(N_PosErrorsINS_NFmc(:,startprocessepoch:endprocessepoch)','y--')
% plot(2*sqrt(varN_PosErrorsINS_NFmc(startprocessepoch:endprocessepoch)),'r--')
% plot(-2*sqrt(varN_PosErrorsINS_NFmc(startprocessepoch:endprocessepoch)),'r--')
% plot(meanN_PosErrorsINS_NFmc(startprocessepoch:endprocessepoch),'b--')
% plot(2*sqrt(PnNF(startprocessepoch:endprocessepoch)),'b--')
% plot(-2*sqrt(PnNF(startprocessepoch:endprocessepoch)),'b--')
% plot(2*sqrt(Pn_minusNF(startprocessepoch:endprocessepoch)),'k--')
% plot(-2*sqrt(Pn_minusNF(startprocessepoch:endprocessepoch)),'k--')
% legend('green Pos Errors', 'red var pos errors', 'blue mean pos errors', 'blue posterior var', 'black apriori var');
% 
% 
% 
% 
% 
% tilefigs;
% 
% 
% 
% 
% %FOR THE ADM
% 
% %unfused
% %plot ensemble mean and variance
% figure;
% hold;
% title 'Ensemble mean and variance unfused ADM'
% plot(N_PosErrorsMODEL_NFmc(:,startprocessepoch:endprocessepoch)','g')
% plot(2*sqrt(varN_PosErrorsMODEL_NFmc(startprocessepoch:endprocessepoch)),'r')
% plot(-2*sqrt(varN_PosErrorsMODEL_NFmc(startprocessepoch:endprocessepoch)),'r')
% plot(meanN_PosErrorsMODEL_NFmc(startprocessepoch:endprocessepoch))
% plot(2*sqrt(PnNF_MODEL(startprocessepoch:endprocessepoch)),'b')
% plot(-2*sqrt(PnNF_MODEL(startprocessepoch:endprocessepoch)),'b')
% plot(2*sqrt(Pn_minusNF_MODEL(startprocessepoch:endprocessepoch)),'k')
% plot(-2*sqrt(Pn_minusNF_MODEL(startprocessepoch:endprocessepoch)),'k')
% legend('green Pos Errors', 'red var pos errors', 'blue mean pos errors', 'blue posterior var', 'black apriori var');
% 
% 
% %fused
% figure;
% hold;
% title 'Ensemble mean and variance fused ADM'
% plot(N_PosErrorsMODEL_Cmc(:,startprocessepoch:endprocessepoch)','g')
% plot(2*sqrt(varN_PosErrorsMODEL_Cmc(startprocessepoch:endprocessepoch)),'r')
% plot(-2*sqrt(varN_PosErrorsMODEL_Cmc(startprocessepoch:endprocessepoch)),'r')
% plot(meanN_PosErrorsMODEL_Cmc(startprocessepoch:endprocessepoch))
% plot(2*sqrt(PnC_MODEL(startprocessepoch:endprocessepoch)),'b')
% plot(-2*sqrt(PnC_MODEL(startprocessepoch:endprocessepoch)),'b')
% plot(2*sqrt(Pn_minusC_MODEL(startprocessepoch:endprocessepoch)),'k')
% plot(-2*sqrt(Pn_minusC_MODEL(startprocessepoch:endprocessepoch)),'k')
% legend('green Pos Errors', 'red var pos errors', 'blue mean pos errors', 'blue posterior var', 'black apriori var');
% 
% 
% 
% %this is for plotting on the same figure as the fused
% plot(N_PosErrorsMODEL_NFmc(:,startprocessepoch:endprocessepoch)','g')
% plot(2*sqrt(varN_PosErrorsMODEL_NFmc(startprocessepoch:endprocessepoch)),'r')
% plot(-2*sqrt(varN_PosErrorsMODEL_NFmc(startprocessepoch:endprocessepoch)),'r')
% plot(meanN_PosErrorsMODEL_NFmc(startprocessepoch:endprocessepoch))
% plot(2*sqrt(PnNF(startprocessepoch:endprocessepoch)),'b')
% plot(-2*sqrt(PnNF(startprocessepoch:endprocessepoch)),'b')
% plot(2*sqrt(Pn_minusNF(startprocessepoch:endprocessepoch)),'k')
% plot(-2*sqrt(Pn_minusNF(startprocessepoch:endprocessepoch)),'k')
% tilefigs;
% 
% 
% 
% 
% 
% 
% 
% 
% 
% % 
% % 
% % plot(N_PosErrorsINS_NFmc','b')
% % hold;
% % plot(N_PosErrorsINS_Cmc','g')
% % 
% % 
% %  plot(N_PosErrorsINS_Cmc,'r')
% %  plot(N_PosErrorsINS_Cmc','r') plot(N_PosErrorsINS_Cmc','k')
% % 
% % 
% % 
% % 
% %  plot(2*sqrt(PnNF),'r')
% %  plot(-2*sqrt(PnNF),'r')
% %  plot(-2*sqrt(PnC),'g')
% % plot(2*sqrt(PnC),'g')
% %  plot(2*sqrt(Pn_minusNF),'b')
% %  plot(-2*sqrt(Pn_minusNF),'b')
% %  plot(2*sqrt(Pn_minusC),'k')
% % plot(-2*sqrt(Pn_minusC),'k')
% % 
% 
% 
% 
% 
% 
% %state consistency test for No Fusion INS and GPS
% 
% %Nmc = 83;
% %see p 238 of Bar Shalom
% DOF = 3*Nmc;
% Qprob = 0.05;  %5 %
% 
% XlowNEES = chi2inv(Qprob/2,DOF)/Nmc;
% XhighNEES = chi2inv(1-Qprob/2,DOF)/Nmc;
% 
% 
% %state consistency test for No Fusion INS and GPS
% 
% XhighCountNEES_INS_NF = 0;
% XlowCountNEES_INS_NF = 0;
% 
% 
% for no = 1:Nmc
% 
%     %calculate normalised state errors. ie NEES Normalised  (state) Estimation Error Squared
%     for ko = startprocessepoch:endprocessepoch
% 
%         Ptemp(1:3,1:3) = P_save_NFmc(no,7:9,7:9,ko);
%         xhatsavetemp(1:3) = x_hat_save_NFmc(no,7:9,ko);
%         NEES_INS_NF(no,ko) = xhatsavetemp*inv(Ptemp)*xhatsavetemp';
% 
% 
%     end
% end
% 
% for ko = startprocessepoch:endprocessepoch
% 
%     meanNEES_INS_NF(ko) = mean(NEES_INS_NF(:,ko));
% 
%     %determine number of points outside bounds
%     if meanNEES_INS_NF(ko) > XhighNEES
%         XhighCountNEES_INS_NF = XhighCountNEES_INS_NF+1;
%     end
% 
%     if meanNEES_INS_NF(ko) < XlowNEES
%         XlowCountNEES_INS_NF = XlowCountNEES_INS_NF+1;
%     end
% end
% 
% 
% 
% %state consistency test for No Fusion INS and GPS
% 
% %Nmc = 83;
% %see p 238 of Bar Shalom
% DOF = 3*Nmc;
% Qprob = 0.05;  %5 %
% 
% XlowNEES = chi2inv(Qprob/2,DOF)/Nmc;
% XhighNEES = chi2inv(1-Qprob/2,DOF)/Nmc;
% 
% 
% 
% 
% 
% 
% 
% %state consistency test for WITH Fusion INS and GPS
% 
% XhighCountNEES_INS_C = 0;
% XlowCountNEES_INS_C = 0;
% 
% 
% for no = 1:Nmc
% 
%     %calculate normalised state errors. ie NEES Normalised  (state) Estimation Error Squared
%     for ko = startprocessepoch:endprocessepoch
% 
%         Ptemp(1:3,1:3) = P_save_Cmc(no,7:9,7:9,ko);
%         xhatsavetemp(1:3) = x_hat_save_Cmc(no,7:9,ko);
%         NEES_INS_C(no,ko) = xhatsavetemp*inv(Ptemp)*xhatsavetemp';
% 
% 
%     end
% end
% 
% for ko = startprocessepoch:endprocessepoch
% 
%     meanNEES_INS_C(ko) = mean(NEES_INS_C(:,ko));
% 
%     %determine number of points outside bounds
%     if meanNEES_INS_C(ko) > XhighNEES
%         XhighCountNEES_INS_C = XhighCountNEES_INS_C+1;
%     end
% 
%     if meanNEES_INS_C(ko) < XlowNEES
%         XlowCountNEES_INS_C = XlowCountNEES_INS_C+1;
%     end
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %state consistency test for No Fusion INS and GPS
% 
% %Nmc = 83;
% %see p 238 of Bar Shalom
% DOF = 9*Nmc;
% Qprob = 0.05;  %5 %
% 
% XlowNEES = chi2inv(Qprob/2,DOF)/Nmc;
% XhighNEES = chi2inv(1-Qprob/2,DOF)/Nmc;
% 
% 
% %state consistency test for No Fusion INS and GPS
% 
% XhighCountNEES_INS_NF = 0;
% XlowCountNEES_INS_NF = 0;
% 
% 
% for no = 1:Nmc
% 
%     %calculate normalised state errors. ie NEES Normalised  (state) Estimation Error Squared
%     for ko = startprocessepoch:endprocessepoch
% 
%         Ptemp(1:9,1:9) = P_save_NFmc(no,1:9,1:9,ko);
%         xhatsavetemp(1:9) = x_hat_save_NFmc(no,1:9,ko);
%         NEES_INS_NF(no,ko) = xhatsavetemp*inv(Ptemp)*xhatsavetemp';
% 
% 
%     end
% end
% 
% for ko = startprocessepoch:endprocessepoch
% 
%     meanNEES_INS_NF(ko) = mean(NEES_INS_NF(:,ko));
% 
%     %determine number of points outside bounds
%     if meanNEES_INS_NF(ko) > XhighNEES
%         XhighCountNEES_INS_NF = XhighCountNEES_INS_NF+1;
%     end
% 
%     if meanNEES_INS_NF(ko) < XlowNEES
%         XlowCountNEES_INS_NF = XlowCountNEES_INS_NF+1;
%     end
% end
% 
% 
% 
% 
%   
%   
% 
% %Plot NEES consistency test results
% 
% figure;
% hold;
% plot(meanNEES_INS_NF,'k-','LineWidth',2); title 'Normalized Estimation Error Squared No Fusion INS State' ; xlabel('Time (s)');  ylabel('NEES');
% plot(ones(1,length(meanNEES_INS_NF))*XhighNEES,'k--');
% plot(ones(1,length(meanNEES_INS_NF))*XlowNEES,'k--');
%   
%   
%   
%   
% %state consistency test with Fusion INS and GPS
% 
% %Nmc = 83;
% %see p 238 of Bar Shalom
% DOF = 9*Nmc;
% Qprob = 0.05;  %5 %
% 
% XlowNEES = chi2inv(Qprob/2,DOF)/Nmc;
% XhighNEES = chi2inv(1-Qprob/2,DOF)/Nmc;
% 
% 
% %state consistency test for No Fusion INS and GPS
% 
% XhighCountNEES_INS_C = 0;
% XlowCountNEES_INS_C = 0;
% 
% 
% for no = 1:Nmc
% 
%     %calculate normalised state errors. ie NEES Normalised  (state) Estimation Error Squared
%     for ko = startprocessepoch:endprocessepoch
% 
%         Ptemp(1:9,1:9) = P_save_Cmc(no,1:9,1:9,ko);
%         xhatsavetemp(1:9) = x_hat_save_Cmc(no,1:9,ko);
%         NEES_INS_C(no,ko) = xhatsavetemp*inv(Ptemp)*xhatsavetemp';
% 
% 
%     end
% end
% 
% for ko = startprocessepoch:endprocessepoch
% 
%     meanNEES_INS_C(ko) = mean(NEES_INS_C(:,ko));
% 
%     %determine number of points outside bounds
%     if meanNEES_INS_C(ko) > XhighNEES
%         XhighCountNEES_INS_C = XhighCountNEES_INS_C+1;
%     end
% 
%     if meanNEES_INS_C(ko) < XlowNEES
%         XlowCountNEES_INS_C = XlowCountNEES_INS_C+1;
%     end
% end
% 
% 
% 
% 
% %Plot NEES consistency test results
% 
% figure;
% hold;
% plot(meanNEES_INS_C,'k-','LineWidth',2); title 'Normalized Estimation Error Squared with Fusion INS State' ; xlabel('Time (s)');  ylabel('NEES');
% plot(ones(1,length(meanNEES_INS_C))*XhighNEES,'k--');
% plot(ones(1,length(meanNEES_INS_C))*XlowNEES,'k--');
%   
%   
% 
%   
%   
%   
% 
% 
% %state consistency test for No Fusion MODEL and GPS
% 
% XhighCountNEES_MODEL_NF = 0;
% XlowCountNEES_MODEL_NF = 0;
% %plot normalised state errors. ie NEES Normalised  (state) Estimation Error Squared
% for ko = startprocessepoch:endprocessepoch    
% NEES_MODEL_NF(ko) = meanx_hat_save_NFmc(24:26,ko)'*inv(meanP_save_NFmc(24:26,24:26,ko))*meanx_hat_save_NFmc(24:26,ko);
% 
% 
% %determine number of points outside bounds
% if NEES_MODEL_NF(ko) > XhighNEES
%     XhighCountNEES_MODEL_NF = XhighCountNEES_MODEL_NF+1;
% end
% 
% if NEES_MODEL_NF(ko) < XlowNEES
%     XlowCountNEES_MODEL_NF = XlowCountNEES_MODEL_NF+1;
% end
% end
% 
% 
% 
% 
% 
% 
% %state consistency test for with Fusion INS and GPS
% 
% XhighCountNEES_INS_C = 0;
% XlowCountNEES_INS_C = 0;
% %plot normalised state errors. ie NEES Normalised  (state) Estimation Error Squared
% for ko = startprocessepoch:endprocessepoch    
% NEES_INS_C(ko) = meanx_hat_save_Cmc(7:9,ko)'*inv(meanP_save_Cmc(7:9,7:9,ko))*meanx_hat_save_Cmc(7:9,ko);
% 
% 
% %determine number of points outside bounds
% if NEES_INS_C(ko) > XhighNEES
%     XhighCountNEES_INS_C = XhighCountNEES_INS_C+1;
% end
% 
% if NEES_INS_C(ko) < XlowNEES
%     XlowCountNEES_INS_C = XlowCountNEES_INS_C+1;
% end
% end
% 
% 
% 
% 
% %state consistency test for with Fusion MODEL and GPS
% 
% XhighCountNEES_MODEL_C = 0;
% XlowCountNEES_MODEL_C = 0;
% %plot normalised state errors. ie NEES Normalised  (state) Estimation Error Squared
% for ko = startprocessepoch:endprocessepoch    
% NEES_MODEL_C(ko) = meanx_hat_save_Cmc(24:26,ko)'*inv(meanP_save_Cmc(24:26,24:26,ko))*meanx_hat_save_Cmc(24:26,ko);
% 
% 
% %determine number of points outside bounds
% if NEES_MODEL_C(ko) > XhighNEES
%     XhighCountNEES_MODEL_C = XhighCountNEES_MODEL_C+1;
% end
% 
% if NEES_MODEL_C(ko) < XlowNEES
%     XlowCountNEES_MODEL_C = XlowCountNEES_MODEL_C+1;
% end
% end
% 
% 
% 
% %Plot NEES consistency test results
% 
% figure;
% hold;
% plot(NEES_INS_NF,'k-','LineWidth',2); title 'Normalized Estimation Error Squared No Fusion INS State' ; xlabel('Time (s)');  ylabel('NEES');
% plot(ones(1,length(NEES_INS_NF))*XhighNEES,'k--');
% plot(ones(1,length(NEES_INS_NF))*XlowNEES,'k--');
% 
% 
% figure;
% hold;
% plot(NEES_MODEL_NF,'k-','LineWidth',2); title 'Normalized Estimation Error Squared No Fusion ADM State' ; xlabel('Time (s)');  ylabel('NEES');
% plot(ones(1,length(NEES_MODEL_NF))*XhighNEES,'k--');
% plot(ones(1,length(NEES_MODEL_NF))*XlowNEES,'k--');
% 
% 
% figure;
% hold;
% plot(NEES_INS_C,'k-','LineWidth',2); title 'Normalized Estimation Error Squared with Fusion INS State' ; xlabel('Time (s)');  ylabel('NEES');
% plot(ones(1,length(NEES_INS_C))*XhighNEES,'k--');
% plot(ones(1,length(NEES_INS_C))*XlowNEES,'k--');
% 
% 
% figure;
% hold;
% plot(NEES_MODEL_C,'k-','LineWidth',2); title 'Normalized Estimation Error Squared with Fusion ADM State' ; xlabel('Time (s)');  ylabel('NEES');
% plot(ones(1,length(NEES_MODEL_C))*XhighNEES,'k--');
% plot(ones(1,length(NEES_MODEL_C))*XlowNEES,'k--');
% 
% 
% tilefigs;
% pause;
% close all;
% 
% 
% 
% 
% 
% %Innovations Test
% 
% 
% %Whiteness Test
% 
% 
% %calculate % improvements based the means:
% 
% HPLImprovement = meanHPL_INS_NFmc(startprocessepoch:endprocessepoch) -meanHPL_INS_Cmc(startprocessepoch:endprocessepoch);
% VPLImprovement = meanVPL_INS_NFmc(startprocessepoch:endprocessepoch) -meanVPL_INS_Cmc(startprocessepoch:endprocessepoch);
% 
% 
% HPLImprovepercent = (HPLImprovement./meanHPL_INS_Cmc(startprocessepoch:endprocessepoch))*100;
% VPLImprovepercent = (VPLImprovement./meanVPL_INS_Cmc(startprocessepoch:endprocessepoch))*100;
% 
% figure;
% plot(HPLImprovepercent(2:length(HPLImprovepercent)));
% 
% figure;
% plot(VPLImprovepercent(2:length(HPLImprovepercent)));
% % 
% % meanHPLImprov = mean(HPLImprovepercent(2:length(HPLImprovepercent)))
% % meanVPLImprov = mean(VPLImprovepercent(2:length(VPLImprovepercent)))
% 
% 
% meanHPLImprov = mean(HPLImprovepercent(2:length(HPLImprovepercent)))
% meanVPLImprov = mean(VPLImprovepercent(2:length(HPLImprovepercent)))
% 
% 
% 
% 
% 
% figure();
% hold;
% plot(HAL(startprocessepoch:endprocessepoch), 'k','LineWidth',2);
% %plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
% plot(meanHPL_INS_Cmc(startprocessepoch:endprocessepoch),'r--','LineWidth',2);  title 'Horizontal Position Error' ; xlabel('Time (s)');  ylabel('H Dist (m)');
% plot(meanHPL_INS_NFmc(startprocessepoch:endprocessepoch),'g--','LineWidth',2);  
% %plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
% plot(meanHPE_INS_Cmc(startprocessepoch:endprocessepoch),'g','LineWidth',2);    %navigation system error
% plot(meanHPE_GPS_mc(startprocessepoch:endprocessepoch),'k--','LineWidth',1);
% 
% legend('HAL','HPL INS-ADM','HPL INS', 'HNSE Fused', 'HNSE GPS');       
% 
% axis([0 endprocessepoch-startprocessepoch 0 60]);
% 
% %axis([0 endepochplot-startepochplot 0 1000]);
% 
% 
% 
% 
% 
% %For Vertical HPL etc
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(VAL(startprocessepoch:endprocessepoch), 'k','LineWidth',2);
% %plot(AccuracyLineV(startepochplot:endepochplot), 'k','LineWidth',1);
% plot(meanVPL_INS_Cmc(startprocessepoch:endprocessepoch),'r--','LineWidth',2);  title 'Vertical Position Error' ; xlabel('Time (s)');  ylabel('V Dist (m)');
% plot(meanVPL_INS_NFmc(startprocessepoch:endprocessepoch),'g--','LineWidth',2);  
% %plot(VPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
% plot(meanVPE_INS_Cmc(startprocessepoch:endprocessepoch),'g','LineWidth',2);    %navigation system error
% plot(meanVPE_GPS_mc(startprocessepoch:endprocessepoch),'k--','LineWidth',1);
% 
% legend('VAL','VPL INS-ADM','VPL INS','VNSE Fused','VNSE GPS');    
% 
% axis([0 endprocessepoch-startprocessepoch 0 60]);
% 
% 
% 
% 
% 
% 
% %Consistency Checks of the filter
% 
% 
% 
% 
% %state errors test
% 
% 
% %whiteness test
% 
% 
% 
% %innovations test
% 
% 
% 
% 
% Pn_minusNF(ko) = (sqrt(meanPminus_save_NFmc(7,7,ko))*Rn)^2;
% 
%     %need to convert to metres
% Pe_minusNF(ko) = (sqrt(meanPminus_save_NFmc(8,8,ko))*Re)^2;
%     %need to convert to metres
% Pd_minusNF(ko) = (sqrt(meanPminus_save_NFmc(9,9,ko)))^2;
% 
% 
% Pn_minusC(ko) = (sqrt(meanPminus_save_Cmc(7,7,ko))*Rn)^2;
% 
%     %need to convert to metres
% Pe_minusC(ko) = (sqrt(meanPminus_save_Cmc(8,8,ko))*Re)^2;
%     %need to convert to metres
% Pd_minusC(ko) = (sqrt(meanPminus_save_Cmc(9,9,ko)))^2;
% 
% 
% 
% 
% %===========================================
% %Do a determinant check of the process models covariance Pminus
% %===========================================
% 
% 
% for ko = startprocessepoch:endprocessepoch
% 
% 
% detPosPminus_save_INS_NFmc(ko) = det(meanPminus_save_NFmc(7:9,7:9,ko)); 
% detPosPminus_save_INS_Cmc(ko) = det(meanPminus_save_Cmc(7:9,7:9,ko)); 
% 
% detPosPminus_save_MODEL_NFmc(ko) = det(meanPminus_save_NFmc(24:26,24:26,ko)); 
% detPosPminus_save_MODEL_Cmc(ko) = det(meanPminus_save_Cmc(24:26,24:26,ko)); 
% 
% 
% detVelPminus_save_INS_NFmc(ko) = det(meanPminus_save_NFmc(4:6,4:6,ko)); 
% detVelPminus_save_INS_Cmc(ko) = det(meanPminus_save_Cmc(4:6,4:6,ko)); 
% 
% detVelPminus_save_MODEL_NFmc(ko)  = det(meanPminus_save_NFmc(21:23,21:23,ko)); 
% detVelPminus_save_MODEL_Cmc(ko) = det(meanPminus_save_Cmc(21:23,21:23,ko)); 
% 
% 
% 
% detAttPminus_save_INS_NFmc(ko) = det(meanPminus_save_NFmc(1:3,1:3,ko)); 
% detAttPminus_save_INS_Cmc(ko) = det(meanPminus_save_Cmc(1:3,1:3,ko)); 
% 
% detAttPminus_save_MODEL_NFmc(ko) = det(meanPminus_save_NFmc(18:20,18:20,ko)); 
% detAttPminus_save_MODEL_Cmc(ko) = det(meanPminus_save_Cmc(18:20,18:20,ko)); 
% 
% 
% 
% 
% %do determinant check of the update covariance
% 
% detPosP_save_INS_NFmc(ko) = det(meanP_save_NFmc(7:9,7:9,ko)); 
% detPosP_save_INS_Cmc(ko) = det(meanP_save_Cmc(7:9,7:9,ko)); 
% 
% 
% detPosP_save_MODEL_NFmc(ko) = det(meanP_save_NFmc(24:26,24:26,ko)); 
% detPosP_save_MODEL_Cmc(ko) = det(meanP_save_Cmc(24:26,24:26,ko)); 
% 
% 
% detVelP_save_INS_NFmc(ko) = det(meanP_save_NFmc(4:6,4:6,ko)); 
% detVelP_save_INS_Cmc(ko) = det(meanP_save_Cmc(4:6,4:6,ko)); 
% 
% detVelP_save_MODEL_NFmc(ko)  = det(meanP_save_NFmc(21:23,21:23,ko)); 
% detVelP_save_MODEL_Cmc(ko) = det(meanP_save_Cmc(21:23,21:23,ko)); 
% 
% 
% detAttP_save_INS_NFmc(ko) = det(meanP_save_NFmc(1:3,1:3,ko)); 
% detAttP_save_INS_Cmc(ko) = det(meanP_save_Cmc(1:3,1:3,ko)); 
% 
% detAttP_save_MODEL_NFmc(ko) = det(meanP_save_NFmc(18:20,18:20,ko)); 
% detAttP_save_MODEL_Cmc(ko) = det(meanP_save_Cmc(18:20,18:20,ko)); 
% 
% 
% 
% end
% 
% 
% 
% %plot determinant results for P minus
% 
% %Positions
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detPosPminus_save_INS_Cmc(startprocessepoch:endprocessepoch),'r','LineWidth',2);  title 'P minus Determinant Position, INS State' ; xlabel('Time (s)'); % ylabel('V Dist (m)');
% plot(detPosPminus_save_INS_NFmc(startprocessepoch:endprocessepoch),'g','LineWidth',2);  
% legend('INS Fused','INS Not Fused');
% 
% 
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detPosPminus_save_MODEL_Cmc(startprocessepoch:endprocessepoch),'k','LineWidth',2);   title 'P minus Determinant Position, ADM State' ; xlabel('Time (s)'); % ylabel('V Dist (m)');
% plot(detPosPminus_save_MODEL_NFmc(startprocessepoch:endprocessepoch),'b','LineWidth',2);  
% legend('ADM Fused','ADM Not Fused');
% 
% 
% 
% 
% 
% %Velocities
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detVelPminus_save_INS_Cmc(startprocessepoch:endprocessepoch),'r','LineWidth',2);  title 'P minus Determinant Position, INS State' ; xlabel('Time (s)'); % ylabel('V Dist (m)');
% plot(detVelPminus_save_INS_NFmc(startprocessepoch:endprocessepoch),'g','LineWidth',2);  
% legend('INS Fused','INS Not Fused');
% 
% 
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detVelPminus_save_MODEL_Cmc(startprocessepoch:endprocessepoch),'k','LineWidth',2);   title 'P minus Determinant Position, ADM State' ; xlabel('Time (s)');  %ylabel('V Dist (m)');
% plot(detVelPminus_save_MODEL_NFmc(startprocessepoch:endprocessepoch),'b','LineWidth',2);  
% legend('ADM Fused','ADM Not Fused');
% 
% 
% 
% %Attitude 
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detAttPminus_save_INS_Cmc(startprocessepoch:endprocessepoch),'r','LineWidth',2);  title 'P minus Determinant Position, INS State' ; xlabel('Time (s)'); % ylabel('V Dist (m)');
% plot(detAttPminus_save_INS_NFmc(startprocessepoch:endprocessepoch),'g','LineWidth',2);  
% legend('INS Fused','INS Not Fused');
% 
% 
% 
% 
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detAttPminus_save_MODEL_Cmc(startprocessepoch:endprocessepoch),'k','LineWidth',2);   title 'P minus Determinant Position, ADM State' ; xlabel('Time (s)'); % ylabel('V Dist (m)');
% plot(detAttPminus_save_MODEL_NFmc(startprocessepoch:endprocessepoch),'b','LineWidth',2);  
% legend('ADM Fused','ADM Not Fused');
% 
% 
% 
% 
% 
% 
% 
% 
% %plot determinant results for P update ie P posteriori
% 
% 
% %Positions
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detPosP_save_INS_Cmc(startprocessepoch:endprocessepoch),'r','LineWidth',2);  title 'P Determinant Position, INS State' ; xlabel('Time (s)'); % ylabel('V Dist (m)');
% plot(detPosP_save_INS_NFmc(startprocessepoch:endprocessepoch),'g','LineWidth',2);  
% legend('INS Fused','INS Not Fused');
% 
% 
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detPosP_save_MODEL_Cmc(startprocessepoch:endprocessepoch),'k','LineWidth',2);   title 'P Determinant Position, ADM State' ; xlabel('Time (s)'); % ylabel('V Dist (m)');
% plot(detPosP_save_MODEL_NFmc(startprocessepoch:endprocessepoch),'b','LineWidth',2);  
% legend('ADM Fused','ADM Not Fused');
% 
% 
% 
% 
% 
% %Velocities
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detVelP_save_INS_Cmc(startprocessepoch:endprocessepoch),'r','LineWidth',2);  title 'P Determinant Position, INS State' ; xlabel('Time (s)');  %ylabel('V Dist (m)');
% plot(detVelP_save_INS_NFmc(startprocessepoch:endprocessepoch),'g','LineWidth',2);  
% legend('INS Fused','INS Not Fused');
% 
% 
% 
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detVelP_save_MODEL_Cmc(startprocessepoch:endprocessepoch),'k','LineWidth',2);   title 'P Determinant Position, ADM State' ; xlabel('Time (s)');  %ylabel('V Dist (m)');
% plot(detVelP_save_MODEL_NFmc(startprocessepoch:endprocessepoch),'b','LineWidth',2);  
% legend('ADM Fused','ADM Not Fused');
% 
% 
% 
% %Attitude 
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detAttP_save_INS_Cmc(startprocessepoch:endprocessepoch),'r','LineWidth',2);  title 'P Determinant Position, INS State' ; xlabel('Time (s)'); % ylabel('V Dist (m)');
% plot(detAttP_save_INS_NFmc(startprocessepoch:endprocessepoch),'g','LineWidth',2);  
% 
% legend('INS Fused','INS Not Fused');
% 
% 
% figure();
% %plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
% hold;
% 
% plot(detAttP_save_MODEL_Cmc(startprocessepoch:endprocessepoch),'k','LineWidth',2);   title 'P Determinant Position, ADM State' ; xlabel('Time (s)');  %ylabel('V Dist (m)');
% plot(detAttP_save_MODEL_NFmc(startprocessepoch:endprocessepoch),'b','LineWidth',2);  
% 
% legend('ADM Fused','ADM Not Fused');
% 
% 
% 
%  tilefigs;
%  
%  
% 
% %Is esnemble mean estimation error unbiased (ie has effectively zero mean),
% %and is the ensemble covariance in close agreement with the theoreteical
% %valuce computed from the KF (see p 273 of Grewal , kalman filtering,
% %theory and practice
% 
% 
% 
% 
% %The INS process model without fusion
% 
% % 
% % meanN_PosErrorsINS_NFmc
% % 
% % 
% % 
% % meanE_PosErrorsINS_NFmc
% % 
% % 
% % meanD_PosErrorsINS_NFmc
% % 
% % 
% % 
% % meanN_VelErrorsINS_NFmc
% % 
% % 
% % 
% % meanE_VelErrorsINS_NFmc
% % 
% % 
% % meanD_VelErrorsINS_NFmc
% 
% 
% 
% 
% 
% %The INS process model with fusion
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %The ADM process model without fusion
% 
% meanN_PosErrorsINS_NFmc
% 
% 
% 
% 
% %The INS process model with fusion
% 
% 
% 
% 
% 
% 
% 
% %The ADM process model with fusion
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% figure;
% hold;
% 
% for hh = 1:33
%     plot(HPLNmcNF(hh,:))
%     
%      plot(HPLNmc(hh,:),'r')
% end
% 
% 
% figure;
% hold;
% for hh = 1:32
%     plot(HPENmcNF(hh,:))
%     plot(HPENmc(hh,:),'r')
% end
% 
% figure;
% hold;
% for hh = 1:50
%     plot(VPENmcNF(hh,:))
%     plot(VPENmc(hh,:),'r')
% end
% 
% 
% tend = cputime;
% 
% 
% 
% figure;
% hold;
% for hh = 1:50
%     plot(HPENmcNF(hh,:))
%     plot(HPENmc(hh,:),'r')
% end
% 
% 
% 
% figure;
% hold;
% 
% for hh = 1:50
%     plot(VPLNmcNF(hh,:))
%     
%      plot(VPLNmc(hh,:),'r')
% end
% 
% 
% %process Montecarlo run
% 
% 
% %Consistency Check
% 
% % 
% % for iM = 1:Nmc
% %     
% %     for i = startepoch:endepoch
% %         
% %         
% %     end
% %     
% %     
% % end
% 
% 
% 
% %Process Accuracy ie RMS 
% 
% %Position
% 
% %velocity
% 
% 
% %attitude
% 
% 
% 
% 
% 
% %Do plots
% 
% 
% 
% 
% %for estimating Qk1k2
% 
% 
% %For FUSED 
% 
% 
% for iM = 1:Nmc
% 
% for i = startepoch:endepoch
%     
% 
%  %position errors
%     
% LatErrorINS(i) = Lat_truth1Hz(i) - x_state_total_minusINSsaveCmc(iM,8,i);
% LonErrorINS(i) = Lon_truth1Hz(i)  - x_state_total_minusINSsaveCmc(iM,9,i);
% HgtErrorINS(i) = Hgt_truth1Hz(i) - x_state_total_minusINSsaveCmc(iM,10,i);
% 
%     
% LatErrorMODEL(i) = Lat_truth1Hz(i) - x_state_total_minusMODELsaveCmc(iM,8,i);
% LonErrorMODEL(i) = Lon_truth1Hz(i)  - x_state_total_minusMODELsaveCmc(iM,9,i);
% HgtErrorMODEL(i) = Hgt_truth1Hz(i) - x_state_total_minusMODELsaveCmc(iM,10,i);
% 
% 
% 
% %velocity errors    
% N_VelErrorINS(i) = V_n_truth1Hz(i) - x_state_total_minusINSsaveCmc(iM,4,i);
% E_VelErrorINS(i) = V_e_truth1Hz(i)  - x_state_total_minusINSsaveCmc(iM,5,i);
% D_VelErrorINS(i) = V_d_truth1Hz(i) - x_state_total_minusINSsaveCmc(iM,6,i);
% 
%     
% N_VelErrorMODEL(i) = V_n_truth1Hz(i) - x_state_total_minusMODELsaveCmc(iM,4,i);
% E_VelErrorMODEL(i) = V_e_truth1Hz(i)  - x_state_total_minusMODELsaveCmc(iM,5,i);
% D_VelErrorMODEL(i) = V_d_truth1Hz(i) - x_state_total_minusMODELsaveCmc(iM,6,i);
% 
% 
% 
% 
% %Attitude tilt errors (for Q matrix
% 
% 
% 
% 
% 
% 
%  q0 = Quaternions_truth1Hz(1,i);
%     q1 = Quaternions_truth1Hz(2,i);
%     q2 = Quaternions_truth1Hz(3,i);
%     q3 = Quaternions_truth1Hz(4,i);
% 
%     quat = [q0,q1,q2,q3];
% 
%     [euler] = QuatToEuler(quat);
% 
%     phitoconvTRUTH = euler(1);
%     thetatoconvTRUTH = euler(2);
%     psitoconvTRUTH = euler(3);
%     
%     
%     
%     
%     
%     q0 = x_state_total_minusINSsaveCmc(iM,1,i);
%     q1 = x_state_total_minusINSsaveCmc(iM,2,i);
%     q2 = x_state_total_minusINSsaveCmc(iM,3,i);
%     q3 = x_state_total_minusINSsaveCmc(iM,4,i);
% 
%     quat = [q0,q1,q2,q3];
% 
%     [euler] = QuatToEuler(quat);
% 
%     phitoconvINS = euler(1);
%     thetatoconvINS = euler(2);
%     psitoconvINS = euler(3);    
%     
%     
%     q0 = x_state_total_minusMODELsaveCmc(iM,1,i);
%     q1 = x_state_total_minusMODELsaveCmc(iM,2,i);
%     q2 = x_state_total_minusMODELsaveCmc(iM,3,i);
%     q3 = x_state_total_minusMODELsaveCmc(iM,4,i);
% 
%     quat = [q0,q1,q2,q3];
% 
%     [euler] = QuatToEuler(quat);
% 
%     phitoconvMODEL = euler(1);
%     thetatoconvMODEL = euler(2);
%     psitoconvMODEL = euler(3);    
%     
%     
%        
%     
%     EulerDiffINS = [phitoconvTRUTH,thetatoconvTRUTH, psitoconvTRUTH] - [phitoconvINS,thetatoconvINS, psitoconvINS];
%     
%       
%     
%     
%     
%     
%     EulerDiffMODEL = [phitoconvTRUTH,thetatoconvTRUTH, psitoconvTRUTH] - [phitoconvMODEL,thetatoconvMODEL, psitoconvMODEL];
%                
%     
%     
%     T_alpha_TilttoEuler = [-cos(psitoconvTRUTH)/cos(thetatoconvTRUTH), -sin(psitoconvTRUTH)/cos(thetatoconvTRUTH), 0;        
%                                 sin(psitoconvTRUTH), -cos(psitoconvTRUTH), 0;
%                                  -cos(psitoconvTRUTH)*tan(thetatoconvTRUTH), -sin(psitoconvTRUTH)*tan(thetatoconvTRUTH), -1;];    
%                              
%        
% 
%     TiltErrorINS = T_alpha_TilttoEuler'*[ EulerDiffINS]';    
% 
% 
%     TiltErrorMODEL = T_alpha_TilttoEuler'*[ EulerDiffMODEL]';    
%     
%     
% 
% % 
% % 
% % 
% %     
% %      x_state_total_minusINS(:,i) = [Quaternions_truth1Hz(1,i),   %q0
% %             Quaternions_truth1Hz(2,i),    %q1
% %             Quaternions_truth1Hz(3,i),     %q2
% %            Quaternions_truth1Hz(4,i),     %q3
% %             V_n_truth1Hz(i),    %Vn
% %             V_e_truth1Hz(i),    %Ve
% %             V_d_truth1Hz(i),   %Vd
% %             Lat_truth1Hz(i),  %lat
% %             Lon_truth1Hz(i), %lon
% %             Hgt_truth1Hz(i)];   
% %         
% %             
% %             
%             
%                  
%             
%    x_state_total_minusINSerrorVectorCmc(iM,1:9,i) =  [TiltErrorINS', N_VelErrorINS(i),  E_VelErrorINS(i),  D_VelErrorINS(i), LatErrorINS(i),LonErrorINS(i),  HgtErrorINS(i)]';
%                         
%    x_state_total_minusMODELerrorVectorCmc(iM,1:9,i) =  [TiltErrorMODEL', N_VelErrorMODEL(i),  E_VelErrorMODEL(i),  D_VelErrorMODEL(i), LatErrorMODEL(i),LonErrorMODEL(i),  HgtErrorMODEL(i)]';
%       
%        
%     
% end
% 
% end
% 
% 
% 
% 
% 
% %calculate correlations and variances matrix between INS process and ADM
% %process
% for iM = 1:Nmc
% 
% for i = startepoch:endepoch
%     
%     
%     for j = 1:9
%         for k = 1:9
%             
%             
%     CovMatrixTemp(iM,j,k,i) = (x_state_total_minusINSerrorVectorCmc(iM,j,i)*x_state_total_minusMODELerrorVectorCmc(iM,k,i)');
%    
%     
%         end
%     end
%     
%     
% end
% 
% end
% 
% 
% %now get the mean 
% %for iM = 1:Nmc
% 
% for i = startepoch:endepoch   
%     
%     for j = 1:9
%         for k = 1:9            
%             
%     CovMatrix(j,k,i) = (1/Nmc)*sum(CovMatrixTemp(:,j,k,i));   
%     
%         end
%     end    
%     
% end
% 
% %end
% 
% 
% 
% 
% 
% 
% %Calculate variances of the errors
% for ko = startepoch:endepoch
% for j = 1:9
% 
%         
% 
% var_x_state_total_minusINSerrorVectorCmc(j,ko) = (1/Nmc)*sum(x_state_total_minusINSerrorVectorCmc(:,j,ko).^2);
% 
% 
% var_x_state_total_minusMODELerrorVectorCmc(j,ko) = (1/Nmc)*sum(x_state_total_minusMODELerrorVectorCmc(:,j,ko).^2);
%     
%     
% 
% end
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% for i = startepoch:endepoch
%     
%     
%     for j = 1:9
%         for k = 1:9
%             
%             
%    % CovMatrixCorrCoef(j,k,i) = CovMatrix(j,k,i)/(sqrt(CovMatrix(j,j,i))*sqrt(CovMatrix(k,k,i)));
%    
%                 
%     CovMatrixCorrCoef(j,k,i) = CovMatrix(j,k,i)/(sqrt(var_x_state_total_minusINSerrorVectorCmc(j,i))*sqrt(var_x_state_total_minusMODELerrorVectorCmc(k,i)));
%    
%     
%         end
%     end
%     
%     
% end
% 
% 
% 
% 
% 
% 
% %For NOT FUSED 
% 
% 
% for iM = 1:Nmc
% 
% for i = startepoch:endepoch
%     
% 
%  %position errors
%     
% LatErrorINS(i) = Lat_truth1Hz(i) - x_state_total_minusINSsaveNFmc(iM,8,i);
% LonErrorINS(i) = Lon_truth1Hz(i)  - x_state_total_minusINSsaveNFmc(iM,9,i);
% HgtErrorINS(i) = Hgt_truth1Hz(i) - x_state_total_minusINSsaveNFmc(iM,10,i);
% 
%     
% LatErrorMODEL(i) = Lat_truth1Hz(i) - x_state_total_minusMODELsaveNFmc(iM,8,i);
% LonErrorMODEL(i) = Lon_truth1Hz(i)  - x_state_total_minusMODELsaveNFmc(iM,9,i);
% HgtErrorMODEL(i) = Hgt_truth1Hz(i) - x_state_total_minusMODELsaveNFmc(iM,10,i);
% 
% 
% 
% %velocity errors    
% N_VelErrorINS(i) = V_n_truth1Hz(i) - x_state_total_minusINSsaveNFmc(iM,4,i);
% E_VelErrorINS(i) = V_e_truth1Hz(i)  - x_state_total_minusINSsaveNFmc(iM,5,i);
% D_VelErrorINS(i) = V_d_truth1Hz(i) - x_state_total_minusINSsaveNFmc(iM,6,i);
% 
%     
% N_VelErrorMODEL(i) = V_n_truth1Hz(i) - x_state_total_minusMODELsaveNFmc(iM,4,i);
% E_VelErrorMODEL(i) = V_e_truth1Hz(i)  - x_state_total_minusMODELsaveNFmc(iM,5,i);
% D_VelErrorMODEL(i) = V_d_truth1Hz(i) - x_state_total_minusMODELsaveNFmc(iM,6,i);
% 
% 
% 
% 
% %Attitude tilt errors (for Q matrix
% 
% 
% 
% 
% 
% 
%  q0 = Quaternions_truth1Hz(1,i);
%     q1 = Quaternions_truth1Hz(2,i);
%     q2 = Quaternions_truth1Hz(3,i);
%     q3 = Quaternions_truth1Hz(4,i);
% 
%     quat = [q0,q1,q2,q3];
% 
%     [euler] = QuatToEuler(quat);
% 
%     phitoconvTRUTH = euler(1);
%     thetatoconvTRUTH = euler(2);
%     psitoconvTRUTH = euler(3);
%     
%     
%     
%     
%     
%     q0 = x_state_total_minusINSsaveNFmc(iM,1,i);
%     q1 = x_state_total_minusINSsaveNFmc(iM,2,i);
%     q2 = x_state_total_minusINSsaveNFmc(iM,3,i);
%     q3 = x_state_total_minusINSsaveNFmc(iM,4,i);
% 
%     quat = [q0,q1,q2,q3];
% 
%     [euler] = QuatToEuler(quat);
% 
%     phitoconvINS = euler(1);
%     thetatoconvINS = euler(2);
%     psitoconvINS = euler(3);    
%     
%     
%     q0 = x_state_total_minusMODELsaveNFmc(iM,1,i);
%     q1 = x_state_total_minusMODELsaveNFmc(iM,2,i);
%     q2 = x_state_total_minusMODELsaveNFmc(iM,3,i);
%     q3 = x_state_total_minusMODELsaveNFmc(iM,4,i);
% 
%     quat = [q0,q1,q2,q3];
% 
%     [euler] = QuatToEuler(quat);
% 
%     phitoconvMODEL = euler(1);
%     thetatoconvMODEL = euler(2);
%     psitoconvMODEL = euler(3);    
%     
%     
%        
%     
%     EulerDiffINS = [phitoconvTRUTH,thetatoconvTRUTH, psitoconvTRUTH] - [phitoconvINS,thetatoconvINS, psitoconvINS];
%     
%       
%     
%     
%     
%     
%     EulerDiffMODEL = [phitoconvTRUTH,thetatoconvTRUTH, psitoconvTRUTH] - [phitoconvMODEL,thetatoconvMODEL, psitoconvMODEL];
%                
%     
%     
%     T_alpha_TilttoEuler = [-cos(psitoconvTRUTH)/cos(thetatoconvTRUTH), -sin(psitoconvTRUTH)/cos(thetatoconvTRUTH), 0;        
%                                 sin(psitoconvTRUTH), -cos(psitoconvTRUTH), 0;
%                                  -cos(psitoconvTRUTH)*tan(thetatoconvTRUTH), -sin(psitoconvTRUTH)*tan(thetatoconvTRUTH), -1;];    
%                              
%        
% 
%     TiltErrorINS = T_alpha_TilttoEuler'*[ EulerDiffINS]';    
% 
% 
%     TiltErrorMODEL = T_alpha_TilttoEuler'*[ EulerDiffMODEL]';    
%     
%     
% 
% % 
% % 
% % 
% %     
% %      x_state_total_minusINS(:,i) = [Quaternions_truth1Hz(1,i),   %q0
% %             Quaternions_truth1Hz(2,i),    %q1
% %             Quaternions_truth1Hz(3,i),     %q2
% %            Quaternions_truth1Hz(4,i),     %q3
% %             V_n_truth1Hz(i),    %Vn
% %             V_e_truth1Hz(i),    %Ve
% %             V_d_truth1Hz(i),   %Vd
% %             Lat_truth1Hz(i),  %lat
% %             Lon_truth1Hz(i), %lon
% %             Hgt_truth1Hz(i)];   
% %         
% %             
% %             
%             
%                  
%             
%    x_state_total_minusINSerrorVectorNFmc(iM,1:9,i) =  [TiltErrorINS', N_VelErrorINS(i),  E_VelErrorINS(i),  D_VelErrorINS(i), LatErrorINS(i),LonErrorINS(i),  HgtErrorINS(i)]';
%                         
%    x_state_total_minusMODELerrorVectorNFmc(iM,1:9,i) =  [TiltErrorMODEL', N_VelErrorMODEL(i),  E_VelErrorMODEL(i),  D_VelErrorMODEL(i), LatErrorMODEL(i),LonErrorMODEL(i),  HgtErrorMODEL(i)]';
%       
%        
%     
% end
% 
% end
% 
% 
% 
% 
% 
% %calculate correlations and variances matrix between INS process and ADM
% %process
% for iM = 1:Nmc
% 
% for i = startepoch:endepoch
%     
%     
%     for j = 1:9
%         for k = 1:9
%             
%             
%     CovMatrixTempNF(iM,j,k,i) = (x_state_total_minusINSerrorVectorNFmc(iM,j,i)*x_state_total_minusMODELerrorVectorNFmc(iM,k,i)');
%    
%     
%         end
%     end
%     
%     
% end
% 
% end
% 
% 
% %now get the mean 
% %for iM = 1:Nmc
% 
% for i = startepoch:endepoch   
%     
%     for j = 1:9
%         for k = 1:9            
%             
%     CovMatrixNF(j,k,i) = (1/Nmc)*sum(CovMatrixTempNF(:,j,k,i));   
%     
%         end
%     end    
%     
% end
% 
% %end
% 
% 
% 
% 
% 
% 
% %Calculate variances of the errors
% for ko = startepoch:endepoch
% for j = 1:9
% 
%         
% 
% var_x_state_total_minusINSerrorVectorNFmc(j,ko) = (1/Nmc)*sum(x_state_total_minusINSerrorVectorNFmc(:,j,ko).^2);
% 
% 
% var_x_state_total_minusMODELerrorVectorNFmc(j,ko) = (1/Nmc)*sum(x_state_total_minusMODELerrorVectorNFmc(:,j,ko).^2);
%     
%     
% 
% end
% 
% end
% 
% 
% 
% 
% 
% 
% 
% 
% 
% for i = startepoch:endepoch
%     
%     
%     for j = 1:9
%         for k = 1:9
%             
%             
%    % CovMatrixCorrCoef(j,k,i) = CovMatrix(j,k,i)/(sqrt(CovMatrix(j,j,i))*sqrt(CovMatrix(k,k,i)));
%    
%                 
%     CovMatrixCorrCoefNF(j,k,i) = CovMatrixNF(j,k,i)/(sqrt(var_x_state_total_minusINSerrorVectorNFmc(j,i))*sqrt(var_x_state_total_minusMODELerrorVectorNFmc(k,i)));
%    
%     
%         end
%     end
%     
%     
% end
% 
% 
% 
% for i = startepoch:endepoch
%     
%     CovMatrixCorrCoefNFHgt(i) = CovMatrixCorrCoefNF(1,1,i);
%     
% end
% 
% 
% 
% 
% 
% for i = startepoch:endepoch
%     
%     CovMatrixCorrCoefHgt(i) = CovMatrixCorrCoef(1,1,i);
%     
% end
% 
% 
% 
% 
% end %if processdata == 1
% 

%quit; %quit matlab when all done



