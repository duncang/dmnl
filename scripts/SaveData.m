%This saves the data when called periodically from MonteCarloWrapper
%function, for use with HPC for example. 

HPL_INS_NFmc(iM,:) = HPL_INS_NF;
HPL_INS_Cmc(iM,:) = HPL_INS_C;
VPL_INS_NFmc(iM,:) = VPL_INS_NF;
VPL_INS_Cmc(iM,:) = VPL_INS_C;



HPL_MODEL_NFmc(iM,:) = HPL_MODEL_NF;
HPL_MODEL_Cmc(iM,:) = HPL_MODEL_C;
VPL_MODEL_NFmc(iM,:) = VPL_MODEL_NF;
VPL_MODEL_Cmc(iM,:) = VPL_MODEL_C;




HPE_INS_NFmc(iM,:) = HPE_INS_NF;
HPE_INS_Cmc(iM,:) = HPE_INS_C;
VPE_INS_NFmc(iM,:) = VPE_INS_NF;
VPE_INS_Cmc(iM,:) = VPE_INS_C;



HPE_MODEL_NFmc(iM,:) = HPE_MODEL_NF;
HPE_MODEL_Cmc(iM,:) = HPE_MODEL_C;
VPE_MODEL_NFmc(iM,:) = VPE_MODEL_NF;
VPE_MODEL_Cmc(iM,:) = VPE_MODEL_C;



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



%for INS fused
 sigma_vNFmc_INS(iM,:) = sigma_v_INSNF;
        sigma_v_maxNFmc_INS(iM,:) = sigma_v_max_INSNF;        
        dmajorNFmc_INS(iM,:) = dmajor_INSNF;
        dmajor_maxNFmc_INS(iM,:) = dmajor_max_INSNF;
    

        
 %for ADM unfused     
        sigma_v_NFmc_M(iM,:)  = sigma_v_MNF;
        sigma_v_max_NFmc_M(iM,:) = sigma_v_max_MNF;        
        dmajor_NFmc_M(iM,:) = dmajor_MNF;
        dmajor_max_NFmc_M(iM,:) = dmajor_max_MNF;
        
        



 HPE_GPS_mc(iM,:) = HPE_GPS_C;    %there's no difference between fused and not fused for the GPS
 VPE_GPS_mc(iM,:) = VPE_GPS_C;




x_state_total_INSBIAS_NFmc(iM,:,:) =  x_state_total_INSBIAS_NF;
x_state_total_INSBIAS_Cmc(iM,:,:) =  x_state_total_INSBIAS;


x_state_total_MODELBIAS_NFmc(iM,:,:) =  x_state_total_MODELBIAS_NF;
x_state_total_MODELBIAS_Cmc(iM,:,:) =  x_state_total_MODELBIAS;



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



corrcoeffAtt_NFmc(iM,:,:) = corrcoeffAtt_NF;
corrcoeffVel_NFmc(iM,:,:) = corrcoeffVel_NF;
corrcoeffPos_NFmc(iM,:,:) = corrcoeffPos_NF;

corrcoeffCLKbias_NFmc(iM,:,:) = corrcoeffCLKbias_NF;
corrcoeffCLKdrift_NFmc(iM,:,:) = corrcoeffCLKdrift_NF;

corrcoeffAccbias_NFmc(iM,:,:) = corrcoeffAccbias_NF;
corrcoeffGyrobias_NFmc(iM,:,:) = corrcoeffGyrobias_NF;


corrcoeffAtt_Cmc(iM,:,:) = corrcoeffAtt;
corrcoeffVel_Cmc(iM,:,:) = corrcoeffVel;
corrcoeffPos_Cmc(iM,:,:) = corrcoeffPos;

corrcoeffCLKbias_Cmc(iM,:,:) = corrcoeffCLKbias;
corrcoeffCLKdrift_Cmc(iM,:,:) = corrcoeffCLKdrift;

corrcoeffAccbias_Cmc(iM,:,:) = corrcoeffAccbias;
corrcoeffGyrobias_Cmc(iM,:,:) = corrcoeffGyrobias;




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
N_save_NFmc(iM,:,:,:) = N_save_NF;

 
 
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
