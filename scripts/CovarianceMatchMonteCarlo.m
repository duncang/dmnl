%applies covariance matching technique for monte carlo





%===============================================
%do covariance matching to see if HPL meets Pmd

%===============================================

%Pr that pr md = Pmd = Pr that HPE > HPL and pr that test stat < detection
%threshold. 




startepochplot = startepoch + 2;
endepochplot = endepoch

starti = startepochplot;
endi = endepochplot;




for i = starti:endi        
   
    %note, Pr test stat < threshold should be high, ie near 1, because its
    %pr test stat is less than, not more than, threshold. PrHPE should be
    %low, because its prob that hpe exceeds hal. 
    
    %for IMU
    
    
  %  lambda_k = max(lambda_ss_outPLOTINS_H_C(i,SubToPlot)); %noncentrality parameter is max of test statistics.
    
    
        lambda_k = (lambda_ss_outPLOTINS_H_C(i,SubToPlot)^2); %noncentrality parameter is max of test statistics. %removed max in case its false alarms, i specifiy the value
    %threshold squared
   TD_squared = TD_outPLOT_H_C(i)^2;
    N = N_save_Cmc(29,i);
    DOF = 2;  %for HPL
  sigma_re = sqrt((sqrt(P_save_NFmc(29,7,7,i)))^2 + (sqrt(P_save_NFmc(29,8,8,i)))^2)*ReConst;
    PL = HPL_INS_NFmc(29,i);
    de = HPE_INS_NFmc(29,i);        
    
    [Pmd_Est_INS_H(i), PrHPE_INS_H(i), PrTestStat_INS_H(i)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
       
    
      
 
    
    
    
    
    
    
    
end


%answers the question, is Pmd satisfied 

    
    
    
    
      lambda_k = max(lambda_ss_outINS_V(i,:)); %noncentrality parameter is max of test statistics.
    %threshold squared
   TD_squared = aTDINS_V(i);
    N = N_save(i);
    DOF = 1;  %for VPL
  sigma_re = sqrt(P_save(9,9,i));
    PL = VPL_INS_C(i);
    de = VPE_INS_C(i);        
    
    [Pmd_Est_INS_V(i), PrVPE_INS_V(i), PrTestStat_INS_V(i)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
    
    
    %for ADM
    
     lambda_k = max(lambda_ss_outMODEL_H(i,:)); %noncentrality parameter is max of test statistics.
    %threshold squared
   TD_squared = aTDINS_H(i);
    N = N_save(i);
    DOF = 2;  %for HPL
  sigma_re = sqrt((sqrt(P_save(24,24,i)))^2 + (sqrt(P_save(25,25,i)))^2)*ReConst;
    PL = HPL_MODEL_C(i);
    de = HPE_MODEL_C(i);        
    
    [Pmd_Est_MODEL_H(i), PrHPE_MODEL_H(i), PrTestStatMODELS_H(i)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
       
    
    
    
      lambda_k = max(lambda_ss_outMODEL_V(i,:)); %noncentrality parameter is max of test statistics.
    %threshold squared
   TD_squared = aTDINS_V(i);
    N = N_save(i);
    DOF = 1;  %for VPL
  sigma_re = sqrt(P_save(26,26,i));
    PL = VPL_MODEL_C(i);
    de = VPE_MODEL_C(i);        
    
    [Pmd_Est_MODEL_V(i), PrVPE_MODEL_V(i), PrTestStat_MODEL_V(i)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
    
    
    
    %for GPS
    
     
    
     lambda_k = max(lambda_ss_outGPS_H(i,:)); %noncentrality parameter is max of test statistics.
    %threshold squared
   TD_squared = aTDGPS_H(i);
    N = N_save(i);
    DOF = 2;  %for HPL

    
sigma_re = sqrt((VarSolutionVec_Subsave(i,1)*SigmaPR^2)^2 + (VarSolutionVec_Subsave(i,2)*SigmaPR^2)^2) 
    
    PL = HPL_GPS(i);
    de = HPE_GPS_C(i);        
    
    [Pmd_Est_GPS_H(i), PrHPE_GPS_H(i), PrTestStat_GPS_H(i)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
       
    
      
    
    
      lambda_k = max(lambda_ss_outGPS_V(i,:)); %noncentrality parameter is max of test statistics.
    %threshold squared
   TD_squared = aTDGPS_V(i);
    N = N_save(i);
    DOF = 1;  %for VPL

    
sigma_re = sqrt(VarSolutionVec_Subsave(i,3)*SigmaPR^2)
      
  
    PL = VPL_GPS(i);
    de = VPE_GPS_C(i);        
    
    [Pmd_Est_GPS_V(i), PrVPE_GPS_V(i), PrTestStat_GPS_V(i)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
    
    
    
    
end



SubToPlot = MaxSlopeSatelliteToDo_mc(29);
%plots

figure();
hold

plot(HPL_INS_Cmc(29,startepochplot:endepochplot),'k', 'LineWidth',2);


plot(HPE_INS_Cmc(29,startepochplot:endepochplot),'k--', 'LineWidth',2);


plot(((lambda_ss_outPLOTINS_H_C(startepochplot:endepochplot,SubToPlot))),'k-.','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_H_C(startepochplot:endepochplot)),'k+','LineWidth',2);    %not sure if this is right or not



title 'Horizontal Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('HPL GPS-IMU-ADM EKF (m)', 'HPE GPS-IMU-ADM EKF (m)','Test statistic','Threshold');
axis([30 80 0 20]);









figure();
hold
plot(((lambda_ss_outPLOTINS_V_C(startepochplot:endepochplot,SubToPlot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_V_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not

plot(VPL_INS_Cmc(29,startepochplot:endepochplot),'b', 'LineWidth',2);
plot(VPE_INS_Cmc(29,startepochplot:endepochplot),'g', 'LineWidth',2);

title 'Horizontal Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('IMU test statistic','Threshold', 'HPL IMU (m)', 'HPE IMU (m)');
axis([0 80 0 40]);











SubToPlot = MaxSlopeSatelliteToDo_mc(69);
%plots

figure();
hold

plot(HPL_INS_Cmc(69,startepochplot:endepochplot),'k', 'LineWidth',2);


plot(HPE_INS_Cmc(69,startepochplot:endepochplot),'k--', 'LineWidth',2);


plot(((lambda_ss_outPLOTINS_H_C(startepochplot:endepochplot,SubToPlot))),'k-.','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_H_C(startepochplot:endepochplot)),'k+','LineWidth',2);    %not sure if this is right or not



title 'Horizontal Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('HPL GPS-IMU-ADM EKF (m)', 'HPE GPS-IMU-ADM EKF (m)','Test statistic','Threshold');
axis([30 80 0 20]);









%say, to test whether or not Pmd is met we calculate the probability of
%missed detection using the covariance matching method , for up until the
%time when the fault is detected. 