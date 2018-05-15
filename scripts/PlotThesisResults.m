%this is plots to use inthe thesis, took it out of the newjuler3 code, so
%its more generic and easier to follow


startepochplot = startepoch+2;  %start it at i = 5
endepochplot = endepoch;
%visibility plot
plot(N_save(startepochplot:endepochplot),'k','LineWidth',2);  title 'Satellite Visibility' ; xlabel('Time (s)');  ylabel('Number of Satellites');
axis([0 endepochplot 0 10]);






%PLOTS for THESIS

% % 
% % plot from 20 : 137, takes about 50 seconds total to converge to steady
% % state....
% % 





HPLImprovement = HPL_INS_NF(startepoch:endepoch) -HPL_INS_C(startepoch:endepoch);
VPLImprovement = VPL_INS_NF(startepoch:endepoch) -VPL_INS_C(startepoch:endepoch);

HPLImprovepercent = (HPLImprovement./HPL_INS_NF(startepoch:endepoch))*100;
VPLImprovepercent = (VPLImprovement./VPL_INS_NF(startepoch:endepoch))*100;


mean(HPLImprovepercent(startepochplot:endepochplot-2))
mean(VPLImprovepercent(startepochplot:endepochplot-2))








%to work out the improvement %'s in the trace of the covariances..


HPLImprovement1 = tracePminusNF_INS(startepoch:endepoch) -tracePminusC_INS(startepoch:endepoch);
HPLImprovement2 = tracePminusNF_ADM(startepoch:endepoch) -tracePminusC_ADM(startepoch:endepoch);
HPLImprovement3 = tracePNF_INS(startepoch:endepoch) -tracePC_INS(startepoch:endepoch);
HPLImprovement4 = tracePNF_ADM(startepoch:endepoch) -tracePC_ADM(startepoch:endepoch);




HPLImprovepercent1 = (HPLImprovement1./tracePminusNF_INS(startepoch:endepoch))*100;
HPLImprovepercent2 = (HPLImprovement2./tracePminusNF_ADM(startepoch:endepoch))*100;
HPLImprovepercent3 = (HPLImprovement3./tracePNF_INS(startepoch:endepoch))*100;
HPLImprovepercent4 = (HPLImprovement4./tracePNF_ADM(startepoch:endepoch))*100;



mean(HPLImprovepercent1(startepochplot:endepochplot-2))
mean(HPLImprovepercent2(startepochplot:endepochplot-2))
mean(HPLImprovepercent3(startepochplot:endepochplot-2))
mean(HPLImprovepercent4(startepochplot:endepochplot-2))















HPLImprovementADM = HPL_MODEL_NF(startepoch:endepoch) -HPL_MODEL_C(startepoch:endepoch);
VPLImprovementADM = VPL_MODEL_NF(startepoch:endepoch) -VPL_MODEL_C(startepoch:endepoch);

HPLImprovepercentADM = (HPLImprovementADM./HPL_MODEL_NF(startepoch:endepoch))*100;
VPLImprovepercentADM = (VPLImprovementADM./VPL_MODEL_NF(startepoch:endepoch))*100;


mean(HPLImprovepercentADM(startepochplot:endepochplot-2))
mean(VPLImprovepercentADM(startepochplot:endepochplot-2))






%mean percent improvement of IMU over GPS. 


HPLImprovement = HPL_GPS(startepoch:endepoch) -HPL_INS_NF(startepoch:endepoch);
VPLImprovement = VPL_GPS(startepoch:endepoch) -VPL_INS_NF(startepoch:endepoch);

HPLImprovepercent = (HPLImprovement./HPL_GPS(startepoch:endepoch))*100;
VPLImprovepercent = (VPLImprovement./VPL_GPS(startepoch:endepoch))*100;


mean(HPLImprovepercent(startepochplot:endepochplot-2))
mean(VPLImprovepercent(startepochplot:endepochplot-2))


%mean percent improvement of ADM over GPS. 


HPLImprovement = HPL_GPS(startepoch:endepoch) -HPL_MODEL_NF(startepoch:endepoch);
VPLImprovement = VPL_GPS(startepoch:endepoch) -VPL_MODEL_NF(startepoch:endepoch);

HPLImprovepercent = (HPLImprovement./HPL_GPS(startepoch:endepoch))*100;
VPLImprovepercent = (VPLImprovement./VPL_GPS(startepoch:endepoch))*100;


mean(HPLImprovepercent(startepochplot:endepochplot-2))
mean(VPLImprovepercent(startepochplot:endepochplot-2))




%mean percent improvement of IMU over ADM. 


HPLImprovement = HPL_MODEL_NF(startepoch:endepoch) -HPL_INS_NF(startepoch:endepoch);
VPLImprovement = VPL_MODEL_NF(startepoch:endepoch) -VPL_INS_NF(startepoch:endepoch);

HPLImprovepercent = (HPLImprovement./HPL_MODEL_NF(startepoch:endepoch))*100;
VPLImprovepercent = (VPLImprovement./VPL_MODEL_NF(startepoch:endepoch))*100;


mean(HPLImprovepercent(startepochplot:endepochplot-2))
mean(VPLImprovepercent(startepochplot:endepochplot-2))












%plot det of P minus


%convert Pminus save into position 


Pminus_save_NFcopy = Pminus_save_NF;
P_minus_savecopy = P_minus_save;

for i = startepoch:endepoch
        

            [Rn, Re] = WGS84_calcRnRe(Lat_truth1Hz(i));
      
            Pminus_save_NFcopy(7:8,7:8,i) = Pminus_save_NFcopy(7:8,7:8,i)*Rn*Re;
P_minus_savecopy(7:8,7:8,i) = P_minus_savecopy(7:8,7:8,i)*Rn*Re;

      
            Pminus_save_NFcopy(24:25,24:25,i) = Pminus_save_NFcopy(24:25,24:25,i)*Rn*Re;
P_minus_savecopy(24:25,24:25,i) = P_minus_savecopy(24:25,24:25,i)*Rn*Re;        
            
end


for moo = startepochplot:endepochplot
    
       
    
    detPminusNF_INS(moo) = det(Pminus_save_NFcopy(7:9,7:9,moo));
    
    detPminusC_INS(moo) = det(P_minus_savecopy(7:9,7:9,moo));
    
    detPminusNF_ADM(moo) = det(Pminus_save_NFcopy(24:26,24:26,moo));
    
    detPminusC_ADM(moo) = det(P_minus_savecopy(24:26,24:26,moo));

end

for moo = startepochplot:endepochplot
    
       
    
    tracePminusNF_INS(moo) = trace(Pminus_save_NFcopy(7:9,7:9,moo));
    
    tracePminusC_INS(moo) = trace(P_minus_savecopy(7:9,7:9,moo));
    
    tracePminusNF_ADM(moo) = trace(Pminus_save_NFcopy(24:26,24:26,moo));
    
    tracePminusC_ADM(moo) = trace(P_minus_savecopy(24:26,24:26,moo));

end



%plot det of P plus


%convert Pminus save into position 

P_save_NFcopy = P_save_NF;
P_savecopy = P_save;


for i = startepoch:endepoch        

            [Rn, Re] = WGS84_calcRnRe(Lat_truth1Hz(i));
      
            P_save_NFcopy(7:8,7:8,i) = P_save_NFcopy(7:8,7:8,i)*Rn*Re;
P_savecopy(7:8,7:8,i) = P_savecopy(7:8,7:8,i)*Rn*Re;

      
            P_save_NFcopy(24:25,24:25,i) = P_save_NFcopy(24:25,24:25,i)*Rn*Re;
P_savecopy(24:25,24:25,i) = P_savecopy(24:25,24:25,i)*Rn*Re;

                       
end



for moo = startepochplot:endepochplot    
       
    
    detPNF_INS(moo) = det(P_save_NFcopy(7:9,7:9,moo));
    
    detPC_INS(moo) = det(P_savecopy(7:9,7:9,moo));    
    
      detPNF_ADM(moo) = det(P_save_NFcopy(24:26,24:26,moo));
    
    detPC_ADM(moo) = det(P_savecopy(24:26,24:26,moo));

end



for moo = startepochplot:endepochplot    
       
    
    tracePNF_INS(moo) = trace(P_save_NFcopy(7:9,7:9,moo));
    
    tracePC_INS(moo) = trace(P_savecopy(7:9,7:9,moo));    
    
      tracePNF_ADM(moo) = trace(P_save_NFcopy(24:26,24:26,moo));
    
    tracePC_ADM(moo) = trace(P_savecopy(24:26,24:26,moo));
    
    
    GPS_var(1:3,1:3,moo) = AA_Full(moo,1:3,1:3)*SigmaPR^2;
    
    
    
    trace_GPS(moo) = trace(GPS_var(1:3,1:3,moo));
    
    
    

end





%calculate trace of GPS 



for moo = startepochplot:endepochplot    

    
    
    stdevPnINS_minus(moo) = sqrt(Pminus_save_NFcopy(7,7,moo));


    stdevPeINS_minus(moo) = sqrt(Pminus_save_NFcopy(8,8,moo));

    stdevPdINS_minus(moo) = sqrt(Pminus_save_NFcopy(9,9,moo));

    
    H_ADM_StdDevMinusINS(moo) = sqrt(stdevPnINS_minus(moo)^2 + stdevPeINS_minus(moo)^2) ;
    V_ADM_StdDevMinusINS(moo) = stdevPdINS_minus(moo);
    
      

    stdevPnADM_minus(moo) = sqrt(Pminus_save_NFcopy(24,24,moo));


    stdevPeADM_minus(moo) = sqrt(Pminus_save_NFcopy(25,25,moo));

    stdevPdADM_minus(moo) = sqrt(Pminus_save_NFcopy(26,26,moo));

    
    H_ADM_StdDevMinusADM(moo) = sqrt(stdevPnADM_minus(moo)^2 + stdevPeADM_minus(moo)^2) ;
    V_ADM_StdDevMinusADM(moo) = stdevPdADM_minus(moo);
    
     
        
    
    stdevPnINS_plus(moo) = sqrt(P_save_NFcopy(7,7,moo));

    stdevPeINS_plus(moo) = sqrt(P_save_NFcopy(8,8,moo));

    stdevPdINS_plus(moo) = sqrt(P_save_NFcopy(9,9,moo));

    
    H_ADM_StdDevplusINS(moo) = sqrt(stdevPnINS_plus(moo)^2 + stdevPeINS_plus(moo)^2) ;
    V_ADM_StdDevplusINS(moo) = stdevPdINS_plus(moo);
    
    
      

    stdevPnADM_plus(moo) = sqrt(P_save_NFcopy(24,24,moo));


    stdevPeADM_plus(moo) = sqrt(P_save_NFcopy(25,25,moo));

    stdevPdADM_plus(moo) = sqrt(P_save_NFcopy(26,26,moo));

    
    H_ADM_StdDevplusADM(moo) = sqrt(stdevPnADM_plus(moo)^2 + stdevPeADM_plus(moo)^2) ;
    V_ADM_StdDevplusADM(moo) = stdevPdADM_plus(moo);
    
           
    
    
end







 




%calculate quickly the std of the p q r Ax Ay Az error




std((omega_xMODEL1Hz(startepochplot:endepochplot) - GyroTruth1Hz(1,startepochplot:endepochplot))*180/pi)


std((omega_yMODEL1Hz(startepochplot:endepochplot) - GyroTruth1Hz(2,startepochplot:endepochplot))*180/pi)


std((omega_zMODEL1Hz(startepochplot:endepochplot) - GyroTruth1Hz(3,startepochplot:endepochplot))*180/pi)




std((ax_b_MODEL1Hz(startepochplot:endepochplot) - AccelTruth1Hz(1,startepochplot:endepochplot))) 


std((ay_b_MODEL1Hz(startepochplot:endepochplot) - AccelTruth1Hz(2,startepochplot:endepochplot))) 


std((az_b_MODEL1Hz(startepochplot:endepochplot) - AccelTruth1Hz(3,startepochplot:endepochplot))) 









% 1.2.1	Comparison between GPS, GPS-IMU, and GPS-ADM
% 
% 
% 
% plot of covariance of GPS, GPS-IMU, and GPS-ADM for horizontal for apriori covariances only

%plot the actual standard deviations or trace here because the determinant is too large for the ADM for using det.. 




figure();
hold;
plot(tracePminusNF_INS(startepochplot:endepochplot), 'b','LineWidth',2); title 'Trace of P^- IMU, ADM, and GPS' ; xlabel('Time (s)');  ylabel('trace(P^-) (m)');
plot(tracePminusNF_ADM(startepochplot:endepochplot), 'r','LineWidth',2); 
plot(trace_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
legend('P^- IMU','P^- ADM','P GPS');       
axis([0 endepochplot-startepochplot 0 40]);


% 
% plot of covariance of GPS, GPS-IMU, and GPS-ADM for vertical for apriori covariances only
% 
% 
% 
% 
% plot of covariance of GPS, GPS-IMU, and GPS-ADM for horizontal
% 
% plot of covariance of GPS, GPS-IMU, and GPS-ADM for vertical
% 
% 

figure();
hold;
plot(tracePNF_INS(startepochplot:endepochplot), 'b','LineWidth',2); title 'Trace of P IMU, ADM, and GPS' ; xlabel('Time (s)');  ylabel('trace(P) (m)');
plot(tracePNF_ADM(startepochplot:endepochplot), 'r','LineWidth',2); 
plot(trace_GPS(startepochplot:endepochplot),'g','LineWidth',2); 

legend('P IMU','P ADM','P GPS');       

axis([0 endepochplot-startepochplot 0 35]);









%plot accuracy things





figure;
hold;
title 'Horizontal Accuracy GPS, IMU, ADM' ; xlabel('Time (s)'); ylabel('Horizontal Position Error 1 Sigma (m)');

plot(dmajor_INSNF(startepochplot:endepochplot),'b','LineWidth',2);
plot(dmajor_MNF(startepochplot:endepochplot),'r','LineWidth',2);
plot(dmajor_GPS(startepochplot:endepochplot),'g','LineWidth',2);
plot(dmajor_max_INSNF(startepochplot:endepochplot),'k','LineWidth',4);

legend('IMU','ADM','GPS');      


figure;
hold;
title 'Vertical Accuracy GPS, IMU, ADM' ; xlabel('Time (s)'); ylabel('Vertical Position Error 1 Sigma (m)');
plot(sigma_v_INSNF(startepochplot:endepochplot),'b','LineWidth',2);
plot(sigma_v_MNF(startepochplot:endepochplot),'r','LineWidth',2);
plot(sigma_v_GPS(startepochplot:endepochplot),'g','LineWidth',2);
plot(sigma_v_max_INSNF(startepochplot:endepochplot),'k','LineWidth',4);
legend('IMU','ADM','GPS');     







% 
% 
% plot of protection levels of GPS, GPS-IMU, and GPS-ADM for horizontal. Also plot the lower bound (ie NSE, see Lee how to do this)
% 
% plot of protection levels of GPS, GPS-IMU, and GPS-ADM for vertical




% plot of protection level of GPS-IMU, and GPS-IMU-ADM for horizontal
% 
% plot of protection level of GPS-IMU, and GPS-IMU-ADM for vertical
% 


figure();
hold;
plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(HPL_INS_NF(startepochplot:endepochplot),'b','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot(HPL_MODEL_NF(startepochplot:endepochplot),'r','LineWidth',2);  
plot(HPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('HAL','HPL IMU','HPL ADM', 'HPL GPS');    

axis([0 endepochplot-startepochplot 0 45]);




figure();
hold;
plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(VPL_INS_NF(startepochplot:endepochplot),'b','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot(VPL_MODEL_NF(startepochplot:endepochplot),'r','LineWidth',2);  
plot(VPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('VAL','VPL IMU','VPL ADM', 'VPL GPS');       

axis([0 endepochplot-startepochplot 0 60]);








%plots for comparing with wind estimation





figure();
hold;
plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
%plot(HPL_INS_NF(startepochplot:endepochplot),'b','LineWidth',2);  

title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot(HPLnowindest(startepochplot:endepochplot),'g','LineWidth',2);  
plot(HPLnoad(startepochplot:endepochplot),'b','LineWidth',2);  
plot(HPLwithad(startepochplot:endepochplot),'r','LineWidth',2);  
%plot(HPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('HAL', 'HPL ADM', 'HPL ADM Wind', 'HPL ADM AD', 'HPL GPS');    

axis([0 endepochplot-startepochplot 0 120]);




figure();
hold;
plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
%plot(VPL_INS_NF(startepochplot:endepochplot),'b','LineWidth',2);  
title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot(VPLnowindest(startepochplot:endepochplot),'g','LineWidth',2);  
plot(VPLnoad(startepochplot:endepochplot),'b','LineWidth',2);  
plot(VPLwithad(startepochplot:endepochplot),'r','LineWidth',2);  

%plot(VPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('VAL', 'VPL ADM', 'VPL ADM Wind', 'VPL ADM AD', 'VPL GPS');     

axis([0 endepochplot-startepochplot 0 120]);









% 
% 
% 1.2.2	Effect of MMF with ADM 
% 
% plot of covariance of GPS-IMU, and GPS-IMU-ADM for horizontal
% 
% plot of covariance of GPS-IMU, and GPS-IMU-ADM for vertical
% 


%apriori covariances

figure();
hold;
plot(tracePminusNF_INS(startepochplot:endepochplot), 'k','LineWidth',2); title 'Trace of P^- IMU' ; xlabel('Time (s)');  ylabel('Trace(P^-) (m)');
plot(tracePminusC_INS(startepochplot:endepochplot), 'g','LineWidth',2);


legend('P^- IMU','P^- IMU with MMF');       
axis([0 endepochplot-startepochplot 10 15]);


figure();
hold;
plot(tracePminusNF_ADM(startepochplot:endepochplot), 'k','LineWidth',2); title 'Trace of P^- ADM' ; xlabel('Time (s)');  ylabel('Trace(P^-) (m)');
plot(tracePminusC_ADM(startepochplot:endepochplot), 'g','LineWidth',2);


legend('P^- ADM','P^- ADM with MMF');    
axis([0 endepochplot-startepochplot 15 30]);




%posteriori covariances


figure();
hold;
plot(tracePNF_INS(startepochplot:endepochplot), 'k','LineWidth',2); title 'Trace of P IMU' ; xlabel('Time (s)');  ylabel('Trace(P) (m)');
plot(tracePC_INS(startepochplot:endepochplot), 'g','LineWidth',2);
legend('P IMU','P IMU with MMF');       

axis([0 endepochplot-startepochplot 3.5 7]);




figure();
hold;
plot(tracePNF_ADM(startepochplot:endepochplot), 'k','LineWidth',2); title 'Trace of P ADM' ; xlabel('Time (s)');  ylabel('Trace(P) (m)');
plot(tracePC_ADM(startepochplot:endepochplot), 'g','LineWidth',2);
legend('P ADM','P ADM with MMF');   

axis([0 endepochplot-startepochplot 0 10]);




figure;
hold;
title 'Horizontal Accuracy with MMF' ; xlabel('Time (s)');
plot(dmajor_INS(startepochplot:endepochplot),'b','LineWidth',2);
plot(dmajor_INSNF(startepochplot:endepochplot),'b','LineWidth',2);
plot(dmajor_MNF(startepochplot:endepochplot),'r','LineWidth',2);
plot(dmajor_GPS(startepochplot:endepochplot),'g','LineWidth',2);
plot(dmajor_max_INSNF(startepochplot:endepochplot),'k','LineWidth',4);
legend('IMU', 'IMU MMF', 'ADM','GPS');   
axis([0 endepochplot-startepochplot 0 3]);


figure;
hold;
title 'Vertical Accuracy with MMF' ; xlabel('Time (s)');
plot(sigma_v_INS(startepochplot:endepochplot),'b','LineWidth',2);
plot(sigma_v_INSNF(startepochplot:endepochplot),'b--','LineWidth',2);
plot(sigma_v_MNF(startepochplot:endepochplot),'r--','LineWidth',2);
plot(sigma_v_GPS(startepochplot:endepochplot),'g--','LineWidth',2);
plot(sigma_v_max_INSNF(startepochplot:endepochplot),'k','LineWidth',4);
legend('IMU', 'IMU MMF', 'ADM','GPS');   
axis([0 endepochplot-startepochplot 0 3]);








% plot of protection level of GPS-IMU, and GPS-IMU-ADM for horizontal
% 
% plot of protection level of GPS-IMU, and GPS-IMU-ADM for vertical
% 


figure();
hold;
%plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(HPL_INS_NF(startepochplot:endepochplot),'k--','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot(HPL_INS_C(startepochplot:endepochplot),'k','LineWidth',2);  
%plot(HPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

%legend('HAL','HPL IMU','HPL IMU-ADM', 'HPL GPS');       


legend('HPL IMU','HPL IMU-ADM');    

axis([0 endepochplot-startepochplot 6 18]);


figure();
hold;
%plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(VPL_INS_NF(startepochplot:endepochplot),'k--','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot(VPL_INS_C(startepochplot:endepochplot),'k','LineWidth',2);  
%plot(VPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

%legend('VAL','VPL IMU','VPL IMU-ADM', 'VPL GPS');       


legend('VPL IMU','VPL IMU-ADM');    

axis([0 endepochplot-startepochplot 9 22]);



























%==========================
%FILTER TUNING FOR THESIS
%===========================



mean(HPE_INS_NF)
mean(HPE_INS_C)


mean(VPE_INS_NF)
mean(VPE_INS_C)


var(HPE_INS_NF)
var(HPE_INS_C)

var(VPE_INS_NF)
var(VPE_INS_C)

%consistency tests moved to ConsistencyTests.m
   
    










%PLOTS for test statistic etc for THESIS



%plot test stat of the most difficult to detect satellite.. 
%thats the one the fault is on. 


%determine max test stat




%find maxteststat 


  sizeaa = size(lambda_ss_outPLOTINS_H_C);
    sizeaatt = sizeaa(2);
    

%this determines which subfilter has the biggest test statistic, for
%plotting
for i = startepochplot:endepochplot        
    for n = 1:sizeaatt   
        
        maxsubfilter_value(i,n) = lambda_ss_outPLOTINS_H_C(i,n);
    end

maxsubfilter(i) = max(maxsubfilter_value(i,:));

for n = 1:sizeaatt
    
    if maxsubfilter_value(i,n) == maxsubfilter(i) 
        
        maxsubfilter_n(i) = n;
        
    end

end

end





%this determines which subfilter has the biggest test statistic, for
%plotting
for i = startepochplot:endepochplot        
    for n = 1:sizeaatt   
        
        maxsubfilter_value(i,n) = lambda_ss_outPLOTINS_V_C(i,n);
    end

maxsubfilter(i) = max(maxsubfilter_value(i,:));

for n = 1:sizeaatt
    
    if maxsubfilter_value(i,n) == maxsubfilter(i) 
        
        maxsubfilter_n(i) = n;
        
    end

end

end





%this determines which subfilter has the biggest test statistic, for
%plotting
for i = startepochplot:endepochplot        
    for n = 1:sizeaatt   
        
        maxsubfilter_value(i,n) = lambda_ss_outPLOTMODEL_H_C(i,n);
    end

maxsubfilter(i) = max(maxsubfilter_value(i,:));

for n = 1:sizeaatt
    
    if maxsubfilter_value(i,n) == maxsubfilter(i) 
        
        maxsubfilter_n(i) = n;
        
    end

end

end





%this determines which subfilter has the biggest test statistic, for
%plotting
for i = startepochplot:endepochplot        
    for n = 1:sizeaatt   
        
        maxsubfilter_value(i,n) = lambda_ss_outPLOTMODEL_V_C(i,n);
    end

maxsubfilter(i) = max(maxsubfilter_value(i,:));

for n = 1:sizeaatt
    
    if maxsubfilter_value(i,n) == maxsubfilter(i) 
        
        maxsubfilter_n(i) = n;
        
    end

end

end



    

%this determines which subfilter has the biggest test statistic, for
%plotting, it also determines where it crosses the line
for i = startepochplot:endepochplot        
    for n = 1:sizeaatt   
        
        maxsubfilter_value(i,n) = lambda_ss_outPLOTGPS_H_C(i,n);
    end

maxsubfilter(i) = max(maxsubfilter_value(i,:));

for n = 1:sizeaatt
    
    if maxsubfilter_value(i,n) == maxsubfilter(i) 
        
        maxsubfilter_n(i) = n;
        
    end

end

end




sizeaa = size(lambda_ss_outPLOTINS_H_C);
sizeaatt = sizeaa(2);
    

%this automatically determines which subfilter has the biggest test statistic, for
%plotting, it also determines where it crosses the line
for i = startepochplot:endepochplot        
    for n = 1:sizeaatt   
        
        maxsubfilter_value(i,n) = lambda_ss_outPLOTGPS_V_C(i,n);
    end

maxsubfilter(i) = max(maxsubfilter_value(i,:));

for n = 1:sizeaatt
    
    if maxsubfilter_value(i,n) == maxsubfilter(i) 
        
        maxsubfilter_n(i) = n;
        
    end

end

end











SubToPlot = MaxSlopeSatelliteToDo;

%SubToPlot = 12;


figure();
hold
plot(((lambda_ss_outPLOTINS_H_C(startepochplot:endepochplot,SubToPlot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_H_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not

plot(HPL_INS_C(startepochplot:endepochplot),'b', 'LineWidth',2);
plot(HPE_INS_C(startepochplot:endepochplot),'g', 'LineWidth',2);

title 'Horizontal Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('IMU test statistic','Threshold', 'HPL IMU (m)', 'HPE IMU (m)');
axis([30 80 0 15]);

%PLOT VERTICAL TEST STATS


figure();
hold
plot(((lambda_ss_outPLOTINS_V_C(startepochplot:endepochplot,SubToPlot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_V_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not

plot(VPL_INS_C(startepochplot:endepochplot),'b', 'LineWidth',2);
plot(VPE_INS_C(startepochplot:endepochplot),'g', 'LineWidth',2);
title 'Vertical Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('IMU test statistic','Threshold', 'VPL IMU (m)', 'VPE IMU (m)');
axis([30 65 0 23]);






figure();
hold
plot(((lambda_ss_outPLOTMODEL_H_C(startepochplot:endepochplot,SubToPlot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_H_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not

plot(HPL_MODEL_C(startepochplot:endepochplot),'b', 'LineWidth',2);
plot(HPE_MODEL_C(startepochplot:endepochplot),'g', 'LineWidth',2);

title 'Horizontal Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('ADM test statistic','Threshold','HPL ADM (m)', 'HPE ADM (m)');
axis([30 90 0 18]);

%PLOT VERTICAL TEST STATS


figure();
hold
plot(((lambda_ss_outPLOTMODEL_V_C(startepochplot:endepochplot,SubToPlot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_V_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not

plot(VPL_MODEL_C(startepochplot:endepochplot),'b', 'LineWidth',2);
plot(VPE_MODEL_C(startepochplot:endepochplot),'g', 'LineWidth',2);

title 'Vertical Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('ADM test statistic','Threshold','VPL ADM (m)', 'VPE ADM (m)');
axis([30 80 0 20]);




%FOR GPS




figure();
hold
plot(((lambda_ss_outPLOTGPS_H_C(startepochplot:endepochplot,SubToPlot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_H_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not

plot(HPL_GPS(startepochplot:endepochplot),'b', 'LineWidth',2);
plot(HPE_GPS_C(startepochplot:endepochplot),'g', 'LineWidth',2);

title 'Horizontal Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('GPS test statistic','Threshold','HPL GPS (m)', 'HPE GPS (m)');
axis([30 120 0 30]);

%PLOT VERTICAL TEST STATS


figure();
hold
plot(((lambda_ss_outPLOTGPS_V_C(startepochplot:endepochplot,SubToPlot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((TD_outPLOT_V_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not

plot(VPL_GPS(startepochplot:endepochplot),'b', 'LineWidth',2);
plot(VPE_GPS_C(startepochplot:endepochplot),'g', 'LineWidth',2);

title 'Vertical Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('GPS test statistic','Threshold','VPL GPS (m)', 'VPE GPS (m)');
axis([30 120 0 40]);





startepochplot = startepoch+2;  %start it at i = 5
endepochplot = endepoch;

SubToPlot = MaxSlopeSatelliteToDo;


%calculate the time when things cross

%errorstarted at i = 30


doneit1_INS_H = 0;
doneit2_INS_H = 0;

doneit1_INS_V = 0;
doneit2_INS_V = 0;


doneit1_MODEL_H = 0;
doneit2_MODEL_H = 0;

doneit1_MODEL_V = 0;
doneit2_MODEL_V = 0;


doneit1_GPS_H = 0;
doneit2_GPS_H = 0;


doneit1_GPS_V = 0;
doneit2_GPS_V = 0;




for pp = startepochplot:endepochplot
    
    
    %IMU H
    if (lambda_ss_outPLOTINS_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_INS_H == 0        
        TimeDetected_INS_H = pp;        
        doneit1_INS_H = 1;        
    end
    
     if HPE_INS_C(pp) > HPL_INS_C(pp) && doneit2_INS_H == 0        
        TimePECrossPL_INS_H = pp;        
        doneit2_INS_H = 1;        
     end      
     
     
     %IMU V
    if (lambda_ss_outPLOTINS_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_INS_V == 0        
        TimeDetected_INS_V = pp;        
        doneit1_INS_V = 1;        
    end
    
    
     if VPE_INS_C(pp) > VPL_INS_C(pp) && doneit2_INS_V == 0        
        TimePECrossPL_INS_V = pp;        
        doneit2_INS_V = 1;        
     end
    
     
     
     
     
         %ADM H
    if (lambda_ss_outPLOTMODEL_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_MODEL_H == 0        
        TimeDetected_MODEL_H = pp;        
        doneit1_MODEL_H = 1;        
    end
    
     if HPE_MODEL_C(pp) > HPL_MODEL_C(pp) && doneit2_MODEL_H == 0        
        TimePECrossPL_MODEL_H = pp;        
        doneit2_MODEL_H = 1;        
     end      
     
     
     %ADM V
    if (lambda_ss_outPLOTMODEL_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_MODEL_V == 0        
        TimeDetected_MODEL_V = pp;        
        doneit1_MODEL_V = 1;        
    end
    
    
     if VPE_MODEL_C(pp) > VPL_MODEL_C(pp) && doneit2_MODEL_V == 0        
        TimePECrossPL_MODEL_V = pp;        
        doneit2_MODEL_V = 1;        
     end
    
     
     
     
         %GPS H
    if (lambda_ss_outPLOTGPS_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_GPS_H == 0        
        TimeDetected_GPS_H = pp;        
        doneit1_GPS_H = 1;        
    end
    
     if HPE_GPS_C(pp) > HPL_GPS(pp) && doneit2_GPS_H == 0        
        TimePECrossPL_GPS_H = pp;        
        doneit2_GPS_H = 1;        
     end      
     
     
     %GPS V
    if (lambda_ss_outPLOTGPS_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_GPS_V == 0        
        TimeDetected_GPS_V = pp;        
        doneit1_GPS_V = 1;        
    end
    
    
     if VPE_GPS_C(pp) > VPL_GPS(pp) && doneit2_GPS_V == 0        
        TimePECrossPL_GPS_V = pp;        
        doneit2_GPS_V = 1;        
     end
    
    
     
        
        
end
    





%with correction for when error started ,at time 30


TimeDetected_INS_H = TimeDetected_INS_H - 30;
TimePECrossPL_INS_H = TimePECrossPL_INS_H - 30;

    TimeDetected_INS_V = TimeDetected_INS_V - 30;
    TimePECrossPL_INS_V = TimePECrossPL_INS_V - 30;
    


TimeDetected_MODEL_H = TimeDetected_MODEL_H - 30;
TimePECrossPL_MODEL_H = TimePECrossPL_MODEL_H - 30;

    TimeDetected_MODEL_V = TimeDetected_MODEL_V - 30;
    TimePECrossPL_MODEL_V = TimePECrossPL_MODEL_V   - 30; 
    
    

TimeDetected_GPS_H = TimeDetected_GPS_H - 30;
TimePECrossPL_GPS_H = TimePECrossPL_GPS_H - 30;

TimeDetected_GPS_V = TimeDetected_GPS_V - 30;
TimePECrossPL_GPS_V = TimePECrossPL_GPS_V - 30;



TimeDetected_INS_H
TimePECrossPL_INS_H

TimeDetected_INS_V
TimePECrossPL_INS_V    


TimeDetected_MODEL_H
TimePECrossPL_MODEL_H

TimeDetected_MODEL_V
TimePECrossPL_MODEL_V
    
    
    
TimeDetected_GPS_H
TimePECrossPL_GPS_H

TimeDetected_GPS_V
TimePECrossPL_GPS_V
    
















%plot of position error

figure();
hold
plot(((HPE_INS_C(startepochplot:endepochplot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((HPL_INS_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not


title 'Horizontal Position Error' ; xlabel('Time (s)');  ylabel('Position Error (m)');
legend('HPE IMU','HPL IMU');
axis([0 endepochplot-startepochplot 0 10]);



figure();
hold
plot(((VPE_INS_C(startepochplot:endepochplot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((VPL_INS_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not


title 'Vertical Position Error' ; xlabel('Time (s)');  ylabel('Position Error (m)');
legend('VPE IMU','VPL IMU');
axis([0 endepochplot-startepochplot 0 20]);







figure();
hold
plot(((HPE_MODEL_C(startepochplot:endepochplot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((HPL_MODEL_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not


title 'Horizontal Position Error' ; xlabel('Time (s)');  ylabel('Position Error (m)');
legend('HPE ADM','HPL ADM');
axis([0 endepochplot-startepochplot 0 20]);



figure();
hold
plot(((VPE_MODEL_C(startepochplot:endepochplot))),'r','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases

plot((VPL_MODEL_C(startepochplot:endepochplot)),'k','LineWidth',2);    %not sure if this is right or not


title 'Vertical Position Error' ; xlabel('Time (s)');  ylabel('Position Error (m)');
legend('VPE ADM','VPL ADM');
axis([0 endepochplot-startepochplot 0 50]);






figure();
%plot((sqrt(lambda_ss_outPLOT_H(startepoch:endepoch,:))),'g-'); title 'Horizontal Test statistic vs time(m)';xlabel('Time');
hold

%plot((sqrt(lambda_ss_outPLOTINS_H(startepoch:endepoch,:))),'b-');
%plot((sqrt(lambda_ss_outPLOTINS_H(startepoch:endepoch,:))),'r-');

plot(((lambda_ss_outPLOTGPS_H_C(startepoch:endepoch,:))),'r-');

%plot((HPL(startepoch:endepoch)),'g--','LineWidth',2);
%plot(HPL_INS(startepoch:endepoch),'b--','LineWidth',2);
plot(HPL_GPS(startepoch:endepoch),'r--','LineWidth',2);

plot((TD_outPLOT_H_C(startepoch:endepoch)),'k--','LineWidth',2);    %not sure if this is right or not





















%first tune the filters individually, ins, adm, then
%see/check that the fused one is consistent.

%Tuning Parameters for Thesis Results
%29 August 
%when tuning, considered the steady state period only
%this configuration seems to give good tuning results for the steady state
% The HPE and VPE for the fused is better than the unfused, when its
% reached steady state too which is good...
% 
% 
% WINS = 
%  8.55670332180864e-005
%      6.16850275068085e-005
%      5.89739324707068e-005
%                   0.000169
%                   0.000324
%                     0.0001
% 
% 
%                     
%                     
%                     
% Wa = 0.133023255813953
%         0.0315926829268293        
%         0.0882632956521739        
%          0.495270588235294         
%                     3.5948                  
%                   104.7116                
% 
% 
% R  = 6.25
%                       6.25
%                       6.25
%                       6.25
%                       6.25
%                       6.25
%                       6.25
%                       6.25
%                          1
%                          1
%                          1
%                          1
%                          1
%                          1
%                          1
%                          1 


%for MMF:

% mean(HPLImprovepercent(20:100))
% 
% ans =
% 
%          0.671359093513095
% 
% >> mean(VPLImprovepercent(20:100))
% 
% ans =
% 
%           1.00033626125258




%Try reducing the variance (104) on the z accelerometer. 
%I think the R values are OK, probably leave them. 

%29.8.08
%When the bias is turned on for the ADM the uncertainty increases and the
%bias goes large and doesnt stay bounded...
%perhaps i need to model the p q r and accelerations differently, still
%have 6 states but they are built up from the equations involving the ADM
%parameters which are modelled as GM processes... ... 





%with higher PR noise in R



figure();
hold;
%plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(HPL_INS_NF(startepochplot:endepochplot),'k--','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot(HPL_INS_C(startepochplot:endepochplot),'k','LineWidth',2);  
%plot(HPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

%legend('HAL','HPL IMU','HPL IMU-ADM', 'HPL GPS');       


legend('HPL IMU','HPL IMU-ADM');    

axis([0 endepochplot-startepochplot 12 16]);


figure();
hold;
%plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(VPL_INS_NF(startepochplot:endepochplot),'k--','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot(VPL_INS_C(startepochplot:endepochplot),'k','LineWidth',2);  
%plot(VPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

%legend('VAL','VPL IMU','VPL IMU-ADM', 'VPL GPS');       


legend('VPL IMU','VPL IMU-ADM');    

axis([0 endepochplot-startepochplot 20 30]);






figure();
hold;
%plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(HPL_MODEL_NF(startepochplot:endepochplot),'k--','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot(HPL_MODEL_C(startepochplot:endepochplot),'k','LineWidth',2);  
%plot(HPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

%legend('HAL','HPL IMU','HPL IMU-ADM', 'HPL GPS');       


legend('HPL ADM','HPL IMU-ADM');    

axis([0 endepochplot-startepochplot 8 13]);


figure();
hold;
%plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(VPL_MODEL_NF(startepochplot:endepochplot),'k--','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot(VPL_MODEL_C(startepochplot:endepochplot),'k','LineWidth',2);  
%plot(VPL_GPS(startepochplot:endepochplot),'g','LineWidth',2); 
%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

%legend('VAL','VPL IMU','VPL IMU-ADM', 'VPL GPS');       


legend('VPL ADM','VPL IMU-ADM');    

axis([0 endepochplot-startepochplot 11 15]);





%to calculate improvement in accuracy




dmajorimprovementINS = dmajor_INSNF(startepoch:endepoch) - dmajor_INS(startepoch:endepoch);

improvement_sigma_v_INS = sigma_v_INSNF(startepoch:endepoch) - sigma_v_INS(startepoch:endepoch);




dmajorImprovepercentINS = (dmajorimprovementINS./dmajor_INSNF(startepoch:endepoch))*100;
sigma_v_ImprovepercentINS = (improvement_sigma_v_INS./sigma_v_INS(startepoch:endepoch))*100;


mean(dmajorImprovepercentINS(startepochplot:endepochplot-2))
mean(sigma_v_ImprovepercentINS(startepochplot:endepochplot-2))









figure;
hold;
title 'Horizontal Accuracy with MMF' ; xlabel('Time (s)');
plot(dmajor_INS(startepochplot:endepochplot),'b','LineWidth',2);
plot(dmajor_INSNF(startepochplot:endepochplot),'b','LineWidth',2);
plot(dmajor_MNF(startepochplot:endepochplot),'r','LineWidth',2);
plot(dmajor_GPS(startepochplot:endepochplot),'g','LineWidth',2);
plot(dmajor_max_INSNF(startepochplot:endepochplot),'k','LineWidth',4);
legend('IMU', 'IMU MMF', 'ADM','GPS');   
axis([0 endepochplot-startepochplot 0 3]);











%plots for pqr2 and pqr 3


%plot of the true wind conditions






figure;
hold;
title 'Airspeed Measurement' ; xlabel('Time (s)'); ylabel('Airspeed (m/s)');
plot(Airspeed_truth1Hz(startepochplot:endepochplot),'r','LineWidth',2);
plot(Airspeed_Sensor(startepochplot:endepochplot),'b','LineWidth',2);
legend('Airspeed Truth', 'Airspeed Measured');   


figure;
hold;
title 'Angle of Attack Measurement' ; xlabel('Time (s)');  ylabel('Angle of Attack (deg)'); 
plot(Alpha_truth1Hz(startepochplot:endepochplot)*180/pi,'r','LineWidth',2);
plot(Alpha_Sensor(startepochplot:endepochplot)*180/pi,'b','LineWidth',2);
legend('Angle of Attack Truth', 'Angle of Attack Measured');   


figure;
hold;
title 'Angle of Sideslip Measurement' ; xlabel('Time (s)');  ylabel('Angle of Sideslip (deg)');
plot(Beta_truth1Hz(startepochplot:endepochplot)*180/pi,'r','LineWidth',2);
plot(BDeta_Sensor(startepochplot:endepochplot)*180/pi,'b','LineWidth',2);
legend('Angle of Sideslip Truth', 'Angle of Sideslip Measured');   





%aircraft attitude plots

figure;
hold;
title 'Aircraft Attitude' ; xlabel('Time (s)');  ylabel('Attitude (deg)');
plot(Roll_truth1Hz(startepochplot:endepochplot)*180/pi,'b','LineWidth',2);
plot(Pitch_truth1Hz(startepochplot:endepochplot)*180/pi,'g','LineWidth',2);
%plot(Yaw_truth1Hz(startepochplot:endepochplot)*180/pi,'r','LineWidth',2);

legend('Roll', 'Pitch', 'Yaw');   









%PLOTS for Journal


%calculate average HPL improvementss for Monte Carlo stuff...




HPLImprovement = (mean(HPL_INS_NFmc(:,startepochplot:endepochplot))) -(mean(HPL_INS_Cmc(:,startepochplot:endepochplot)));
VPLImprovement = (mean(VPL_INS_NFmc(:,startepochplot:endepochplot))) -(mean(VPL_INS_Cmc(:,startepochplot:endepochplot)));


HPLImprovepercent = (HPLImprovement./(mean(HPL_INS_NFmc(:,startepochplot:endepochplot))))*100;
VPLImprovepercent = (VPLImprovement./(mean(VPL_INS_NFmc(:,startepochplot:endepochplot))))*100;


mean(HPLImprovepercent)
mean(VPLImprovepercent)




HPLImprovement = (mean(HPL_MODEL_NFmc(:,startepochplot:endepochplot))) -(mean(HPL_MODEL_Cmc(:,startepochplot:endepochplot)));
VPLImprovement = (mean(VPL_MODEL_NFmc(:,startepochplot:endepochplot))) -(mean(VPL_MODEL_Cmc(:,startepochplot:endepochplot)));


HPLImprovepercent = (HPLImprovement./(mean(HPL_MODEL_NFmc(:,startepochplot:endepochplot))))*100;
VPLImprovepercent = (VPLImprovement./(mean(VPL_MODEL_NFmc(:,startepochplot:endepochplot))))*100;


mean(HPLImprovepercent)
mean(VPLImprovepercent)




%for percent improvement of IMU over GPS:




HPLImprovement = (mean(HPL_GPS_mc(:,startepochplot:endepochplot))) -(mean(HPL_INS_NFmc(:,startepochplot:endepochplot)));
VPLImprovement = (mean(VPL_GPS_mc(:,startepochplot:endepochplot))) -(mean(VPL_INS_NFmc(:,startepochplot:endepochplot)));


HPLImprovepercent = (HPLImprovement./(mean(HPL_GPS_mc(:,startepochplot:endepochplot))))*100;
VPLImprovepercent = (VPLImprovement./(mean(VPL_GPS_mc(:,startepochplot:endepochplot))))*100;


mean(HPLImprovepercent)
mean(VPLImprovepercent)


%mean percent improvement of ADM over GPS. 

HPLImprovement = (mean(HPL_GPS_mc(:,startepochplot:endepochplot))) -(mean(HPL_MODEL_NFmc(:,startepochplot:endepochplot)));
VPLImprovement = (mean(VPL_GPS_mc(:,startepochplot:endepochplot))) -(mean(VPL_MODEL_NFmc(:,startepochplot:endepochplot)));


HPLImprovepercent = (HPLImprovement./(mean(HPL_GPS_mc(:,startepochplot:endepochplot))))*100;
VPLImprovepercent = (VPLImprovement./(mean(VPL_GPS_mc(:,startepochplot:endepochplot))))*100;


mean(HPLImprovepercent)
mean(VPLImprovepercent)

%mean percent improvement of IMU over ADM. 





HPLImprovement = (mean(HPL_MODEL_NFmc(:,startepochplot:endepochplot))) -(mean(HPL_INS_NFmc(:,startepochplot:endepochplot)));
VPLImprovement = (mean(VPL_MODEL_NFmc(:,startepochplot:endepochplot))) -(mean(VPL_INS_NFmc(:,startepochplot:endepochplot)));


HPLImprovepercent = (HPLImprovement./(mean(HPL_MODEL_NFmc(:,startepochplot:endepochplot))))*100;
VPLImprovepercent = (VPLImprovement./(mean(VPL_MODEL_NFmc(:,startepochplot:endepochplot))))*100;

mean(HPLImprovepercent)
mean(VPLImprovepercent)








%For journal (with using no color etc)
%plot means 


figure();
hold;
plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(mean(HPL_INS_NF(startepochplot:endepochplot),'k--','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot(HPL_INS_C(startepochplot:endepochplot),'k','LineWidth',2);  
%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('HAL','HPL IMU','HPL IMU-ADM');       

axis([0 endepochplot-startepochplot 11 41]);





figure();
hold;
plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(VPL_INS_NF(startepochplot:endepochplot),'k--','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot(VPL_INS_C(startepochplot:endepochplot),'k','LineWidth',2);  
%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('VAL','VPL IMU','VPL IMU-ADM');       

axis([0 endepochplot-startepochplot 16 51]);





















%calculate the time when things cross, for Monte Carlo stuff!!!

%errorstarted at i = 30



startepochplot = startepoch+2;  %start it at i = 5
endepochplot = endepoch;

for iM = 1:Nmc
    
       
%     
% 
% TD_outPLOT_H_C = sqrt(aTDINS_H);
%     
%     lambda_ss_outPLOTINS_H_C = sqrt(lambda_ss_outINS_H_NFmc(pp,:,:);
% 
% %lambda_ss_outPLOTINS_H = lambda_ss_outINS_H;  %converting from radians to metres
% 
% 
% 
% 
% TD_outPLOT_V_C = sqrt(aTDINS_V);
%     
% lambda_ss_outPLOTINS_V_C = sqrt(lambda_ss_outINS_V);  
%     
%     
%     
    
      
      
%iM = 29;



%for UNFUSED Results


TD_outPLOT_H_C = sqrt(aTDINS_Hmc(iM,:));
    
    lambda_ss_outPLOTINS_H_C(:,:) = sqrt(lambda_ss_outINS_H_NFmc(iM,:,:));


TD_outPLOT_V_C = sqrt(aTDINS_Vmc(iM,:));
    
    lambda_ss_outPLOTINS_V_C(:,:)= sqrt(lambda_ss_outINS_V_NFmc(iM,:,:));
    
    
    
      
    

    
    lambda_ss_outPLOTMODEL_H_C(:,:) = sqrt(lambda_ss_outMODEL_H_NFmc(iM,:,:));    
    lambda_ss_outPLOTMODEL_V_C(:,:) = sqrt(lambda_ss_outMODEL_V_NFmc(iM,:,:));         
        

    
    lambda_ss_outPLOTGPS_H_C(:,:) = sqrt(lambda_ss_outGPS_H_NFmc(iM,:,:));    
    lambda_ss_outPLOTGPS_V_C(:,:) = sqrt(lambda_ss_outGPS_V_NFmc(iM,:,:));
    




      
      
    
    
    
    SubToPlot = MaxSlopeSatelliteToDo_mc(iM,1);

  %  SubToPlot = 12;
    

doneit1_INS_H = 0;
doneit2_INS_H = 0;

doneit1_INS_V = 0;
doneit2_INS_V = 0;


doneit1_MODEL_H = 0;
doneit2_MODEL_H = 0;

doneit1_MODEL_V = 0;
doneit2_MODEL_V = 0;


doneit1_GPS_H = 0;
doneit2_GPS_H = 0;


doneit1_GPS_V = 0;
doneit2_GPS_V = 0;











for pp = startepochplot:endepochplot
    
    
    
    
    
    
    %IMU H
    if (lambda_ss_outPLOTINS_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_INS_H == 0        
        TimeDetected_INS_H = pp;        
        doneit1_INS_H = 1;        
    end
    
     if HPE_INS_NFmc(iM,pp) > HPL_INS_NFmc(iM,pp) && doneit2_INS_H == 0        
        TimePECrossPL_INS_H = pp;        
        doneit2_INS_H = 1;        
     end      
     
     
     %IMU V
    if (lambda_ss_outPLOTINS_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_INS_V == 0        
        TimeDetected_INS_V = pp;        
        doneit1_INS_V = 1;        
    end
    
    
     if VPE_INS_NFmc(iM,pp) > VPL_INS_NFmc(iM,pp) && doneit2_INS_V == 0        
        TimePECrossPL_INS_V = pp;        
        doneit2_INS_V = 1;        
     end
    
     
     
     
     
      
        lambda_k = (lambda_ss_outPLOTINS_H_C(pp,SubToPlot)^2); %noncentrality parameter is max of test statistics. %removed max in case its false alarms, i specifiy the value
    %threshold squared
   TD_squared = TD_outPLOT_H_C(pp)^2;
    N = N_save_Cmc(iM,pp);
    DOF = 2;  %for HPL
  sigma_re = sqrt((sqrt(P_save_NFmc(iM,7,7,pp)))^2 + (sqrt(P_save_NFmc(iM,8,8,pp)))^2)*ReConst;
    PL = HPL_INS_NFmc(iM,pp);
    de = HPE_INS_NFmc(iM,pp);        
    
    [Pmd_Est_INS_Hmc(iM,pp), PrHPE_INS_Hmc(iM,pp), PrTestStat_INS_Hmc(iM,pp)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
       
    
    
        lambda_k = (lambda_ss_outPLOTINS_V_C(pp,SubToPlot)^2); %noncentrality parameter is max of test statistics. %removed max in case its false alarms, i specifiy the value
    %threshold squared
   TD_squared = TD_outPLOT_V_C(pp)^2;
    N = N_save_Cmc(iM,pp);
    DOF = 1;  %for VPL
  sigma_re = sqrt(P_save_NFmc(iM,9,9,pp));
    PL = VPL_INS_NFmc(iM,pp);
    de = VPE_INS_NFmc(iM,pp);        
    
    [Pmd_Est_INS_Vmc(iM,pp), PrHPE_INS_Vmc(iM,pp), PrTestStat_INS_Vmc(iM,pp)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
       
    
     
     
     
     
     
     
     
     
     
     
     
     
     
         %ADM H
    if (lambda_ss_outPLOTMODEL_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_MODEL_H == 0        
        TimeDetected_MODEL_H = pp;        
        doneit1_MODEL_H = 1;        
    end
    
     if HPE_MODEL_NFmc(iM,pp) > HPL_MODEL_NFmc(iM,pp) && doneit2_MODEL_H == 0        
        TimePECrossPL_MODEL_H = pp;        
        doneit2_MODEL_H = 1;        
     end      
     
     
     %ADM V
    if (lambda_ss_outPLOTMODEL_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_MODEL_V == 0        
        TimeDetected_MODEL_V = pp;        
        doneit1_MODEL_V = 1;        
    end
    
    
     if VPE_MODEL_NFmc(iM,pp) > VPL_MODEL_NFmc(iM,pp) && doneit2_MODEL_V == 0        
        TimePECrossPL_MODEL_V = pp;        
        doneit2_MODEL_V = 1;        
     end
    
     
     
     
         %GPS H
    if (lambda_ss_outPLOTGPS_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_GPS_H == 0        
        TimeDetected_GPS_H = pp;        
        doneit1_GPS_H = 1;        
    end
    
     if HPE_GPS_mc(iM,pp) > HPL_GPS_mc(iM,pp) && doneit2_GPS_H == 0        
        TimePECrossPL_GPS_H = pp;        
        doneit2_GPS_H = 1;        
     end      
     
     
     %GPS V
    if (lambda_ss_outPLOTGPS_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_GPS_V == 0        
        TimeDetected_GPS_V = pp;        
        doneit1_GPS_V = 1;        
    end
    
    
     if VPE_GPS_mc(iM,pp) > VPL_GPS_mc(iM,pp) && doneit2_GPS_V == 0        
        TimePECrossPL_GPS_V = pp;        
        doneit2_GPS_V = 1;        
     end
    
    
        
        
end
    


TimeDetected_INS_H_mc(iM) = TimeDetected_INS_H;
TimePECrossPL_INS_H_mc(iM) = TimePECrossPL_INS_H;

TimeDetected_INS_V_mc(iM) = TimeDetected_INS_V;
TimePECrossPL_INS_V_mc(iM) = TimePECrossPL_INS_V;


TimeDetected_MODEL_H_mc(iM) = TimeDetected_MODEL_H;
TimePECrossPL_MODEL_H_mc(iM) = TimePECrossPL_MODEL_H;

TimeDetected_MODEL_V_mc(iM) = TimeDetected_MODEL_V;
TimePECrossPL_MODEL_V_mc(iM) = TimePECrossPL_MODEL_V;
    
    
    
TimeDetected_GPS_H_mc(iM) = TimeDetected_GPS_H;
TimePECrossPL_GPS_H_mc(iM) = TimePECrossPL_GPS_H;

TimeDetected_GPS_V_mc(iM) = TimeDetected_GPS_V;
TimePECrossPL_GPS_V_mc(iM) = TimePECrossPL_GPS_V;
    
     
     
     



%with correction for when error started ,at time 30
% 
% 
% TimeDetected_INS_H = TimeDetected_INS_H - 30;
% TimePECrossPL_INS_H = TimePECrossPL_INS_H - 30;
% 
%     TimeDetected_INS_V = TimeDetected_INS_V - 30;
%     TimePECrossPL_INS_V = TimePECrossPL_INS_V - 30;
%     
% 
% 
% TimeDetected_MODEL_H = TimeDetected_MODEL_H - 30;
% TimePECrossPL_MODEL_H = TimePECrossPL_MODEL_H - 30;
% 
%     TimeDetected_MODEL_V = TimeDetected_MODEL_V - 30;
%     TimePECrossPL_MODEL_V = TimePECrossPL_MODEL_V   - 30; 
%     
%     
% 
% TimeDetected_GPS_H = TimeDetected_GPS_H - 30;
% TimePECrossPL_GPS_H = TimePECrossPL_GPS_H - 30;
% 
% TimeDetected_GPS_V = TimeDetected_GPS_V - 30;
% TimePECrossPL_GPS_V = TimePECrossPL_GPS_V - 30;



end %for pp = 1:iM






TimeDetected_INS_H
TimePECrossPL_INS_H

TimeDetected_INS_V
TimePECrossPL_INS_V    


TimeDetected_MODEL_H
TimePECrossPL_MODEL_H

TimeDetected_MODEL_V
TimePECrossPL_MODEL_V
    
    
    
TimeDetected_GPS_H
TimePECrossPL_GPS_H

TimeDetected_GPS_V
TimePECrossPL_GPS_V
    
% 









%USE THIS CODE HERE FOR THE FUSED ie MMF!! performance stats:





startepochplot = startepoch+2;  %start it at i = 5
endepochplot = endepoch;

for iM = 1:Nmc
    
       

      
      
iM = 29;



%iM = 69;



%for MMF


TD_outPLOT_H_C = sqrt(aTDINS_Hmc(iM,:));
    
    lambda_ss_outPLOTINS_H_C(:,:) = sqrt(lambda_ss_outINS_Hmc(iM,:,:));


TD_outPLOT_V_C = sqrt(aTDINS_Vmc(iM,:));
    
    lambda_ss_outPLOTINS_V_C(:,:)= sqrt(lambda_ss_outINS_Vmc(iM,:,:));
    
    
    
      
    

    
    lambda_ss_outPLOTMODEL_H_C(:,:) = sqrt(lambda_ss_outMODEL_Hmc(iM,:,:));    
    lambda_ss_outPLOTMODEL_V_C(:,:) = sqrt(lambda_ss_outMODEL_Vmc(iM,:,:));         
        

    
    lambda_ss_outPLOTGPS_H_C(:,:) = sqrt(lambda_ss_outGPS_Hmc(iM,:,:));    
    lambda_ss_outPLOTGPS_V_C(:,:) = sqrt(lambda_ss_outGPS_Vmc(iM,:,:));
    




      
      
    
    
    
    SubToPlot = MaxSlopeSatelliteToDo_mc(iM,1);

  %  SubToPlot = 12;
    

doneit1_INS_H = 0;
doneit2_INS_H = 0;

doneit1_INS_V = 0;
doneit2_INS_V = 0;


doneit1_MODEL_H = 0;
doneit2_MODEL_H = 0;

doneit1_MODEL_V = 0;
doneit2_MODEL_V = 0;


doneit1_GPS_H = 0;
doneit2_GPS_H = 0;


doneit1_GPS_V = 0;
doneit2_GPS_V = 0;











for pp = startepochplot:endepochplot
    
    
    
    
    
    %IMU H
    if (lambda_ss_outPLOTINS_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_INS_H == 0        
        TimeDetected_INS_H = pp;        
        doneit1_INS_H = 1;        
    end
    
     if HPE_INS_Cmc(iM,pp) > HPL_INS_Cmc(iM,pp) && doneit2_INS_H == 0        
        TimePECrossPL_INS_H = pp;        
        doneit2_INS_H = 1;        
     end      
     
     
     
    
    
    
     
     
      
   lambda_k = (lambda_ss_outPLOTINS_H_C(pp,SubToPlot)^2); %noncentrality parameter is max of test statistics. %removed max in case its false alarms, i specifiy the value
    %threshold squared
   TD_squared = TD_outPLOT_H_C(pp)^2;
    N = N_save_Cmc(iM,pp);
    DOF = 2;  %for HPL
  sigma_re = sqrt((sqrt(P_save_Cmc(iM,7,7,pp)))^2 + (sqrt(P_save_Cmc(iM,8,8,pp)))^2)*ReConst;
    PL = HPL_INS_Cmc(iM,pp);
    de = HPE_INS_Cmc(iM,pp);        
    
    [Pmd_Est_INS_Hmc_C(iM,pp), PrHPE_INS_Hmc_C(iM,pp), PrTestStat_INS_Hmc_C(iM,pp)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
       
    
    
        lambda_k = (lambda_ss_outPLOTINS_V_C(pp,SubToPlot)^2); %noncentrality parameter is max of test statistics. %removed max in case its false alarms, i specifiy the value
    %threshold squared
   TD_squared = TD_outPLOT_V_C(pp)^2;
    N = N_save_Cmc(iM,pp);
    DOF = 1;  %for VPL
  sigma_re = sqrt(P_save_Cmc(iM,9,9,pp));
    PL = VPL_INS_Cmc(iM,pp);
    de = VPE_INS_Cmc(iM,pp);        
    
    [Pmd_Est_INS_Vmc_C(iM,pp), PrHPE_INS_Vmc_C(iM,pp), PrTestStat_INS_Vmc_C(iM,pp)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
       
    
    
    
     
     
     
     
     
     
     
     
     
     
     
     
     
     %IMU V
    if (lambda_ss_outPLOTINS_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_INS_V == 0        
        TimeDetected_INS_V = pp;        
        doneit1_INS_V = 1;        
    end
    
    
     if VPE_INS_Cmc(iM,pp) > VPL_INS_Cmc(iM,pp) && doneit2_INS_V == 0        
        TimePECrossPL_INS_V = pp;        
        doneit2_INS_V = 1;        
     end
    
     
     
     
     
         %ADM H
    if (lambda_ss_outPLOTMODEL_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_MODEL_H == 0        
        TimeDetected_MODEL_H = pp;        
        doneit1_MODEL_H = 1;        
    end
    
     if HPE_MODEL_Cmc(iM,pp) > HPL_MODEL_Cmc(iM,pp) && doneit2_MODEL_H == 0        
        TimePECrossPL_MODEL_H = pp;        
        doneit2_MODEL_H = 1;        
     end      
     
     
     %ADM V
    if (lambda_ss_outPLOTMODEL_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_MODEL_V == 0        
        TimeDetected_MODEL_V = pp;        
        doneit1_MODEL_V = 1;        
    end
    
    
     if VPE_MODEL_Cmc(iM,pp) > VPL_MODEL_Cmc(iM,pp) && doneit2_MODEL_V == 0        
        TimePECrossPL_MODEL_V = pp;        
        doneit2_MODEL_V = 1;        
     end
    
     
     
     
         %GPS H
    if (lambda_ss_outPLOTGPS_H_C(pp,SubToPlot) > TD_outPLOT_H_C(pp)) && doneit1_GPS_H == 0        
        TimeDetected_GPS_H = pp;        
        doneit1_GPS_H = 1;        
    end
    
     if HPE_GPS_mc(iM,pp) > HPL_GPS_mc(iM,pp) && doneit2_GPS_H == 0        
        TimePECrossPL_GPS_H = pp;        
        doneit2_GPS_H = 1;        
     end      
     
     
     %GPS V
    if (lambda_ss_outPLOTGPS_V_C(pp,SubToPlot) > TD_outPLOT_V_C(pp)) && doneit1_GPS_V == 0        
        TimeDetected_GPS_V = pp;        
        doneit1_GPS_V = 1;        
    end
    
    
     if VPE_GPS_mc(iM,pp) > VPL_GPS_mc(iM,pp) && doneit2_GPS_V == 0        
        TimePECrossPL_GPS_V = pp;        
        doneit2_GPS_V = 1;        
     end
    
    
        
        
end
    


TimeDetected_INS_H_mc(iM) = TimeDetected_INS_H;
TimePECrossPL_INS_H_mc(iM) = TimePECrossPL_INS_H;

TimeDetected_INS_V_mc(iM) = TimeDetected_INS_V;
TimePECrossPL_INS_V_mc(iM) = TimePECrossPL_INS_V;


TimeDetected_MODEL_H_mc(iM) = TimeDetected_MODEL_H;
TimePECrossPL_MODEL_H_mc(iM) = TimePECrossPL_MODEL_H;

TimeDetected_MODEL_V_mc(iM) = TimeDetected_MODEL_V;
TimePECrossPL_MODEL_V_mc(iM) = TimePECrossPL_MODEL_V;
    
    
    
TimeDetected_GPS_H_mc(iM) = TimeDetected_GPS_H;
TimePECrossPL_GPS_H_mc(iM) = TimePECrossPL_GPS_H;

TimeDetected_GPS_V_mc(iM) = TimeDetected_GPS_V;
TimePECrossPL_GPS_V_mc(iM) = TimePECrossPL_GPS_V;
    
     
     
     



end %for iM = 1:Nmc

figure;
hold;


%process
for iM = 1:Nmc
    
    
    plot(Pmd_Est_INS_Hmc(iM,1:TimeDetected_INS_H_mc(iM)));
    
    plot(Pmd_H,'r*')
    
 %plot(Pmd_Est_INS_Hmc(iM,1:40));
    
    
end



figure;
hold;



    plot(Pmd_Est_INS_Vmc(iM,1:TimeDetected_INS_V_mc(iM)),'g');
    
    plot(Pmd_V,'g*')
    
 %plot(Pmd_Est_INS_Hmc(iM,1:40));
    
 
 
    plot(PrHPE_INS_Vmc_C(iM,1:TimeDetected_INS_V_mc(iM)),'k');
       plot(PrTestStat_INS_Vmc_C(iM,1:TimeDetected_INS_V_mc(iM)),'b');




%process



%count how many geometries dont mean Pmd requirement


%for unfused

countabove_H = 0;
for iM = 1:Nmc
    
    for pp = startepochplot:TimeDetected_INS_H_mc(iM)
    
    

       
       
       if Pmd_Est_INS_Hmc(iM,pp) >Pmd_H
           
           countabove_H = countabove_H + 1;        
           break                     %put a break in here and it counts the total number of geometries when Pmd is not met..
           
       end
       
       
       end
      
end






countabove_V = 0;
for iM = 1:Nmc
    
    for pp = startepochplot:TimeDetected_INS_V_mc(iM)   
       
       
       if Pmd_Est_INS_Vmc(iM,pp) >Pmd_V
           
           countabove_V = countabove_V + 1;        
           
           
       end
       
       
       end
      
end


    




TimeDetected_INS_H
TimePECrossPL_INS_H

TimeDetected_INS_V
TimePECrossPL_INS_V    


TimeDetected_MODEL_H
TimePECrossPL_MODEL_H

TimeDetected_MODEL_V
TimePECrossPL_MODEL_V
    
    
    
TimeDetected_GPS_H
TimePECrossPL_GPS_H

TimeDetected_GPS_V
TimePECrossPL_GPS_V
    
% 
























%Plots for IEEE Journal , resubmitted version




%For journal (with using no color etc)
%plot means 




figure();
hold;
%plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(mean(HPL_GPS_mc(:,startepochplot:endepochplot)),'k-.','LineWidth',2);  
plot(mean(HPL_INS_NFmc(:,startepochplot:endepochplot)),'k--','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot(mean(HPL_INS_Cmc(:,startepochplot:endepochplot)),'k','LineWidth',2);  

%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('HPL GPS', 'HPL GPS-IMU EKF','HPL GPS-IMU-ADM EKF');       

axis([0 endepochplot-startepochplot 9 30]);





figure();
hold;
%plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(mean(VPL_GPS_mc(:,startepochplot:endepochplot)),'k-.','LineWidth',2);  
plot(mean(VPL_INS_NFmc(:,startepochplot:endepochplot)),'k--','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot(mean(VPL_INS_Cmc(:,startepochplot:endepochplot)),'k','LineWidth',2); 

%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);


legend('VPL GPS', 'VPL GPS-IMU EKF','VPL GPS-IMU-ADM EKF');          

axis([0 endepochplot-startepochplot 10 45]);






%Plot Worst case

%when there was an average of ..6 satellites...


VPL_GPS_mc

find((max(VPL_GPS_mc(:,100))))

[I,J] = find ((N_save_Cmc(:,100)==6))

MC runs when theres 6 sats:



 29
    55
    56
    57

    use 29
    
    
    
    
    
    
plot(N_save_Cmc(29,startepochplot:endepochplot),'k','LineWidth',2);  title 'Satellite Visibility' ; xlabel('Time (s)');  ylabel('Number of Satellites');
axis([0 endepochplot 0 10]);
    
    
    
figure();
hold;
plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot((HPL_GPS_mc(29,startepochplot:endepochplot)),'k-.','LineWidth',2);  
plot((HPL_INS_NFmc(29,startepochplot:endepochplot)),'k--','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot((HPL_INS_Cmc(29,startepochplot:endepochplot)),'k','LineWidth',2);  

%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('HAL', 'HPL GPS', 'HPL GPS-IMU EKF','HPL GPS-IMU-ADM EKF');       

axis([0 endepochplot-startepochplot 9 48]);





figure();
hold;
plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot((VPL_GPS_mc(29,startepochplot:endepochplot)),'k-.','LineWidth',2);  
plot((VPL_INS_NFmc(29,startepochplot:endepochplot)),'k--','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot((VPL_INS_Cmc(29,startepochplot:endepochplot)),'k','LineWidth',2); 

%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);


legend('VAL', 'VPL GPS', 'VPL GPS-IMU EKF','VPL GPS-IMU-ADM EKF');          

axis([0 endepochplot-startepochplot 6 80]);




%mean improvements for the worst case case



HPLImprovement = ((HPL_INS_NFmc(29,startepochplot:endepochplot))) -((HPL_INS_Cmc(29,startepochplot:endepochplot)));
VPLImprovement = ((VPL_INS_NFmc(29,startepochplot:endepochplot))) -((VPL_INS_Cmc(29,startepochplot:endepochplot)));


HPLImprovepercent = (HPLImprovement./((HPL_INS_NFmc(29,startepochplot:endepochplot))))*100;
VPLImprovepercent = (VPLImprovement./((VPL_INS_NFmc(29,startepochplot:endepochplot))))*100;


mean(HPLImprovepercent)
mean(VPLImprovepercent)








HPLImprovement = ((HPL_GPS_mc(29,startepochplot:endepochplot))) -((HPL_INS_NFmc(29,startepochplot:endepochplot)));
VPLImprovement = ((VPL_GPS_mc(29,startepochplot:endepochplot))) -((VPL_INS_NFmc(29,startepochplot:endepochplot)));


HPLImprovepercent = (HPLImprovement./((HPL_GPS_mc(29,startepochplot:endepochplot))))*100;
VPLImprovepercent = (VPLImprovement./((VPL_GPS_mc(29,startepochplot:endepochplot))))*100;


mean(HPLImprovepercent)
mean(VPLImprovepercent)





%plot test stats verse Threshold for subfilter 29



%NOTE!!! REMEMBER HAVE TO RUN THE PREVIOUS Montecarloi code so it already
%puts the right MC run in the TD_outPLOT_H_C etc...
%run this code first:

%       
% iM = 29;
% 
% 
% 
% %iM = 69;
% 
% 
% 
% %for MMF
% 
% 
% TD_outPLOT_H_C = sqrt(aTDINS_Hmc(iM,:));
%     
%     lambda_ss_outPLOTINS_H_C(:,:) = sqrt(lambda_ss_outINS_Hmc(iM,:,:));







SubToPlot = MaxSlopeSatelliteToDo_mc(29);

%SubToPlot = 12;





%NOTE to self - 14.12.08, I haven't included these figures in the journal
%paper. 

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

plot(VPL_INS_Cmc(29,startepochplot:endepochplot),'k', 'LineWidth',2);
plot(VPE_INS_Cmc(29,startepochplot:endepochplot),'k--', 'LineWidth',2);
plot(((lambda_ss_outPLOTINS_V_C(startepochplot:endepochplot,SubToPlot))),'k-.','LineWidth',2);  %even though fault is on SV 12, its the 9th one that increases
plot((TD_outPLOT_V_C(startepochplot:endepochplot)),'k+','LineWidth',2);    %not sure if this is right or not

title 'Vertical Test Statistic vs Threshold' ; xlabel('Time (s)');  ylabel('Test Statistic');
legend('VPL GPS-IMU-ADM EKF (m)', 'VPE GPS-IMU-ADM EKF (m)','Test statistic','Threshold');
axis([30 80 0 20]);






%for working out % of availability:





Unavailable1 = 0;
Unavailable2= 0;
Unavailable3= 0;
Unavailable4= 0;
Unavailable5= 0;




for iM = 1:Nmc
    
      
for pp = startepochplot:endepochplot 


    %IMU
    % unfused
    if HPL_INS_NFmc(iM,pp) > HAL(pp) || VPL_INS_NFmc(iM,pp) > VAL(pp)  
            Unavailable1 = Unavailable1 + 1;        
    end    
    
    %fused
    if HPL_INS_Cmc(iM,pp) > HAL(pp) || VPL_INS_Cmc(iM,pp) > VAL(pp)      
                Unavailable2 = Unavailable2 + 1;        
    end    
    
    
    
      %ADM
    % unfused
    if HPL_MODEL_NFmc(iM,pp) > HAL(pp) || VPL_MODEL_NFmc(iM,pp) > VAL(pp)  
            Unavailable3 = Unavailable3 + 1;        
    end    
    
    
    %fused
    if HPL_MODEL_Cmc(iM,pp) > HAL(pp) || VPL_MODEL_Cmc(iM,pp) > VAL(pp)      
                Unavailable4 = Unavailable4 + 1;        
    end    
    
    
    %GPS
      
      if HPL_GPS_mc(iM,pp) > HAL(pp) || VPL_GPS_mc(iM,pp) > VAL(pp)      
                Unavailable5 = Unavailable5 + 1;        
    end    
    
    
    
end


end


%so percentages of times unavailable is

totalgeometres = Nmc*(endepochplot-startepochplot);

percentunavailable = (Unavailable5/totalgeometres)*100;









%Plots in colour now for thesis, using geometryjournal.mat results..


%for mean results:



figure();
hold;
%plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(mean(HPL_GPS_mc(:,startepochplot:endepochplot)),'g','LineWidth',2);
plot(mean(HPL_MODEL_NFmc(:,startepochplot:endepochplot)),'r','LineWidth',2);
plot(mean(HPL_INS_NFmc(:,startepochplot:endepochplot)),'b','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot(mean(HPL_INS_Cmc(:,startepochplot:endepochplot)),'k','LineWidth',2);  

%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('HPL GPS', 'HPL ADM', 'HPL IMU','HPL IMU-ADM');       

axis([0 endepochplot-startepochplot 9 30]);



figure();
hold;
%plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot(mean(VPL_GPS_mc(:,startepochplot:endepochplot)),'g','LineWidth',2);
plot(mean(VPL_MODEL_NFmc(:,startepochplot:endepochplot)),'r','LineWidth',2);
plot(mean(VPL_INS_NFmc(:,startepochplot:endepochplot)),'b','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot(mean(VPL_INS_Cmc(:,startepochplot:endepochplot)),'k','LineWidth',2);  

%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('VPL GPS', 'VPL ADM', 'VPL IMU','VPL IMU-ADM');       

axis([0 endepochplot-startepochplot 10 48]);





%for worst case:




figure();
hold;
plot(HAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot((HPL_GPS_mc(29,startepochplot:endepochplot)),'g','LineWidth',2);
plot((HPL_MODEL_NFmc(29,startepochplot:endepochplot)),'r','LineWidth',2);
plot((HPL_INS_NFmc(29,startepochplot:endepochplot)),'b','LineWidth',2);  title 'Comparison of Horizontal Protection Levels' ; xlabel('Time (s)');  ylabel('Horizontal Protection Level (m)');
plot((HPL_INS_Cmc(29,startepochplot:endepochplot)),'k','LineWidth',2);  

%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('HAL', 'HPL GPS', 'HPL ADM', 'HPL IMU','HPL IMU-ADM');       

axis([0 endepochplot-startepochplot 9 55]);



figure();
hold;
plot(VAL(startepochplot:endepochplot), 'k','LineWidth',4);
%plot(AccuracyLineH(startepochplot:endepochplot), 'k','LineWidth',1);
plot((VPL_GPS_mc(29,startepochplot:endepochplot)),'g','LineWidth',2);
plot((VPL_MODEL_NFmc(29,startepochplot:endepochplot)),'r','LineWidth',2);
plot((VPL_INS_NFmc(29,startepochplot:endepochplot)),'b','LineWidth',2);  title 'Comparison of Vertical Protection Levels' ; xlabel('Time (s)');  ylabel('Vertical Protection Level (m)');
plot((VPL_INS_Cmc(29,startepochplot:endepochplot)),'k','LineWidth',2);  

%plot(HPL_INS_R(startepochplot:endepochplot),'k--','LineWidth',2); 
%plot(HPE_INS_C(startepochplot:endepochplot),'k','LineWidth',2);    %navigation system error
%plot(HPE_GPS_C(startepochplot:endepochplot),'g--','LineWidth',1);

legend('VAL', 'VPL GPS', 'VPL ADM', 'VPL IMU','VPL IMU-ADM');       

axis([0 endepochplot-startepochplot 10 90]);
















