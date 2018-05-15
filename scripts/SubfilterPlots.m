






numsub = 13;   %the subfilter we want to plot

%Plot NED position Errors for a certain subfilter

figure();
plot(N_PosErrorsINS_CtestSub(numsub,startepochplot:endepochplot),'b-','LineWidth',2); title 'North Position Error(m)';xlabel('Time');
hold;
% plot(N_PosErrorsMODEL,'r-', 'LineWidth',2);
% plot(N_PosErrors,'k-', 'LineWidth',2);
plot(N_PosErrorsGPS_C(startepochplot:endepochplot),'g-'); 

legend('EKF INS','EKF MODEL','EKF INS MODEL','LSQ GPS');

plot(0 + 2*stdev7pnINS_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev7pnINS_CSub(numsub,startepochplot:endepochplot),'b--');
% plot(0 + 2*stdev7pnMODEL(startepochplot:endepochplot),'r--'); 
% plot(0 - 2*stdev7pnMODEL(startepochplot:endepochplot),'r--');
% plot(0 + 2*stdev7pn(startepochplot:endepochplot),'k--');
% plot(0 - 2*stdev7pn(startepochplot:endepochplot),'k--');



figure();
plot(E_PosErrorsINS_CtestSub(numsub,startepochplot:endepochplot),'b-','LineWidth',2); title 'East Position Error(m)';xlabel('Time');
hold;
% plot(E_PosErrorsMODEL,'r-', 'LineWidth',2);
% plot(E_PosErrors,'k-', 'LineWidth',2);
plot(E_PosErrorsGPS_C(startepochplot:endepochplot),'g-'); 

legend('EKF INS','EKF MODEL','EKF INS MODEL','LSQ GPS');

plot(0 + 2*stdev8peINS_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev8peINS_CSub(numsub,startepochplot:endepochplot),'b--');
% plot(0 + 2*stdev8peMODEL(startepochplot:endepochplot),'r--');
% plot(0 - 2*stdev8peMODEL(startepochplot:endepochplot),'r--');
% plot(0 + 2*stdev8pe(startepochplot:endepochplot),'k--');
% plot(0 - 2*stdev8pe(startepochplot:endepochplot),'k--');

figure();
plot(D_PosErrorsINS_CtestSub(numsub,startepochplot:endepochplot),'b-','LineWidth',2); title 'Down Position Error(m)';xlabel('Time');
hold;
% plot(D_PosErrorsMODEL,'r-', 'LineWidth',2);
% plot(D_PosErrors,'k-', 'LineWidth',2);
plot(D_PosErrorsGPS_C(startepochplot:endepochplot),'g-'); 

legend('EKF INS','EKF MODEL','EKF INS MODEL','LSQ GPS');

plot(0 + 2*stdev9pdINS_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev9pdINS_CSub(numsub,startepochplot:endepochplot),'b--');
% plot(0 + 2*stdev9pdMODEL(startepochplot:endepochplot),'r--');
% plot(0 - 2*stdev9pdMODEL(startepochplot:endepochplot),'r--');
% plot(0 + 2*stdev9pd(startepochplot:endepochplot),'k--');
% plot(0 - 2*stdev9pd(startepochplot:endepochplot),'k--');

tilefigs;
pause;
close all;






%for model subfilters


figure();
plot(N_PosErrorsMODEL_CtestSub(numsub,startepochplot:endepochplot),'b-','LineWidth',2); title 'North Position Error(m)';xlabel('Time');
hold;
% plot(N_PosErrorsMODEL,'r-', 'LineWidth',2);
% plot(N_PosErrors,'k-', 'LineWidth',2);
plot(N_PosErrorsGPS_C(startepochplot:endepochplot),'g-'); 

legend('EKF INS','EKF MODEL','EKF INS MODEL','LSQ GPS');

plot(0 + 2*stdev7pnMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev7pnMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
% plot(0 + 2*stdev7pnMODEL(startepochplot:endepochplot),'r--'); 
% plot(0 - 2*stdev7pnMODEL(startepochplot:endepochplot),'r--');
% plot(0 + 2*stdev7pn(startepochplot:endepochplot),'k--');
% plot(0 - 2*stdev7pn(startepochplot:endepochplot),'k--');



figure();
plot(E_PosErrorsMODEL_CtestSub(numsub,startepochplot:endepochplot),'b-','LineWidth',2); title 'East Position Error(m)';xlabel('Time');
hold;
% plot(E_PosErrorsMODEL,'r-', 'LineWidth',2);
% plot(E_PosErrors,'k-', 'LineWidth',2);
plot(E_PosErrorsGPS_C(startepochplot:endepochplot),'g-'); 

legend('EKF INS','EKF MODEL','EKF INS MODEL','LSQ GPS');

plot(0 + 2*stdev8peMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev8peMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
% plot(0 + 2*stdev8peMODEL(startepochplot:endepochplot),'r--');
% plot(0 - 2*stdev8peMODEL(startepochplot:endepochplot),'r--');
% plot(0 + 2*stdev8pe(startepochplot:endepochplot),'k--');
% plot(0 - 2*stdev8pe(startepochplot:endepochplot),'k--');

figure();
plot(D_PosErrorsMODEL_CtestSub(numsub,startepochplot:endepochplot),'b-','LineWidth',2); title 'Down Position Error(m)';xlabel('Time');
hold;
% plot(D_PosErrorsMODEL,'r-', 'LineWidth',2);
% plot(D_PosErrors,'k-', 'LineWidth',2);
plot(D_PosErrorsGPS_C(startepochplot:endepochplot),'g-'); 

legend('EKF INS','EKF MODEL','EKF INS MODEL','LSQ GPS');

plot(0 + 2*stdev9pdMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev9pdMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
% plot(0 + 2*stdev9pdMODEL(startepochplot:endepochplot),'r--');
% plot(0 - 2*stdev9pdMODEL(startepochplot:endepochplot),'r--');
% plot(0 + 2*stdev9pd(startepochplot:endepochplot),'k--');
% plot(0 - 2*stdev9pd(startepochplot:endepochplot),'k--');

tilefigs;
pause;
close all;




%plot subfilters on top of full filter solution, blue is a specific
%subfilter

numsub = 22;

figure; hold;
 plot(N_PosErrorsINS_Ctest','r');

plot(N_PosErrorsINS_CtestSub','g');
plot(N_PosErrorsINS_CtestSub(numsub,:)','b');

figure; hold;
 plot(E_PosErrorsINS_Ctest','r');

plot(E_PosErrorsINS_CtestSub','g');
plot(E_PosErrorsINS_CtestSub(numsub,:)','b');


figure; hold;
 plot(D_PosErrorsINS_Ctest','r');

plot(D_PosErrorsINS_CtestSub','g');
plot(D_PosErrorsINS_CtestSub(numsub,:)','b');


tilefigs;
pause;
close all;




figure; hold;
 plot(N_PosErrorsMODEL_Ctest','r');

plot(N_PosErrorsMODEL_CtestSub','g');
plot(N_PosErrorsMODEL_CtestSub(numsub,:)','b');

figure; hold;
 plot(E_PosErrorsMODEL_Ctest','r');

plot(E_PosErrorsMODEL_CtestSub','g');
plot(E_PosErrorsMODEL_CtestSub(numsub,:)','b');


figure; hold;
 plot(D_PosErrorsMODEL_Ctest','r');

plot(D_PosErrorsMODEL_CtestSub','g');
plot(D_PosErrorsMODEL_CtestSub(numsub,:)','b');


tilefigs;
pause;
close all;












figure; hold;
 plot(N_PosErrorsINS_Ctest','r');

plot(N_PosErrorsINS_CtestSub','g');
plot(N_PosErrorsINS_CtestSub(numsub,:)','b');



plot(0 + 2*stdev7pnINS_C(:,startepochplot:endepochplot)','r--');
plot(0 - 2*stdev7pnINS_C(:,startepochplot:endepochplot)','r--');


plot(0 + 2*stdev7pnINS_CSub(:,startepochplot:endepochplot)','g--');
plot(0 - 2*stdev7pnINS_CSub(:,startepochplot:endepochplot)','g--');


plot(0 + 2*stdev7pnINS_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev7pnINS_CSub(numsub,startepochplot:endepochplot),'b--');




figure; hold;
 plot(E_PosErrorsINS_Ctest','r');

plot(E_PosErrorsINS_CtestSub','g');
plot(E_PosErrorsINS_CtestSub(numsub,:)','b');



plot(0 + 2*stdev8peINS_C(:,startepochplot:endepochplot)','r--');
plot(0 - 2*stdev8peINS_C(:,startepochplot:endepochplot)','r--');



plot(0 + 2*stdev8peINS_CSub(:,startepochplot:endepochplot)','g--');
plot(0 - 2*stdev8peINS_CSub(:,startepochplot:endepochplot)','g--');


plot(0 + 2*stdev8peINS_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev8peINS_CSub(numsub,startepochplot:endepochplot),'b--');



figure; hold;
 plot(D_PosErrorsINS_Ctest','r');

plot(D_PosErrorsINS_CtestSub','g');
plot(D_PosErrorsINS_CtestSub(numsub,:)','b');




plot(0 + 2*stdev9pdINS_C(:,startepochplot:endepochplot)','r--');
plot(0 - 2*stdev9pdINS_C(:,startepochplot:endepochplot)','r--');


plot(0 + 2*stdev9pdINS_CSub(:,startepochplot:endepochplot)','g--');
plot(0 - 2*stdev9pdINS_CSub(:,startepochplot:endepochplot)','g--');


plot(0 + 2*stdev9pdINS_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev9pdINS_CSub(numsub,startepochplot:endepochplot),'b--');



tilefigs;
pause;
close all;








figure; hold;
 plot(N_PosErrorsMODEL_Ctest','r');

plot(N_PosErrorsMODEL_CtestSub','g');
plot(N_PosErrorsMODEL_CtestSub(numsub,:)','b');


plot(0 + 2*stdev7pnMODEL_CSub(:,startepochplot:endepochplot)','g--');
plot(0 - 2*stdev7pnMODEL_CSub(:,startepochplot:endepochplot)','g--');


plot(0 + 2*stdev7pnMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev7pnMODEL_CSub(numsub,startepochplot:endepochplot),'b--');




figure; hold;
 plot(E_PosErrorsMODEL_Ctest','r');

plot(E_PosErrorsMODEL_CtestSub','g');
plot(E_PosErrorsMODEL_CtestSub(numsub,:)','b');

plot(0 + 2*stdev8peMODEL_CSub(:,startepochplot:endepochplot)','g--');
plot(0 - 2*stdev8peMODEL_CSub(:,startepochplot:endepochplot)','g--');


plot(0 + 2*stdev8peMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev8peMODEL_CSub(numsub,startepochplot:endepochplot),'b--');



figure; hold;
 plot(D_PosErrorsMODEL_Ctest','r');

plot(D_PosErrorsMODEL_CtestSub','g');
plot(D_PosErrorsMODEL_CtestSub(numsub,:)','b');

plot(0 + 2*stdev9pdMODEL_CSub(:,startepochplot:endepochplot)','g--');
plot(0 - 2*stdev9pdMODEL_CSub(:,startepochplot:endepochplot)','g--');


plot(0 + 2*stdev9pdMODEL_CSub(numsub,startepochplot:endepochplot),'b--');
plot(0 - 2*stdev9pdMODEL_CSub(numsub,startepochplot:endepochplot),'b--');



tilefigs;
pause;
close all;





%plot the actual difference between full and subfilter and covariance
%ellipse







for i = startepoch:endepoch

for kk = 1:N_save(i)
    
for n = Excluded_PRNSub(kk,i)
    
    
     P_save_subtempINSPlot(1:2,1) = P_out_CSubTotal(n,7:8,7,i) ;  %note, in radians lat
     P_save_subtempINSPlot(1:2,2) = P_out_CSubTotal(n,7:8,8,i) ;  %note, in radians lon
          


B_ssPlot(n,1:2,1:2,i) = (P_save_subtempINSPlot - P_save(7:8,7:8,i))*ReConst^2;



end
end
end



n = 8; 
i = 117;





figure; hold;



for i = 5:endepoch
    
 B_ssPlotEllipse(1:2,1:2) = B_ssPlot(n,1:2,1:2,i);
 
 Beta_ssPlotEllipse = zeros(2,2);
 
 Beta_ssPlotEllipse(1,1) = Beta_ss_saveINS_H(n,1,i)*ReConst;
  Beta_ssPlotEllipse(2,2) = Beta_ss_saveINS_H(n,2,i)*ReConst;
 
error_ellipse(B_ssPlotEllipse);
%plot( sqrt((Beta_ss_saveINS_H(n,1,i)*ReConst)^2+(Beta_ss_saveINS_H(n,2,i)*ReConst)^2),'r*');


error_ellipse(Beta_ssPlotEllipse);

end




%For ADM plots:





for i = startepoch:endepoch

for kk = 1:N_save(i)
    
for n = Excluded_PRNSub(kk,i)
    
    
     %P_save_subtempMODELPlot(1,1) = P_out_CSubTotal(n,7:8,7,i) ;  %note, in radians lat
     P_save_subtempMODELPlot = P_out_CSubTotal(n,26,26,i) ;  %note, in radians lon
          


B_ssPlot(n,1,i) = (P_save_subtempMODELPlot - P_save(26,26,i));



end
end
end













