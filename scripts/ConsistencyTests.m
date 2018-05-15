


%this is to test the sigma pr  values used in the least squares, by
%multiplying sigmah types dop..



for i = 4:120
    
    sigmah(i) = std(HPE_GPS_mc(:,i));
    
    sigmav(i) = std(VPE_GPS_mc(:,i));
    
  
    
    %actually i think i should multiply sigmah + sigma v times pdop?
    sigmaprH(i) =  sigmah(i)/DOP_Observed_save(i,3);
    
    sigmaprV(i) =  sigmav(i)/DOP_Observed_save(i,4);
    
    % sigmaprH(i) =  sigmah(i)*DOP_Observed_save(i,3);
    
    sigmapr(i) =  sqrt(sigmav(i)^2 + sigmah(i)^2)/DOP_Observed_save(i,1);
    
    
end

%yes sigma pr is about 2.5 metres on average, sometimes lower, sometimes
%higher. 
















%Note I didnt do consistency check for bias states, it makes the whole
%thing too large, even for the GPS-INS. So I just did the main states to
%check they were consistent. ie the PVA. Left clock out too hehehe. 



%NOTE, WITH ADM, if INV P is badly scaled, try using PINV, that turns it
%from being inconsistent to consistent, weird. 
%=============================
%consistency Tests for thesis results
%=============================

%NEES


%INS STATE 

DOF = 9;

%DOF = 17;

%DOF = 11;

XhighNEES = chi2inv(1- 0.05/2,DOF);
XlowNEES = chi2inv(0.05/2,DOF);


%convert roll pitch yaw errors to tilt errors


XhighCountNEES = 0;
XlowCountNEES = 0;
%plot normalised state errors. ie NEES Normalised  (state) Estimation Error Squared



for i = 5:70        
    
    
    %need to convert attitude to tilt error because thats what the
    %covariance is in 
    
    q0 = x_state_total_plusINS(1,i);
    q1 = x_state_total_plusINS(2,i);
    q2 = x_state_total_plusINS(3,i);
    q3 = x_state_total_plusINS(4,i);

    quat = [q0,q1,q2,q3];

    [euler] = QuatToEuler(quat);

    phitoconv = euler(1);
    thetatoconv = euler(2);
    psitoconv = euler(3);


phitoconv = Roll_truth1Hz(i)*pi/180;
    thetatoconv = Pitch_truth1Hz(i)*pi/180;
    psitoconv = Yaw_truth1Hz(i)*pi/180;


    T_alpha_TilttoEuler = [-cos(psitoconv)/cos(thetatoconv), -sin(psitoconv)/cos(thetatoconv), 0;
        sin(psitoconv), -cos(psitoconv), 0;
        -cos(psitoconv)*tan(thetatoconv), -sin(psitoconv)*tan(thetatoconv), -1;];


    
    Tilt = T_alpha_TilttoEuler'*[ RollErrorINS_C(i)*pi/180,  PitchErrorINS_C(i)*pi/180,  YawErrorINS_C(i)*pi/180]';  %these were in degrees, need to convert back to radians


tiltx(i) = Tilt(1);
tilty(i) = Tilt(2);
tiltz(i) = Tilt(3);    

    

%for the 9 states

if DOF == 9
  ErrorVector = [tiltx(i); tilty(i); tiltz(i); N_VelErrorsINS_C(i);E_VelErrorsINS_C(i);D_VelErrorsINS_C(i); Lat_PosErrorINS_C(i);Lon_PosErrorINS_C(i);Hgt_PosErrorINS_C(i)];
     
  NEES_C(i) = ErrorVector'*inv(P_save(1:9,1:9,i))*ErrorVector;
end


%for 11 (with clock)
if DOF == 11
  ErrorVector = [tiltx(i); tilty(i); tiltz(i); N_VelErrorsINS_C(i);E_VelErrorsINS_C(i);D_VelErrorsINS_C(i); Lat_PosErrorINS_C(i);Lon_PosErrorINS_C(i);Hgt_PosErrorINS_C(i);clkbiaserror_C(i); clkdrifterror_C(i)];
     
  NEES_C(i) = ErrorVector'*inv(P_save(1:11,1:11,i))*ErrorVector;
end



%for the 17 states

if DOF == 17
ErrorVector = [tiltx(i); tilty(i); tiltz(i); N_VelErrorsINS_C(i);E_VelErrorsINS_C(i);D_VelErrorsINS_C(i); Lat_PosErrorINS_C(i);Lon_PosErrorINS_C(i);Hgt_PosErrorINS_C(i); clkbiaserror_C(i); clkdrifterror_C(i);BiasErrorXGyro_INS_C(i)*pi/180;BiasErrorYGyro_INS_C(i)*pi/180;BiasErrorZGyro_INS_C(i)*pi/180;BiasErrorXAccel_INS_C(i);BiasErrorYAccel_INS_C(i);BiasErrorZAccel_INS_C(i)];
    
NEES_C(i) = ErrorVector'*inv(P_save(1:17,1:17,i))*ErrorVector;
end


% 
% BiasErrorXGyro_INS_C = (XGyroBiasTruth_INS_C - x_state_total_INSBIAS(1,startepochplot:endepochplot))*180/pi;
% BiasErrorYGyro_INS_C = (YGyroBiasTruth_INS_C - x_state_total_INSBIAS(2,startepochplot:endepochplot))*180/pi;
% BiasErrorZGyro_INS_C = (ZGyroBiasTruth_INS_C - x_state_total_INSBIAS(3,startepochplot:endepochplot))*180/pi;
% BiasErrorXAccel_INS_C = (XAccelBiasTruth_INS_C - x_state_total_INSBIAS(4,startepochplot:endepochplot));
% BiasErrorYAccel_INS_C = (YAccelBiasTruth_INS_C - x_state_total_INSBIAS(5,startepochplot:endepochplot));
% BiasErrorZAccel_INS_C = (ZAccelBiasTruth_INS_C - x_state_total_INSBIAS(6,startepochplot:endepochplot));
% 
% 
% 
% 
% 
%   clkbiaserror_C(i) 
% clkdrifterror_C(i) 




%determine number of points outside bounds
if NEES_C(i) > XhighNEES
    XhighCountNEES = XhighCountNEES+1;
end

if NEES_C(i) < XlowNEES
    XlowCountNEES = XlowCountNEES+1;
end


end




%this gives the upper bound 95%
close all;

plot(NEES_C,'k-','LineWidth',2); title 'Normalized Estimation Error Squared' ; xlabel('Time (s)');  ylabel('NEES');
hold;
plot(ones(1,length(NEES_C))*XhighNEES,'k--');
plot(ones(1,length(NEES_C))*XlowNEES,'k--');







%innovations for INS

XhighCountNIS = 0;
XlowCountNIS = 0;
for i = 5:70    
    
    clear v_innovation
    
    nz(i) = N_save(i)*2 ;
    v_innovation_INS =  z_save(i,1:nz(i))';
    %v_innovation =  z_save(i,1:10)';
     
    NIS_INS(i) = v_innovation_INS'*inv(V_C_save(1:nz(i),1:nz(i),i))*v_innovation_INS;
    
    %because R is singular might cause problms with invV
    %NIS(i) = v_innovation'*inv(V_C_save(1:10,1:10,i))*v_innovation;
    
% XhighNIS = chi2inv(0.95,nz(i));
% XlowNIS = chi2inv(0.05,nz(i));


XhighNIS = chi2inv(1- 0.05/2,nz(i));
XlowNIS = chi2inv(0.05/2,nz(i));



    
%determine number of points outside bounds
if NIS_INS(i) > XhighNIS
    XhighCountNIS = XhighCountNIS+1;
end

if NIS_INS(i) < XlowNIS
    XlowCountNIS = XlowCountNIS+1;
end    
    

end



%Xhigh = chi2inv(0.95,10)
%Xlow = chi2inv(0.05,10)

plot(NIS_INS,'k-','LineWidth',2); title 'Normalized Innovations Squared' ; xlabel('Time (s)');  ylabel('NIS');
hold;
plot(ones(1,length(NIS_INS))*XhighNIS,'k--');
plot(ones(1,length(NIS_INS))*XlowNIS,'k--');






%===============================================
%do covariance matching to see if HPL meets Pmd

%===============================================

%Pr that pr md = Pmd = Pr that HPE > HPL and pr that test stat < detection
%threshold. 

starti = 5;
endi = 100;

for i = starti:endi        
   
    
    
    %for IMU
    
    
    lambda_k = max(lambda_ss_outINS_H(i,:)); %noncentrality parameter is max of test statistics.
    %threshold squared
   TD_squared = aTDINS_H(i);
    N = N_save(i);
    DOF = 2;  %for HPL
  sigma_re = sqrt((sqrt(P_save(7,7,i)))^2 + (sqrt(P_save(8,8,i)))^2)*ReConst;
    PL = HPL_INS_C(i);
    de = HPE_INS_C(i);        
    
    [Pmd_Est_INS_H(i), PrHPE_INS_H(i), PrTestStat_INS_H(i)] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de);
       
    
    
    
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

    
    
       
    


% %using noncentral chisquare distribution
% %lambda_k = sqrt(lambda_k);
% 
% 
% 
% %PrTestStat(i) = ncx2cdf(TD_squared,N_save(i),lambda_k);
% 
% PrTestStat(i) = ncx2pdf(TD_squared,2,lambda_k);
% 
% 
% PrTestStat(i) = 1 - PrTestStat(i)*N_save(i); 
% 
% 
% %PrTestStat(i) = PrTestStat(i);%*N_save(i);     %need to times by N here because threshold is calculated with pfd, or is it divide? 
% 
% %PrTestStat(i) = ncx2cdf(TD_squared,2,lambda_k);   %DOF should be 2 i think?
% 
% % 
% % NCX2CDF Noncentral chi-square cumulative distribution function (cdf).
% %      Returns the noncentral chi-square cdf with V 
% %     degrees of freedom and noncentrality parameter, DELTA, at the values 
% %     in X.
% 
% 
% %derive prob that HPE is > HPL
% 
% 
% %sigma_re  =  std(HPE_INS_C(5:100));  %one sigma value of horizontal error
% 
% 
% 
% %sigma_re  =  HPE_INS_C(i);  %one sigma value of horizontal error
% 
% 
% sigma_re = sqrt((sqrt(P_save(7,7,i)))^2 + (sqrt(P_save(8,8,i)))^2);
% 
% sigma_re = sigma_re*ReConst; %convert to metres
% 
% 
% %radial bias error due to bias alone
% %de = 5;  %this is pbias or r bias? how to calculate?
% 
% 
% de = pbias_norm_H;
% 
% de = HPE_Bj_max_save(i); %i think its this one, note this isnt quite correct, will need to get the value at each time..
% 
% de = HPE_INS_C(i); 
% 
% PrHPE(i) = normcdf(HPL_INS_C(i)-de,0,sigma_re);
% 
% PrHPE(i) = 1-PrHPE(i) ;  %do i need this? Because matlab goes from -inf to x, i want tfrom x to inf..i want the tail region on the main region
% 
% %PrTestStat(i) = 1-PrTestStat(i);
% 
% 
% Pmd_Est_H(i) = PrHPE(i)*PrTestStat(i); 
% 
% 
% end
% 
% 
% 






for i = starti:endi    
    
%     clear v_innovation_INS
%     
%     nz(i) = N_save(i)*2 ;
%     v_innovation_INS =  z_save(i,1:nz(i))';    
    
    
%1. derive prob of no alarm, ie that test stat < threshold


%get noncentrality parameter


%get innovation solely in presence of ramp error

%dr = v_innovation_INS;

%lambda_k = dr'*inv(V_C_save(1:nz(i),1:nz(i),i))*dr;   %noncentrality parameter

lambda_k = max(lambda_ss_outGPS_H(i,:)); %noncentrality parameter is max of test statistics.

%using noncentral chisquare distribution
%lambda_k = sqrt(lambda_k);

%threshold squared
TD_squared = aTDGPS_H(i);


%PrTestStat(i) = ncx2cdf(TD_squared,N_save(i),lambda_k);

PrTestStat(i) = ncx2pdf(TD_squared,2,lambda_k);


PrTestStat(i) = 1 - PrTestStat(i)*N_save(i); 


%PrTestStat(i) = PrTestStat(i);%*N_save(i);     %need to times by N here because threshold is calculated with pfd, or is it divide? 

%PrTestStat(i) = ncx2cdf(TD_squared,2,lambda_k);   %DOF should be 2 i think?

% 
% NCX2CDF Noncentral chi-square cumulative distribution function (cdf).
%      Returns the noncentral chi-square cdf with V 
%     degrees of freedom and noncentrality parameter, DELTA, at the values 
%     in X.


%derive prob that HPE is > HPL


%sigma_re  =  std(HPE_INS_C(5:100));  %one sigma value of horizontal error



sigma_re  =  HPE_GPS_C(i);  %one sigma value of horizontal error


%sigma_re = sqrt((sqrt(P_save(7,7,i)))^2 + (sqrt(P_save(8,8,i)))^2);



sigma_re = sqrt((VarSolutionVec_Subsave(i,1)*2.5^2)^2 + (VarSolutionVec_Subsave(i,2)*2.5^2)^2) 



%radial bias error due to bias alone
%de = 5;  %this is pbias or r bias? how to calculate?


de = pbias_norm_H;

de = HPE_Bj_max_save(i); %i think its this one, note this isnt quite correct, will need to get the value at each time..

de = 3.64922;

PrHPE(i) = normcdf(HPL_GPS(i)-de,0,sigma_re);

PrHPE(i) = 1-PrHPE(i) ;  %do i need this? Because matlab goes from -inf to x, i want tfrom x to inf..i want the tail region on the main region

%PrTestStat(i) = 1-PrTestStat(i);


Pmd_Est_H(i) = PrHPE(i)*PrTestStat(i); 


end







%For vertical

lambda_k = max(lambda_ss_outINS_V(i,:)); %noncentrality parameter is max of test statistics.

%using noncentral chisquare distribution
%lambda_k = sqrt(lambda_k);

%threshold squared
TD_squared = aTDINS_V(i);


PrTestStat_V(i) = ncx2cdf(TD_squared,N_save(i),lambda_k);




PrTestStat_V(i) = 1 - PrTestStat_V(i); 


PrTestStat_V(i) = PrTestStat_V(i);%*N_save(i);  %need to times by N here becuase with NSS method pfd is divided by N for each subfilter.  or should this be divide by N?





sigma_re  =  VPE_INS_C(i);  %one sigma value of horizontal error



de = pbias_norm_V;


de = VPE_Bj_max_save(i); %i think its this one, note this isnt quite correct, will need to get the value at each time..


PrVPE(i) = normcdf(VPL_INS_C(i)-de,0,sigma_re);

PrVPE(i) = 1-PrVPE(i) ;  %do i need this? Because matlab goes from -inf to x, i want tfrom x to inf..i want the tail region on the main region


Pmd_Est_V(i) = PrVPE(i)*PrTestStat_V(i); 




%FOR GPS




lambda_k = max(lambda_ss_outGPS_H(i,:)); %noncentrality parameter is max of test statistics.

%using noncentral chisquare distribution
%lambda_k = sqrt(lambda_k);

%threshold squared
TD_squared = aTDGPS_H(i);





PrTestStat_GPS(i) = ncx2cdf(TD_squared,N_save(i),lambda_k);


%PrTestStat(i) = ncx2pdf(TD_squared,N_save(i),lambda_k);

PrTestStat_GPS(i) = 1 - PrTestStat_GPS(i); 


PrTestStat_GPS(i) = PrTestStat_GPS(i);%*N_save(i);

%PrTestStat(i) = ncx2cdf(TD_squared,2,lambda_k);   %DOF should be 2 i think?

% 
% NCX2CDF Noncentral chi-square cumulative distribution function (cdf).
%      Returns the noncentral chi-square cdf with V 
%     degrees of freedom and noncentrality parameter, DELTA, at the values 
%     in X.




%derive prob that HPE is > HPL


%sigma_re  =  std(HPE_INS_C(5:100));  %one sigma value of horizontal error



sigma_re  =  HPE_GPS_C(i);  %one sigma value of horizontal error


%sigma_re = sqrt((sqrt(P_save(7,7,i)))^2 + (sqrt(P_save(8,8,i)))^2)

%sigma_re = sigma_re*ReConst; %convert to metres


%radial bias error due to bias alone
%de = 5;  %this is pbias or r bias? how to calculate?


de = pbias_norm_H;


de = HPE_Bj; %i think its this one, note this isnt quite correct, will need to get the value at each time..


PrHPE_GPS(i) = normcdf(HPL_GPS(i)-de,0,sigma_re);

PrHPE_GPS(i) = 1-PrHPE_GPS(i) ;  %do i need this? Because matlab goes from -inf to x, i want tfrom x to inf..i want the tail region on the main region

%PrTestStat(i) = 1-PrTestStat(i);


Pmd_Est_H_GPS(i) = PrHPE_GPS(i)*PrTestStat_GPS(i); 





%for GPS vertical




lambda_k = max(lambda_ss_outGPS_V(i,:)); %noncentrality parameter is max of test statistics.

%using noncentral chisquare distribution
%lambda_k = sqrt(lambda_k);

%threshold squared
TD_squared = aTDGPS_V(i);





PrTestStat_GPS_V(i) = ncx2cdf(TD_squared,N_save(i),lambda_k);


%PrTestStat(i) = ncx2pdf(TD_squared,N_save(i),lambda_k);

PrTestStat_GPS_V(i) = 1 - PrTestStat_GPS_V(i); 


PrTestStat_GPS_V(i) = PrTestStat_GPS_V(i);%*N_save(i);

%PrTestStat(i) = ncx2cdf(TD_squared,2,lambda_k);   %DOF should be 2 i think?

% 
% NCX2CDF Noncentral chi-square cumulative distribution function (cdf).
%      Returns the noncentral chi-square cdf with V 
%     degrees of freedom and noncentrality parameter, DELTA, at the values 
%     in X.




%derive prob that HPE is > HPL


%sigma_re  =  std(HPE_INS_C(5:100));  %one sigma value of horizontal error



sigma_re  =  VPE_GPS_C(i);  %one sigma value of horizontal error


%sigma_re = sqrt((sqrt(P_save(7,7,i)))^2 + (sqrt(P_save(8,8,i)))^2)

%sigma_re = sigma_re*ReConst; %convert to metres


%radial bias error due to bias alone
%de = 5;  %this is pbias or r bias? how to calculate?


de = pbias_norm_V;


de = VPE_Bj; %i think its this one, note this isnt quite correct, will need to get the value at each time..


PrVPE_GPS(i) = normcdf(VPL_GPS(i)-de,0,sigma_re);

PrVPE_GPS(i) = 1-PrVPE_GPS(i) ;  %do i need this? Because matlab goes from -inf to x, i want tfrom x to inf..i want the tail region on the main region

%PrTestStat(i) = 1-PrTestStat(i);


Pmd_Est_V_GPS(i) = PrVPE_GPS(i)*PrTestStat_GPS_V(i); 













end



%do the test for HPL_H0 as well. ??



 P_H0_H = 1e-10;
 P_H0_V = 1.9e-9;
        
        
  %this if for H1, though i dont have the H1 subscript in here
  
 Pmd_H = 0.0025;
 Pmd_V = 0.045;










%MODEL STATE 

DOF = 6;

%DOF = 11;  %can't do NEES with the biases, because don't know the true value of the biases

XhighNEES = chi2inv(1- 0.05/2,DOF);
XlowNEES = chi2inv(0.05/2,DOF);


%convert roll pitch yaw errors to tilt errors


XhighCountNEES = 0;
XlowCountNEES = 0;
%plot normalised state errors. ie NEES Normalised  (state) Estimation Error Squared



for i = 5:70        
    
    
    q0 = x_state_total_plusMODEL(1,i);
    q1 = x_state_total_plusMODEL(2,i);
    q2 = x_state_total_plusMODEL(3,i);
    q3 = x_state_total_plusMODEL(4,i);

    quat = [q0,q1,q2,q3];

    [euler] = QuatToEuler(quat);

    phitoconv = euler(1);
    thetatoconv = euler(2);
    psitoconv = euler(3);


phitoconv = Roll_truth1Hz(i)*pi/180;
    thetatoconv = Pitch_truth1Hz(i)*pi/180;
    psitoconv = Yaw_truth1Hz(i)*pi/180;


    T_alpha_TilttoEuler = [-cos(psitoconv)/cos(thetatoconv), -sin(psitoconv)/cos(thetatoconv), 0;
        sin(psitoconv), -cos(psitoconv), 0;
        -cos(psitoconv)*tan(thetatoconv), -sin(psitoconv)*tan(thetatoconv), -1;];


    
    Tilt = T_alpha_TilttoEuler'*[ RollErrorMODEL_C(i)*pi/180,  PitchErrorMODEL_C(i)*pi/180,  YawErrorMODEL_C(i)*pi/180]';


tiltxMODEL(i) = Tilt(1);
tiltyMODEL(i) = Tilt(2);
tiltzMODEL(i) = Tilt(3);    
    
    



if DOF == 9

ErrorVector = [tiltxMODEL(i); tiltyMODEL(i); tiltzMODEL(i); N_VelErrorsMODEL_C(i);E_VelErrorsMODEL_C(i);D_VelErrorsMODEL_C(i); Lat_PosErrorMODEL_C(i);Lon_PosErrorMODEL_C(i);Hgt_PosErrorMODEL_C(i)];
    
NEES_CMODEL(i) = ErrorVector'*inv(P_save(18:26,18:26,i))*ErrorVector;
end



%because I don't know what the true bias is for the ADM, just do NEES for
%the 11 states

if DOF == 11
  
ErrorVector = [tiltxMODEL(i); tiltyMODEL(i); tiltzMODEL(i); N_VelErrorsMODEL_C(i);E_VelErrorsMODEL_C(i);D_VelErrorsMODEL_C(i); Lat_PosErrorMODEL_C(i);Lon_PosErrorMODEL_C(i);Hgt_PosErrorMODEL_C(i);clkbiaserrorMODEL_C(i); clkdrifterrorMODEL_C(i)];
   
NEES_CMODEL(i) = ErrorVector'*inv(P_save(18:28,18:28,i))*ErrorVector;
end




%this is for the GPS KF

if DOF == 6
  
ErrorVector = [ N_VelErrorsMODEL_C(i);E_VelErrorsMODEL_C(i);D_VelErrorsMODEL_C(i); Lat_PosErrorMODEL_C(i);Lon_PosErrorMODEL_C(i);Hgt_PosErrorMODEL_C(i)];
   
NEES_CMODEL(i) = ErrorVector'*inv(P_save(21:26,21:26,i))*ErrorVector;
end






%determine number of points outside bounds
if NEES_CMODEL(i) > XhighNEES
    XhighCountNEES = XhighCountNEES+1;
end

if NEES_CMODEL(i) < XlowNEES
    XlowCountNEES = XlowCountNEES+1;
end


end



%this gives the upper bound 95%
close all;

plot(NEES_CMODEL,'k-','LineWidth',2); title 'Normalized Estimation Error Squared ADM' ; xlabel('Time (s)');  ylabel('NEES');
hold;
plot(ones(1,length(NEES_CMODEL))*XhighNEES,'k--');
plot(ones(1,length(NEES_CMODEL))*XlowNEES,'k--');











%===========================
%test innovations
%===========================

%p.424 of bar shalom
%Monitor innovations 


 clear v_innovation   
 
 
for i = 5:70        
 
       
    %nz(i) = N_save(i)*4 + NumberCriticalStates;
    
    
    %v_innovation =  z_save(i,1:nz(i))';
    
        nz(i) = N_save(i)*2;
    
    
    if Fused == 0
               
        
    
    
     v_innovation =  z_save(i,2*N_save(i)+1:N_save(i)*4 )';
    %v_innovation =  z_save(i,1:10)';
     
    %NIS_MODEL(i) =
    %v_innovation'*pinv(V_C_save(1:nz(i),1:nz(i),i))*v_innovation;
    
     NIS_MODEL(i) = v_innovation'*inv(V_C_save(2*N_save(i)+1:N_save(i)*4 ,2*N_save(i)+1:N_save(i)*4 ,i))*v_innovation;
    
    end
    
    
     if Fused == 1
         
     v_innovation =  z_save(i,2*N_save(i)+3+1:N_save(i)*4 + NumberCriticalStates)';
    
     NIS_MODEL(i) = v_innovation'*pinv(V_C_save(2*N_save(i)+3+1:N_save(i)*4 + NumberCriticalStates,2*N_save(i)+3+1:N_save(i)*4 + NumberCriticalStates,i))*v_innovation;
    
    end
    
    
    
    %because R is singular might cause problms with invV
    %NIS(i) = v_innovation'*inv(V_C_save(1:10,1:10,i))*v_innovation;
    
%    
% XhighNIS = chi2inv(0.95,nz(i));
% XlowNIS = chi2inv(0.05,nz(i));


%moving average of the normalized innovatiosn squared over a sliding window
%of s sample times

% evs(i) = 0; 
% 
% if i > 8
% s = 5; 
% for scount = 1:s
% evs(i) = evs(i) + NIS_MODEL(i+1-scount);
% 
% end
% 
% else
%  s = 1;
%  evs(i)  =  NIS_MODEL(i); 
% end
%     


end %for i = 5:100        



XhighNIS = chi2inv(1- 0.05/2,nz(i));
XlowNIS = chi2inv(0.05/2,nz(i));




%if it's outside the bounds then try to force it down by changing Q

%if evs(i) > (XhighNIS/2)+XlowNIS
%if NIS_MODEL(i) > XhighNIS


%end

%Xhigh = chi2inv(0.95,10)
%Xlow = chi2inv(0.05,10)

plot(NIS_MODEL,'k-','LineWidth',2); title 'Normalized Innovations Squared' ; xlabel('Time (s)');  ylabel('NIS');
hold;
plot(ones(1,length(NIS_MODEL))*XhighNIS,'k--');
plot(ones(1,length(NIS_MODEL))*XlowNIS,'k--');










%save copy for comparing fused with unfused consistency
if Fused == 0
    
NEES_C_NF = NEES_C;  %for INS
NEES_CMODEL_NF = NEES_CMODEL;



NIS_INS_NF = NIS_INS;
NIS_MODEL_NF = NIS_MODEL;

end








%For Combined Filter ie consistency of the whole thing..



%put NEES code in here..




XhighCountNIS = 0;
XlowCountNIS = 0;
for i = 5:100    
    
    clear v_innovation
    
    nz(i) = N_save(i)*4 + NumberCriticalStates;
    v_innovation =  z_save(i,1:nz(i))';
    %v_innovation =  z_save(i,1:10)';
     
    NIS_C(i) = v_innovation'*pinv(V_C_save(1:nz(i),1:nz(i),i))*v_innovation;
    
    %because R is singular might cause problms with invV
    %NIS(i) = v_innovation'*inv(V_C_save(1:10,1:10,i))*v_innovation;
    
XhighNIS = chi2inv(0.95,nz(i));
XlowNIS = chi2inv(0.05,nz(i));
    
%determine number of points outside bounds
if NIS_C(i) > XhighNIS
    XhighCountNIS = XhighCountNIS+1;
end

if NIS_C(i) < XlowNIS
    XlowCountNIS = XlowCountNIS+1;
end    
    

end

















%Xhigh = chi2inv(0.95,10)
%Xlow = chi2inv(0.05,10)

plot(NIS_C,'k-','LineWidth',2); title 'Normalized Innovations Squared' ; xlabel('Time (s)');  ylabel('NIS');
hold;
plot(ones(1,length(NIS_C))*XhighNIS,'k--');
plot(ones(1,length(NIS_C))*XlowNIS,'k--');



%V_C_save((1:N_save(i)*4+NumberCriticalStates),(1:N_save(i)*4+NumberCriticalStates),i)



%if biased then undertake bias test.















%end consistency tests for thesis..
%===========================================



%whiteness test

%z_save is really the innovation (measured - predicted observations)


for i = beginepoch:endepoch

nz = N_save(i)*4 + NumberCriticalStates;
v_innovation =  z_save(i,1:nz);
 
 
end
%degrees of freedom in innovations (this will vary according to number of
%measurements

 
 
 %number of monte carlo runs
 Nmc = 1;
 
 pk_j = (1/sqrt(nz))*v_innovation'*((v_innovation*v_innovation')^-0.5)*((v_innovation*v_innovation')^-0.5)*v_innovation;

 
  num = 1;
 for k = 1:nz
     
     
     for j = 1:nz
         
         pk_j(num) = (1/sqrt(nz))*v_innovation(k)'*((v_innovation(k)*v_innovation(k)')^-0.5)*((v_innovation(j)*v_innovation(j)')^-0.5)*v_innovation(j);
         num = num+1
     end
     
 end 
 
 
 for L = 1:10
 for k = beginepoch:endepoch
     v_innovationk =  z_save(k,L);
     
     for j = beginepoch:endepoch
         
      v_innovationj =  z_save(j,L);
     pk_j(L,:) = v_innovationk*v_innovationj*((v_innovationk^2)*(v_innovationj^2))^-0.5;
 
     end
 end
 end 
 

  for k = 1:nz     
     
     for j = 1:nz
         
         
         pk_j(k,j) = v_innovation(k)*v_innovation(j)*((v_innovation(k)^2)*(v_innovation(j)^2))^-0.5;
         
     end
     
  end
  
 
 
%to do a 50 run monte carlo average




%plot NEES

%plot NIS (normalized innovation squared)


%plot innovation autocorrelation see p 241 of Bar Shalom

total = sqrt((t212minust112_save(:,1)*Rn).^2+ (t212minust112_save(:,2)*Rn).^2 + (t212minust112_save(:,3)).^2)

plot(total);







%=========================================================================
%=========================================================================
%=========================================================================

%CONSISTENCY TESTS For GPSINSMODEL_ONLY_EKF100HzNewJulierPQRstates
%this one includes testing the p q r states for consistency

%=========================================================================
%=========================================================================
%=========================================================================


%MODEL STATE 

%DOF = 9;


%without clock:
DOF = 12;  %can't do NEES with the biases, because don't know the true value of the biases




%DOF = 14;  %can't do NEES with the biases, because don't know the true value of the biases





XhighNEES = chi2inv(1- 0.05/2,DOF);
XlowNEES = chi2inv(0.05/2,DOF);


%convert roll pitch yaw errors to tilt errors


XhighCountNEES = 0;
XlowCountNEES = 0;
%plot normalised state errors. ie NEES Normalised  (state) Estimation Error Squared



for i = 5:100        
    
    
    q0 = x_state_total_plusMODEL(1,i);
    q1 = x_state_total_plusMODEL(2,i);
    q2 = x_state_total_plusMODEL(3,i);
    q3 = x_state_total_plusMODEL(4,i);

    quat = [q0,q1,q2,q3];

    [euler] = QuatToEuler(quat);

    phitoconv = euler(1);
    thetatoconv = euler(2);
    psitoconv = euler(3);


phitoconv = Roll_truth1Hz(i)*pi/180;
    thetatoconv = Pitch_truth1Hz(i)*pi/180;
    psitoconv = Yaw_truth1Hz(i)*pi/180;


    T_alpha_TilttoEuler = [-cos(psitoconv)/cos(thetatoconv), -sin(psitoconv)/cos(thetatoconv), 0;
        sin(psitoconv), -cos(psitoconv), 0;
        -cos(psitoconv)*tan(thetatoconv), -sin(psitoconv)*tan(thetatoconv), -1;];


    
    Tilt = T_alpha_TilttoEuler'*[ RollErrorMODEL_C(i)*pi/180,  PitchErrorMODEL_C(i)*pi/180,  YawErrorMODEL_C(i)*pi/180]';


tiltxMODEL(i) = Tilt(1);
tiltyMODEL(i) = Tilt(2);
tiltzMODEL(i) = Tilt(3);    
    
    



if DOF == 9

ErrorVector = [tiltxMODEL(i); tiltyMODEL(i); tiltzMODEL(i); N_VelErrorsMODEL_C(i);E_VelErrorsMODEL_C(i);D_VelErrorsMODEL_C(i); Lat_PosErrorMODEL_C(i);Lon_PosErrorMODEL_C(i);Hgt_PosErrorMODEL_C(i)];
    
NEES_CMODEL(i) = ErrorVector'*inv(P_save(18:26,18:26,i))*ErrorVector;
end



%because I don't know what the true bias is for the ADM, just do NEES for
%the 11 states


if DOF == 12
  
ErrorVector = [tiltxMODEL(i); tiltyMODEL(i); tiltzMODEL(i); N_VelErrorsMODEL_C(i);E_VelErrorsMODEL_C(i);D_VelErrorsMODEL_C(i); Lat_PosErrorMODEL_C(i);Lon_PosErrorMODEL_C(i);Hgt_PosErrorMODEL_C(i);p_error_model(i);q_error_model(i);r_error_model(i);  ];
   
NEES_CMODEL(i) = ErrorVector'*inv(P_save(18:29,18:29,i))*ErrorVector;
end





if DOF == 14
  
ErrorVector = [tiltxMODEL(i); tiltyMODEL(i); tiltzMODEL(i); N_VelErrorsMODEL_C(i);E_VelErrorsMODEL_C(i);D_VelErrorsMODEL_C(i); Lat_PosErrorMODEL_C(i);Lon_PosErrorMODEL_C(i);Hgt_PosErrorMODEL_C(i);p_error_model(i);q_error_model(i);r_error_model(i);clkbiaserrorMODEL_C(i); clkdrifterrorMODEL_C(i); ];
   
NEES_CMODEL(i) = ErrorVector'*inv(P_save(18:31,18:31,i))*ErrorVector;
end





%determine number of points outside bounds
if NEES_CMODEL(i) > XhighNEES
    XhighCountNEES = XhighCountNEES+1;
end

if NEES_CMODEL(i) < XlowNEES
    XlowCountNEES = XlowCountNEES+1;
end


end



%this gives the upper bound 95%
close all;

plot(NEES_CMODEL,'k-','LineWidth',2); title 'Normalized Estimation Error Squared ADM' ; xlabel('Time (s)');  ylabel('NEES');
hold;
plot(ones(1,length(NEES_CMODEL))*XhighNEES,'k--');
plot(ones(1,length(NEES_CMODEL))*XlowNEES,'k--');



%save copy for comparing fused with unfused consistency
if Fused == 0
    
NEES_C_NF = NEES_C;  %for INS
NEES_CMODEL_NF = NEES_CMODEL;


end








%innovations for INS

XhighCountNIS = 0;
XlowCountNIS = 0;
for i = 5:100    
    
    clear v_innovation
    
    nz(i) = N_save(i)*2 ;
    v_innovation_INS =  z_save(i,1:nz(i))';
    %v_innovation =  z_save(i,1:10)';
     
    NIS_INS(i) = v_innovation_INS'*inv(V_C_save(1:nz(i),1:nz(i),i))*v_innovation_INS;
    
    %because R is singular might cause problms with invV
    %NIS(i) = v_innovation'*inv(V_C_save(1:10,1:10,i))*v_innovation;
    
% XhighNIS = chi2inv(0.95,nz(i));
% XlowNIS = chi2inv(0.05,nz(i));


XhighNIS = chi2inv(1- 0.05/2,nz(i));
XlowNIS = chi2inv(0.05/2,nz(i));



    
%determine number of points outside bounds
if NIS_INS(i) > XhighNIS
    XhighCountNIS = XhighCountNIS+1;
end

if NIS_INS(i) < XlowNIS
    XlowCountNIS = XlowCountNIS+1;
end    
    

end



%Xhigh = chi2inv(0.95,10)
%Xlow = chi2inv(0.05,10)

plot(NIS_INS,'k-','LineWidth',2); title 'Normalized Innovations Squared' ; xlabel('Time (s)');  ylabel('NIS');
hold;
plot(ones(1,length(NIS_INS))*XhighNIS,'k--');
plot(ones(1,length(NIS_INS))*XlowNIS,'k--');









