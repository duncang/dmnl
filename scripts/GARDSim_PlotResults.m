
% plot results
% 
% figure();
% plot(t_save,UKF_x_hat_save(:,1)*180/pi);
% hold on;
% plot(t_save,(UKF_x_hat_save(:,1)+sqrt(UKF_P_save(:,1)))*180/pi,'r');
% plot(t_save,(UKF_x_hat_save(:,1)-sqrt(UKF_P_save(:,1)))*180/pi,'r');
% xlabel('Time (sec)');
% ylabel('Latitude Out');
% 
% 
% % plot results
% figure();
% plot(t_save,UKF_x_hat_save(:,2)*180/pi);
% hold on;
% plot(t_save,(UKF_x_hat_save(:,2)+sqrt(UKF_P_save(:,2)))*180/pi,'r');
% plot(t_save,(UKF_x_hat_save(:,2)-sqrt(UKF_P_save(:,2)))*180/pi,'r');
% xlabel('Time (sec)');
% ylabel('Longitude Out');
% 
% % plot results
% figure();
% plot(t_save,UKF_x_hat_save(:,3));
% hold on;
% plot(t_save,(UKF_x_hat_save(:,3)+sqrt(UKF_P_save(:,3))),'r');
% plot(t_save,(UKF_x_hat_save(:,3)-sqrt(UKF_P_save(:,3))),'r');
% xlabel('Time (sec)');
% ylabel('Height Out');

figure();
hold on;
plot(TimeINS,(Lat_error),'b');
grid on;
hold on;
plot(TimeINS,2*sqrt(UKF_P_save(:,1))*RMh,'g');
plot(TimeINS,-2*sqrt(UKF_P_save(:,1))*RMh,'g');
xlabel('Simulation Time (sec)');
ylabel('North Error (m)');
legend('North Error','2-\sigma Error Bound');
%axis([StartTime StopTime -5 5]);
hold off;

figure();
hold on;
plot(TimeINS,(Long_error),'b');
grid on;
hold on;
plot(TimeINS,2*sqrt(UKF_P_save(:,2))*RPh*cos(InitialPosition(1)),'g');
plot(TimeINS,-2*sqrt(UKF_P_save(:,2))*RPh*cos(InitialPosition(1)),'g');
xlabel('Simulation Time (sec)');
ylabel('East Error (m)');
legend('East Error','2-\sigma Error Bound');
%axis([TimeINS(1) TimeINS(NumberEpochsINS) -20 20]);
%axis([StartTime StopTime -5 5]);
hold off;

%% Horizontal error
North_var_m = (sqrt(UKF_P_save(:,1))*RMh).^2;
East_var_m = (sqrt(UKF_P_save(:,2))*RPh*cos( InitialPosition(1))).^2;
figure();
hold on; grid on;
plot(TimeINS, sqrt(Lat_error.^2+Long_error.^2));
plot(TimeINS,2*sqrt(North_var_m+East_var_m),'g');
xlabel('Simulation Time (sec)');
ylabel('Error (m)');
legend('True Error','2-\sigma Error Bound');
axis([StartTime StopTime 0 5]);
title('UKF Horizontal Error Bounds');
hold off;


%% vertical error
figure();
plot(TimeINS,abs(Height_error),'b');
hold on;
plot(TimeINS,2*sqrt(UKF_P_save(:,3)),'g');
plot(TimeINS,-2*sqrt(UKF_P_save(:,3)),'g');
xlabel('Simulation Time (sec)');
ylabel('Error (m)');
legend('True Error','2-\sigma Error Bound');
axis([StartTime StopTime 0 5]);
title('UKF Vertical Error Bounds');
hold off;
grid on;


% plot velocity
figure(); 
subplot(3,1,1),plot(TimeINS,Vel_error(:,1));hold on; grid on;ylabel('North (m/s)');title('Velocity Error');
subplot(3,1,1),plot(TimeINS,sqrt(UKF_P_save(:,4)),'g');
subplot(3,1,1),plot(TimeINS,-sqrt(UKF_P_save(:,4)),'g');
axis([StartTime StopTime -0.5 0.5]);
subplot(3,1,2),plot(TimeINS,Vel_error(:,2));hold on; grid on;ylabel('East (m/s)');
subplot(3,1,2),plot(TimeINS,sqrt(UKF_P_save(:,5)),'g');
subplot(3,1,2),plot(TimeINS,-sqrt(UKF_P_save(:,5)),'g');
axis([StartTime StopTime -0.5 0.5]);
subplot(3,1,3),plot(TimeINS,Vel_error(:,3));hold on; grid on;ylabel('Down (m/s)');
subplot(3,1,3),plot(TimeINS,sqrt(UKF_P_save(:,6)),'g');
subplot(3,1,3),plot(TimeINS,-sqrt(UKF_P_save(:,6)),'g');
axis([StartTime StopTime -0.5 0.5]);
xlabel('Simulation Time (seconds)');

% bias estimates
figure();
subplot(3,1,1),plot(TimeINS,(UKF_x_hat_save(:,14)-GyroBiasTruth(1,:)')*r2d,'r');grid on;hold on;ylabel('X Gyro Bias (deg/sec)');
subplot(3,1,1),plot(TimeINS,sqrt(UKF_P_save(:,14))*r2d,'g');
subplot(3,1,1),plot(TimeINS,-sqrt(UKF_P_save(:,14))*r2d,'g');
axis([StartTime StopTime -0.01 0.01]);
subplot(3,1,2),plot(TimeINS,(UKF_x_hat_save(:,15)-GyroBiasTruth(2,:)')*r2d,'r');grid on;hold on;ylabel('Y Gyro Bias (deg/sec)');
subplot(3,1,2),plot(TimeINS,sqrt(UKF_P_save(:,15))*r2d,'g');
subplot(3,1,2),plot(TimeINS,-sqrt(UKF_P_save(:,15))*r2d,'g');
axis([StartTime StopTime -0.01 0.01]);
subplot(3,1,3),plot(TimeINS,(UKF_x_hat_save(:,16)-GyroBiasTruth(3,:)')*r2d,'r');grid on;hold on;ylabel('Z Gyro Bias (deg/sec)');
subplot(3,1,3),plot(TimeINS,sqrt(UKF_P_save(:,16))*r2d,'g');
subplot(3,1,3),plot(TimeINS,-sqrt(UKF_P_save(:,16))*r2d,'g');
axis([StartTime StopTime -0.01 0.01]);
xlabel('Time (Sec)');


% accelerometer bias estimates
figure();
subplot(3,1,1),plot(TimeINS,UKF_x_hat_save(:,11)-AccelBiasTruth(1,:)','r');grid on;hold on;ylabel('X Acc Bias (m/s/s)');
subplot(3,1,1),plot(TimeINS,sqrt(UKF_P_save(:,11))*r2d,'g');
subplot(3,1,1),plot(TimeINS,-sqrt(UKF_P_save(:,11))*r2d,'g');
axis([StartTime StopTime -1 1]);
subplot(3,1,2),plot(TimeINS,UKF_x_hat_save(:,12)-AccelBiasTruth(2,:)','r');grid on;hold on;ylabel('Y Acc Bias (m/s/s)');
subplot(3,1,2),plot(TimeINS,sqrt(UKF_P_save(:,12))*r2d,'g');
subplot(3,1,2),plot(TimeINS,-sqrt(UKF_P_save(:,12))*r2d,'g');
axis([StartTime StopTime -1 1]);
subplot(3,1,3),plot(TimeINS,UKF_x_hat_save(:,13)-AccelBiasTruth(3,:)','r');grid on;hold on;ylabel('Z Acc Bias (m/s/s)');
subplot(3,1,3),plot(TimeINS,sqrt(UKF_P_save(:,13))*r2d,'g');
subplot(3,1,3),plot(TimeINS,-sqrt(UKF_P_save(:,13))*r2d,'g');
axis([StartTime StopTime -1 1]);
xlabel('Time (Sec)');

%% attitude errors
figure();
subplot(3,1,1),plot(TimeINS,phi_err*r2d,'r');grid on;hold on;ylabel('Roll (deg)');
%subplot(3,1,1),plot(TimeINS,2*sqrt(P_out_diag(7,:)*180/pi),'g');
%subplot(3,1,1),plot(TimeINS,-2*sqrt(P_out_diag(7,:)*180/pi),'g');
axis([StartTime StopTime -0.1 0.1]);
subplot(3,1,2),plot(TimeINS,theta_err*r2d,'r');grid on;hold on;ylabel('Pitch (deg)');
%subplot(3,1,2),plot(TimeINS,2*sqrt(P_out_diag(8,:)*180/pi),'g');
%subplot(3,1,2),plot(TimeINS,-2*sqrt(P_out_diag(8,:)*180/pi),'g');
axis([StartTime StopTime -0.1 0.1]);
subplot(3,1,3),plot(TimeINS,psi_err*r2d,'r');grid on;hold on;ylabel('Yaw (deg)');
%subplot(3,1,3),plot(TimeINS,2*sqrt(P_out_diag(9,:)*180/pi),'g');
%subplot(3,1,3),plot(TimeINS,-2*sqrt(P_out_diag(9,:)*180/pi),'g');
axis([StartTime StopTime -2 2]);
xlabel('Time (Sec)');

% plot attitude
% figure();
% plot(t_save,phi_q(:)*180/pi); hold on;
% plot(att_truth(1,:),att_truth(2,:)*180/pi,'r');
% xlabel('Time (sec)');
% ylabel('Roll Angle (deg)');
% figure();
% plot(t_save,theta_q(:)*180/pi); hold on;
% plot(att_truth(1,:),att_truth(3,:)*180/pi,'r');
% xlabel('Time (sec)');
% ylabel('Pitch Angle (deg)');
% figure();
% plot(t_save,psi_q(:)*180/pi); hold on;
% plot(att_truth(1,:),att_truth(4,:)*180/pi,'r');
% xlabel('Time (sec)');
% ylabel('Yaw Angle (deg)');

% figure();
% for i=1:85
%     plot((xs_save_lat(6002:7001,i)-xs_save_lat(6002:7001,1))*RM,'r');hold on;
% end
% xlabel('Time (sec)')
% ylabel('Latitude Error (metres)');
% grid on;

% 
% figure();
% plot(TimeGPS,Clock_error); hold on; grid on;
% plot(TimeINS,2*sqrt(UKF_P_save(:,17)),'g');
% plot(TimeINS,-2*sqrt(UKF_P_save(:,17)),'g');
% axis([StartTime StopTime -10 10]);
% xlabel('Simulation Time (sec)');
% ylabel('Error (metres)');
% title('Clock Estimation Error');


figure();
plot(TimeINS,UKF_x_hat_save(:,18)); hold on; grid on;
plot(TimeINS,2*sqrt(UKF_P_save(:,18)),'g');
plot(TimeINS,-2*sqrt(UKF_P_save(:,18)),'g');
axis([StartTime StopTime -1 1]);
xlabel('Simulation Time (sec)');
ylabel('Error (metres/sec)');
title('Frequency Estimation Error');


figure();hold on;grid on;
plot(TimeINS,sqrt(Lat_error.^2+Long_error.^2));
plot(TimeGPS,UKF_HPL*sqrt(RMh*RPh),'r','LineWidth',2);
axis([StartTime StopTime 0 50]);
legend('Horizontal Error','HPL');
xlabel('Simulation Time (sec)');
ylabel('Error (metres)');


figure();hold on;grid on;
plot(TimeINS,sqrt(Height_error.^2));
plot(TimeGPS,UKF_VPL,'r','LineWidth',2);
axis([StartTime StopTime 0 50]);
legend('Vertical Error','VPL');
xlabel('Simulation Time (sec)');
ylabel('Error (metres)');

figure();
hold on; grid on;
plot(TimeGPS,UKF_HPL*sqrt(RMh*RPh),'LineWidth',2);
plot(TimeGPS,UKF_VPL,'r','LineWidth',2);
legend('HPL','VPL');
%axis([StartTime StopTime 0 10]);
xlabel('Simulation Time (sec)');
ylabel('Protection Level (m)');
title('UKF GPS-INS Protection Level - RNAV Approach');


% plot lat error vs lsq
figure();
hold on; grid on;
plot(TimeINS,abs(Lat_error),'b');
plot(TimeGPS,abs(LSQ_LLH_error(:,1))*RM,'r');
legend('UKF INS','LSQ');
axis([StartTime StopTime 0 1]);
title('Latitude Error');

figure();
hold on; grid on;
plot(TimeINS,abs(Long_error),'b');
plot(TimeGPS,abs(LSQ_LLH_error(:,2))*RP*cos(-0.47),'r');
legend('UKF INS','LSQ');
axis([StartTime StopTime 0 1]);
title('Longitude Error');

figure();
hold on; grid on;
plot(TimeINS,Height_error,'b');
plot(TimeGPS,abs(LSQ_LLH_error(:,3)),'r');
legend('UKF INS','LSQ');
axis([StartTime StopTime 0 2]);
title('Height Error');


