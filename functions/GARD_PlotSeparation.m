%% plot separation

% lat-long
figure(); hold on; grid on;
plot(UKF_x_hat_save(:,2),UKF_x_hat_save(:,1),'r')
plot(UKF_x_hat_save_Sub(:,2,1),UKF_x_hat_save_Sub(:,1,1),'b')
plot(UKF_x_hat_save_Sub(:,2,2),UKF_x_hat_save_Sub(:,1,2),'g')
plot(UKF_x_hat_save_Sub(:,2,3),UKF_x_hat_save_Sub(:,1,3),'k')
plot(UKF_x_hat_save_Sub(:,2,4),UKF_x_hat_save_Sub(:,1,4),'m')
plot(UKF_x_hat_save_Sub(:,2,5),UKF_x_hat_save_Sub(:,1,5),'y')
plot(UKF_x_hat_save_Sub(:,2,6),UKF_x_hat_save_Sub(:,1,6),'c')
plot(UKF_x_hat_save_Sub(:,2,7),UKF_x_hat_save_Sub(:,1,7),'b--')
plot(UKF_x_hat_save_Sub(:,2,8),UKF_x_hat_save_Sub(:,1,8),'g--')

% latitude
figure(); hold on; grid on;
plot(TimeINS,UKF_x_hat_save(:,1),'r')
plot(TimeGPS,UKF_x_hat_save_Sub(:,1,1))
plot(TimeGPS,UKF_x_hat_save_Sub(:,1,2))
plot(TimeGPS,UKF_x_hat_save_Sub(:,1,3))
plot(TimeGPS,UKF_x_hat_save_Sub(:,1,4))
plot(TimeGPS,UKF_x_hat_save_Sub(:,1,5))
plot(TimeGPS,UKF_x_hat_save_Sub(:,1,6))
plot(TimeGPS,UKF_x_hat_save_Sub(:,1,7))
plot(TimeGPS,UKF_x_hat_save_Sub(:,1,8))

% longitude
figure(); hold on; grid on;
plot(TimeINS,UKF_x_hat_save(:,2),'r')
plot(TimeGPS,UKF_x_hat_save_Sub(:,2,1))
plot(TimeGPS,UKF_x_hat_save_Sub(:,2,2))
plot(TimeGPS,UKF_x_hat_save_Sub(:,2,3))
plot(TimeGPS,UKF_x_hat_save_Sub(:,2,4))
plot(TimeGPS,UKF_x_hat_save_Sub(:,2,5))
plot(TimeGPS,UKF_x_hat_save_Sub(:,2,6))
plot(TimeGPS,UKF_x_hat_save_Sub(:,2,7))
plot(TimeGPS,UKF_x_hat_save_Sub(:,2,8))

% height
figure(); hold on; grid on;
plot(TimeINS,UKF_x_hat_save(:,3),'r')
plot(TimeGPS,UKF_x_hat_save_Sub(:,3,1))
plot(TimeGPS,UKF_x_hat_save_Sub(:,3,2))
plot(TimeGPS,UKF_x_hat_save_Sub(:,3,3))
plot(TimeGPS,UKF_x_hat_save_Sub(:,3,4))
plot(TimeGPS,UKF_x_hat_save_Sub(:,3,5))
plot(TimeGPS,UKF_x_hat_save_Sub(:,3,6))
plot(TimeGPS,UKF_x_hat_save_Sub(:,3,7))
plot(TimeGPS,UKF_x_hat_save_Sub(:,3,8))


%% plot variances
figure(); hold on; grid on;
plot(TimeINS,sqrt(UKF_P_save(:,1))*RMh,'r');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,1,1))*RMh,'b');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,1,2))*RMh,'g');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,1,3))*RMh,'k');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,1,4))*RMh,'m');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,1,5))*RMh,'y');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,1,6))*RMh,'c');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,1,7))*RMh);
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,1,8))*RMh);

%% plot variances
figure(); hold on; grid on;
plot(TimeINS,sqrt(UKF_P_save(:,2))*RPh,'r');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,2,1))*RPh);
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,2,2))*RPh);
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,2,3))*RPh);
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,2,4))*RPh);
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,2,5))*RPh);
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,2,6))*RPh);
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,2,7))*RPh);
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,2,8))*RPh);

%% plot variances
figure(); hold on; grid on;
plot(TimeINS,sqrt(UKF_P_save(:,3)),'r');
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,3,1)));
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,3,2)));
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,3,3)));
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,3,4)));
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,3,5)));
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,3,6)));
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,3,7)));
plot(TimeGPS,sqrt(UKF_P_save_Sub(:,3,8)));

%% calculate separation errors



