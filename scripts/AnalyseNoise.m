%This function computes and plots autocorrelation function for a noise
%signal for 3 gyro and 3 accel. noises.

%troy bruggemann 

%note: reason why sigma^2 of generated noise may be smaller than the signal
%is because my tau value isnt exactly right?? or something to do wiht hte
%timing. The shapes right just not the size of the sigma. 

p_noise = omega_xMODEL100Hz(startepochHighRate+100:endepochHighRate-1500) - GyroTruth100Hz(1,startepochHighRate+100:endepochHighRate-1500);

%Get the autocorrelation function 

maxlags = 1000;
[autocorr,lags] = xcorr(p_noise,maxlags);
autocorr_plot = autocorr/length(p_noise);
plot(lags, autocorr_plot);

p_var = max(autocorr_plot); %variance
value = 0.368*p_var; 
A = autocorr_plot(maxlags+1:length(autocorr_plot));
index_Vec = find( value - 0.001*value < A & A <  value + 0.001*value);
p_tau_check = A(min(index_Vec));
p_tau = min(index_Vec);

%convert to seconds
p_tau_sec = p_tau/100;

%find value along the slope

%there will be inaccuracies due to not being able to model white noise
%properly
p_noise_Gen(startepochHighRate+100) = 0;
for i = startepochHighRate+101: endepochHighRate-1500
[p_noise_Gen(i),variance,randU] = GaussMarkov_Process(p_noise_Gen(i-1),1/p_tau_sec,sqrt(p_var),0.01);
end

maxlags = 1000;
[autocorrCheck,lags] = xcorr(p_noise_Gen,maxlags);
autocorr_plotCheck = autocorr/length(autocorrCheck);

plot(lags, autocorr_plot); title ('Autocorr of blue - p noise rad/s, red - modelled noise GM1 rad/s')
hold;
plot(lags, autocorr_plotCheck,'r'); 

pause;
close all;



q_noise = omega_yMODEL100Hz(startepochHighRate+100:endepochHighRate-1500) - GyroTruth100Hz(2,startepochHighRate+100:endepochHighRate-1500);

%Get the autocorrelation function 

maxlags = 1000;
[autocorr,lags] = xcorr(q_noise,maxlags);
autocorr_plot = autocorr/length(q_noise);
plot(lags, autocorr_plot);

q_var = max(autocorr_plot); %variance
value = 0.368*q_var; 
A = autocorr_plot(maxlags+1:length(autocorr_plot));
index_Vec = find( value - 0.01*value < A & A <  value + 0.01*value);
q_tau_check = A(min(index_Vec));
q_tau = min(index_Vec);

%convert to seconds
q_tau_sec = q_tau/100;

%find value along the slope

%there will be inaccuracies due to not being able to model white noise
%properly
q_noise_Gen(startepochHighRate+100) = 0;
for i = startepochHighRate+101: endepochHighRate-1500
[q_noise_Gen(i),variance,randU] = GaussMarkov_Process(q_noise_Gen(i-1),1/q_tau_sec,sqrt(q_var),0.01);
end

maxlags = 1000;
[autocorrCheck,lags] = xcorr(q_noise_Gen,maxlags);
autocorr_plotCheck = autocorr/length(autocorrCheck);

plot(lags, autocorr_plot); title ('Autocorr of blue - q noise rad/s, red - modelled noise GM1 rad/s')
hold;
plot(lags, autocorr_plotCheck,'r'); 

pause;
close all;

 

%r noise

r_noise = omega_zMODEL100Hz(startepochHighRate+100:endepochHighRate-1500) - GyroTruth100Hz(3,startepochHighRate+100:endepochHighRate-1500);

%Get the autocorrelation function 

maxlags = 20000;
[autocorr,lags] = xcorr(r_noise,maxlags);
autocorr_plot = autocorr/length(r_noise);
plot(lags, autocorr_plot);

r_var = max(autocorr_plot); %variance
value = 0.368*r_var; 
A = autocorr_plot(maxlags+1:length(autocorr_plot));
index_Vec = find( value - 0.01*value < A & A <  value + 0.01*value);
r_tau_check = A(min(index_Vec));
r_tau = min(index_Vec);

%convert to seconds
r_tau_sec = r_tau/100;

%find value along the slope

%there will be inaccuracies due to not being able to model white noise
%properly
r_noise_Gen(startepochHighRate+100) = 0;
for i = startepochHighRate+101: endepochHighRate-1500
[r_noise_Gen(i),variance,randU] = GaussMarkov_Process(r_noise_Gen(i-1),1/r_tau_sec,sqrt(r_var),0.01);
end

maxlags = 1000;
[autocorrCheck,lags] = xcorr(r_noise_Gen,maxlags);
autocorr_plotCheck = autocorr/length(autocorrCheck);

plot(lags, autocorr_plot); title ('Autocorr of blue - r noise rad/s, red - modelled noise GM1 rad/s')
hold;
plot(lags, autocorr_plotCheck,'r'); 

pause;
close all;



%Ax



%ax noise

ax_noise = ax_b_MODEL100Hz(startepochHighRate+100:endepochHighRate-1500) - AccelTruth100Hz(1,startepochHighRate+100:endepochHighRate-1500);

%Get the autocorrelation function 

maxlags = 20000;
[autocorr,lags] = xcorr(ax_noise,maxlags);
autocorr_plot = autocorr/length(ax_noise);
plot(lags, autocorr_plot);

ax_var = max(autocorr_plot); %variance
value = 0.368*ax_var; 
A = autocorr_plot(maxlags+1:length(autocorr_plot));
index_Vec = find( value - 0.01*value < A & A <  value + 0.01*value);
ax_tau_check = A(min(index_Vec));
ax_tau = min(index_Vec);

%convert to seconds
ax_tau_sec = ax_tau/100;

%find value along the slope

%there will be inaccuracies due to not being able to model white noise
%properly
ax_noise_Gen(startepochHighRate+100) = 0;
for i = startepochHighRate+101: endepochHighRate-1500
[ax_noise_Gen(i),variance,randU] = GaussMarkov_Process(ax_noise_Gen(i-1),1/ax_tau_sec,sqrt(ax_var),0.01);
end

maxlags = 1000;
[autocorrCheck,lags] = xcorr(ax_noise_Gen,maxlags);
autocorr_plotCheck = autocorr/length(autocorrCheck);

plot(lags, autocorr_plot); title ('Autocorr of blue - ax noise m/s^2, red - modelled noise GM1 m/s^2')
hold;
plot(lags, autocorr_plotCheck,'r'); 

pause;
close all;




%ay noise

ay_noise = ay_b_MODEL100Hz(startepochHighRate+100:endepochHighRate-1500) - AccelTruth100Hz(2,startepochHighRate+100:endepochHighRate-1500);

%Get the autocorrelation function 

maxlags = 20000;
[autocorr,lags] = xcorr(ay_noise,maxlags);
autocorr_plot = autocorr/length(ay_noise);
plot(lags, autocorr_plot);

ay_var = max(autocorr_plot); %variance
value = 0.368*ay_var; 
A = autocorr_plot(maxlags+1:length(autocorr_plot));
index_Vec = find( value - 0.01*value < A & A <  value + 0.01*value);
ay_tau_check = A(min(index_Vec));
ay_tau = min(index_Vec);

%convert to seconds
ay_tau_sec = ay_tau/100;

%find value along the slope

%there will be inaccuracies due to not being able to model white noise
%properly
ay_noise_Gen(startepochHighRate+100) = 0;
for i = startepochHighRate+101: endepochHighRate-1500
[ay_noise_Gen(i),variance,randU] = GaussMarkov_Process(ay_noise_Gen(i-1),1/ay_tau_sec,sqrt(ay_var),0.01);
end

maxlags = 1000;
[autocorrCheck,lags] = xcorr(ay_noise_Gen,maxlags);
autocorr_plotCheck = autocorr/length(autocorrCheck);

plot(lags, autocorr_plot); title ('Autocorr of blue - ay noise m/s^2, red - modelled noise GM1 m/s^2')
hold;
plot(lags, autocorr_plotCheck,'r'); 

pause;
close all;





%az noise

az_noise = az_b_MODEL100Hz(startepochHighRate+100:endepochHighRate-1500) - AccelTruth100Hz(3,startepochHighRate+100:endepochHighRate-1500);

%Get the autocorrelation function 

maxlags = 2000;
[autocorr,lags] = xcorr(az_noise,maxlags);
autocorr_plot = autocorr/length(az_noise);
plot(lags, autocorr_plot);

az_var = max(autocorr_plot); %variance
value = 0.368*az_var; 
A = autocorr_plot(maxlags+1:length(autocorr_plot));
index_Vec = find( value - 0.01*value < A & A <  value + 0.01*value);
az_tau_check = A(min(index_Vec));
az_tau = min(index_Vec);

%convert to seconds
az_tau_sec = az_tau/100;

%find value along the slope

%there will be inaccuracies due to not being able to model white noise
%properly
az_noise_Gen(startepochHighRate+100) = 0;
for i = startepochHighRate+101: endepochHighRate-1500
[az_noise_Gen(i),variance,randU] = GaussMarkov_Process(az_noise_Gen(i-1),1/az_tau_sec,sqrt(az_var),1);
end

maxlags = 1000;
[autocorrCheck,lags] = xcorr(az_noise_Gen,maxlags);
autocorr_plotCheck = autocorr/length(autocorrCheck);

plot(lags, autocorr_plot); title ('Autocorr of blue - az noise m/s^2, red - modelled noise GM1 m/s^2')
hold;
plot(lags, autocorr_plotCheck,'r'); 

pause;
close all;





%display


p_std = sqrt(p_var)
p_tau_sec

q_std = sqrt(q_var)
q_tau_sec

r_std = sqrt(r_var)
r_tau_sec

ax_std = sqrt(ax_var)
ax_tau_sec

ay_std = sqrt(ay_var)
ay_tau_sec

az_std = sqrt(az_var)
az_tau_sec






randn('state',2);
az_tau_sec = 300;

NumberPoints = 1000;

az_var = 25;

dT = 1;

randU(1) = randn(1);
az_noise_Gen(1) = sqrt(az_var)*randU(1)*dT; 
for i = 2: NumberPoints
[az_noise_Gen(i) Variance(i),randU(i)] = GaussMarkov_Process(az_noise_Gen(i-1),1/az_tau_sec,sqrt(az_var),dT);
end


maxlags = NumberPoints;
[autocorrCheck,lags] = xcorr(az_noise_Gen,maxlags);
autocorrPlotCheck = autocorrCheck/(length(az_noise_Gen));  %scale the magnitude

%plot(lags, autocorr_plot); title ('Autocorr of blue - az noise m/s^2, red - modelled noise GM1 m/s^2')
%hold;

% 
% %rescale the autocorrelation for the different time step
% for i = 1:length(autocorrPlotCheck)/(1/dT)
%     
%     autocorrPlotCheckScaled(i) = autocorrPlotCheck
%     
% end
%     



plot(lags/(1/dT), autocorrPlotCheck,'r');   %x axis is seconds, as scaled by 1/dT 




da = Rx_est(az_noise_Gen,50000);







[CurveFit, BestFit,  normX, Xest] = LeastSquaresBestFitGMProcess(az_noise_Gen, 5, dT, 0,randU);






