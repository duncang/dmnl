
clear all
i = 2

dtINS = 0.01;
AN1(1,1) = 0;
AN2(1,1) = 0;
%tau = 100; %300 secs



 %-0.0512932943875506

 
 %tau = 1/0.0512932943875506;  %this is 19.45 sec

 tau = 300
stddev = 2;
for i = 2:20000
    
    
     [AN1(1,i), AN2(1,i)] = GaussMarkov_Process(AN1(1,i-1),AN2(1,i-1), 1/tau,stddev,dtINS );
          %[AN(2,i)] = GaussMarkov_Process(AN(2,i-1),1/tau,stddev,dtINS );
     
end


maxlags = 200;

[autocorr1,lags] = xcorr(AN1(1,:),maxlags);
[autocorr2,lags] = xcorr(AN2(1,:),maxlags);

%[autocorr,lags] = xcorr(AN(1,:),maxlags,'coeff');

%need to scale 


autocorr1p = autocorr1/length(AN1(1,:));
autocorr2p = autocorr2/length(AN2(1,:));


plot(lags, autocorr2p);


plot(lags, autocorr1/length(AN1)*100);



%try this algorithm for doing the autocovariance
figure();
AC = Rx_est(AN(1,:),100);

plot(AC);







%=============latest===========

i = 2;


dtINS = 0.01;

stddev = 0.53;

AN1(1,1) = 0.53*randn(1)*dtINS;

AN2(1,1) = 0.013*randn(1)*dtINS;

%beta = 1/300;


beta = 1/300;
%for i = 2:200000
   
for i = 2:10000 
    
     [AN1(1,i),nah, nnah2] = GaussMarkov_Process(AN1(1,i-1), beta,stddev,dtINS );
         
     %[AN1(1,i), AN2(1,i)] = GaussMarkov_Processfortesting(AN1(1,i-1),AN2(1,i-1),1/300,stddev,dtINS );
end


%extract samples


    
for j = 1:180
    
 
    
    AN11(j) = AN1(1,j*100+1)

    
    
end






maxlags = 200000;

[autocorr1,lags] = xcorr(AN1(1,:)*100,maxlags);

autocorr1p = autocorr1/length(AN1(1,:));

%autocorr1p = autocorr1;


plot(lags, autocorr1p);



maxlags = 100000;

[autocorr2,lags] = xcorr(AN2(1,:),maxlags);

autocorr2p = autocorr2/length(AN2(1,:));

%autocorr1p = autocorr1;


plot(lags, autocorr2p);






maxlags = 2000;

[autocorr1,lags] = xcorr(AN1(1,:),maxlags);

autocorr1p = autocorr1/length(AN1(1,:));

%autocorr1p = autocorr1;


plot(lags, autocorr1p);


figure();
AC = Rx_est(AN1(1,:),2000);

plot(AC);





maxlags = 3000;

[autocorr1,lags] = xcorr(GPSPr_noise,maxlags);

autocorr1p = autocorr1/length(GPSPr_noise);
plot(lags, autocorr1p);


figure();
AC = Rx_est(GPSPr_noise,100);

plot(AC);








%to test white noise process

maxlags = 200;

[autocorr1,lags] = xcorr(Wk,maxlags);

%[autocorr,lags] = xcorr(AN(1,:),maxlags,'coeff');

%need to scale 

autocorr1p = autocorr1/length(Wk);



plot(lags, autocorr1p);










