function [Xout,variance,randU] = GaussMarkov_Process(Xin, Beta,WN_Std_dev,dT);

%by Troy Bruggemann

%This generates sequence of discremte samples of a Gauss -markov process
%IMPORTANT - must multiply the result by 1/dT so the magnitude is right!!!

randU = randn(1);  %need to use randn not rand I think

%scale by the variance to make gaussian random sample with variance of
%sigma^2(1-e^(-2*beta*dt) and mean 0

variance = (WN_Std_dev^2)*(1-exp(-2*Beta*dT));


Wk = sqrt(variance)*randU*dT;   % 15.3.08, Troy - added *dT  - note, must * result Xout by 1/dT 


Xout = exp(-Beta*dT)*Xin + Wk;





%test white noise generation
% 
% WN_Std_dev = 3; 
% dT = 1;
% Beta = 10000;
% for i = 1:1000
%     
%    randU = randn(1); %need to use randn not rand
% 
% 
% %scale by the variance to make gaussian random sample with variance of
% %sigma^2(1-e^(-2*beta*dt) and mean 0
% 
% 
%     variance = WN_Std_dev^2*(1-exp(-2*Beta*dT));
% 
% 
%     Wk(i) = sqrt(variance)*randU;
% 
% 
% end


    
    
    
    


%XoutTemp = Ws;

%integrate
%Xout = Xin + XoutTemp*dT;

%Xout = Xin*(1-Beta*dT) + WN_Std_dev*randn(1)*dT;


%this is what is calculated if you integrate the x dot:
%Xout = Xin*(1-Beta*dT) + sqrt(2*Beta*WN_Std_dev^2)*randnnn*dT;


%this is the application of the method in brown and hwang
%see brown and hwang

%zero mean white noise
  %Ws = WN_Std_dev*randn(1)*dT;
  
  % Ws = sqrt(2*Beta1*WN_Std_dev^2)*randnnn;
   
%   Xout2 = Ws + exp(-Beta1*dT)*Xin2; 



   
   %using either of these gave an autocorrelation function with magnitude
   %of about sigma^2 and value of tau at about 
   %they aren't perfect but will have error due to process not being white
   %or integration error etc. 
   
   %I have found that the first method works better than doing it Brown and
   %HWang way, which is consistently higher than it should be. 