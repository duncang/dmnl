function [Xout1, Xout2] = GaussMarkov_Processfortesting(Xin1,Xin2, Beta1,WN_Std_dev,dT);

%by Troy Bruggemann
% XoutTemp = Ws -(Beta)*Xin; 


randnnn = randn(1);
%XoutTemp = Ws;

%integrate
%Xout = Xin + XoutTemp*dT;

%Xout = Xin*(1-Beta*dT) + WN_Std_dev*randn(1)*dT;


%this is what is calculated if you integrate the x dot:
Xout1 = Xin1*(1-Beta1*dT) + sqrt(2*Beta1*WN_Std_dev^2)*randnnn*dT;


%this is the application of the method in brown and hwang
%see brown and hwang

%zero mean white noise
  %Ws = WN_Std_dev*randn(1)*dT;
  
   Ws = sqrt(2*Beta1*WN_Std_dev^2)*randnnn;
   
   Xout2 = Ws + exp(-Beta1*dT)*Xin2; 
   
   %using either of these gave an autocorrelation function with magnitude
   %of about sigma^2 and value of tau at about 
   %they aren't perfect but will have error due to process not being white
   %or integration error etc. 
   
   %I have found that the first method works better than doing it Brown and
   %HWang way, which is consistently higher than it should be. 