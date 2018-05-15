function [ClockBias_t,ClockDrift_t] = GARD_ClockModel(PrevBias);
%Version 1.00
%Troy Bruggemann 27 September 2005

%This function is to simulate the output of a receiver clock such as a TCXO, based on a given allan variance.
% This function, on given inputs, generates a new clock state based on
% those inputs.
% this function is not the estimated receiver clock model of the receiver, but is to
% simulate a real receiver clock.

%INPUT
%PrevBias = the previous clock bias in SECONDS, this will be updated with the new value

%OUTPUT
%ClockBias_t - Receiver clock bias in SECONDS
%ClockDrift_t - Receiver clock bias in SECONDS per SECOND

%----------------------------------------------------------------------------%

%assumes the root allan variance changes between min and max values with gaussian distribution

%from page 191 of Kayton, figure 5.6

%could possibly add frequency terms for aging , as well. This is just using the stability.
%and make the allan variance a function of temperature, between the min and max values, for example.


RootAllanVarMax = 1/1e9; %1 ppb for tau = 1 second = worst case as given in spec sheet for Rakon TXO200B. Observed as 5/1e9 from output of estimated clock bias from Orion.

RootAllanVarMin = 0.8/1e9; % best case 

a = RootAllanVarMax; b = RootAllanVarMin;  %range of values for random number %note. b has to be the minimum value one.

%randomly varies the root allan variance between minimum and max values (as given by spec sheet for oscillator), gaussian distributed zero mean.

RootAllanVar = a + (b-a) * rand(1);

AllanVar = RootAllanVar^2;

%AllanVar = 5e-9^2;

%amount to add to clock

BiasdiffBias = sqrt(2*AllanVar);
%add random amount to it, standard deviation of 15.5 nanoseconds/second

Biasdiff = BiasdiffBias;

a = 50.5e-9; b = -50.5e-9;  %range of values for random number %note. b has to be the negative value one.

RandomNoise = a + (b-a) * rand(1);

Biasdiff = Biasdiff + RandomNoise;

NewBias = PrevBias + Biasdiff;

%simply differentiate to give a clock drift
%Drift(i+1) = Drift(i) + 0.5*(Bias(i+1) - Bias(i));


%might need to improve this calculation of drift, or calculate drift and derive bias from the drift by integrating maybe.

NewDrift  = 0.5*(NewBias - PrevBias);

%end

ClockBias_t = NewBias;

ClockDrift_t = NewDrift;






%----------Background Theory--------

%want standard deviation of bias to be 50 nanoseconds
%want a ramp up of 5.7 nanoseconds per second

%want drift of Mean = 5.34 nanoseconds/sec
%Std dev = 15.5 nanoseconds/sec

%the gain per second

%stddev_y^2 = 1/2* (delta_y)^2

%if a wrist watch gained one second per day, then y = 1 second / 86400 seconds = 1.157e-5, and J is equal to 1 day






%
% %Way to calculate allan variance
%
%
%
%
%
% RootAllanVar = 1/1e9; %1 ppb for tau = 1 second = worst case as given in spec sheet.
%
%
% %work out allan variance of estimated receiver clock bias
%
%
% %data is in 1 second sample length
%
% time = FinalSolution_Observed(100:200,4)/c; %tau  = 1 second.
%
%
% %average within each sample
%
% difft = diff(time);
%
%
% squared = difft.^2;
%
% allanvar = 0.5*mean(squared);
%
% %square root of allan variance
%
% allanvarsqrt = sqrt(allanvar);


%------------------------------







%------------- Original Clock model Function taken from Kayton -----------------
%TBiasPrev - previous clock state
%FBiasPrev -  Frequency offset
%FDriftPrev
%Clock instability measures.

%type used, eg TCXO etc.
%Temperature - Current Temperature
%Vibration - Current vibration in g's.


%OUTPUT
%Clock Bias (frequency offset)
%Clock Drift (frequency Drift)
%Clock Bias seconds
%Clock Bias seconds/second


%Equation for effect of Temperature on clock
%What is relation between clock instability and temp."


%Equation for effect of vibration on clock
%what is relation between clock and vibration?


%Model Crystal clcok error with 2nd order state model. (from Kayton)



%Parameters used from Rakon TXO200B oscillator (similar to that used in SuperSTar II receivers)

%Delta_f/Fnom = Delta_t/TNom


%Determine noise from clock instability.



%
%
% %Nominal Frequency
% fo = 10e6;  % 10.0 MHz @ 25 degrees C (Rakon Data sheet)
%
% RootAllanVar = 1/1e9; %1 ppb for tau = 1 second = worst case as given in spec sheet.
%
%
% %Values change depending on time interval used..but left as 1000 as this
% %seems close enough to how a real clock behaves.
% %Model Crystal clock errors (see Kayton pp83)
%  t = 0:1000;   %Over timespan to calculate for in seconds.
%  %white noise
%  ud = fo*RootAllanVar*0.5*randn(size(t));   %multiply by fo to convert to Hz.
%
%  %integrate
%  ClockDrift_f = cumtrapz(ud -mean(ud)); %Hz/Second
%  ub = RootAllanVar*0.5*randn(size(t));
%  temp = ClockDrift_f + ub;
%  ClockBias_f = cumtrapz(temp - mean(temp));  %Hz
%
%  %Convert from frequency to time.
%  ClockBias_t = ClockBias_f/fo;
%  ClockDrift_t = ClockDrift_f/fo;



%Can relate this to doppler shift due to the clock at the L1.. simple 2ppm
%of 1575.42 MHz..







