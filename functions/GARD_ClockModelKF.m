function [ClockBias_tTCXO,ClockDrift_tTCXO,ClockBias_tCSAC,ClockDrift_tCSAC] = GARD_ClockModelKF(timeinterval,Delta_tt, U1, U2);
%Version 1.00
%Troy Bruggemann 8 December 2005

%INPUT
%timeinterval - the number of seconds to generate a clock bias and drift for

%OUTPUT
%ClockBias_t - Receiver clock bias in SECONDS
%ClockDrift_t - Receiver clock bias in SECONDS per SECOND

%----------------------------------------------------------------------------%


%This function closely resembles the one by the thesis

%Liang, L., "NUMERICAL SIMULATION OF GPS CODE TRACKING LOOPS," in Aeronautics & Astronautics. West Lafayette: Purdue University, 2004.

%-----------function starts here


 %This is from Chapter 10 'introduction to random signals and applied kalman filtering' second edition RG Brown and P Y C Hwang and Kayton p 88-89


%Nominal Frequency 
fo = 10e6;  % 10.0 MHz @ 25 degrees C (Rakon Data sheet)
 
%to simulate changes of the clock with temperature you could make fo a function of temperature. 


%RootAllanVar = 1/1e9; %1 ppb for tau = 1 second = worst case as given in spec sheet.

 %timeinterval = 1001;
 %t = timeinterval-1;   %Over number of epochs (sample points) to calculate for in seconds.
 t = timeinterval;   %Over number of epochs (sample points) to calculate for in seconds.
 Delta_t = Delta_tt; %this is the time between samples in seconds ie. the step time. this will be tau0 in the allan variance calculation.
 
 
 
 
 
 %for IGNSS
 
 %Typical Crystal Oscillator from Brown
 h_zeroTCXO = 9.43e-20;    %White frequency noise 
 h_minusoneTCXO = 1.8e-19; % Flicker frequency noise
 h_minustwoTCXO = 3.8e-21; %Random walk frequency noise
%  
 
 %CSAC 'best performance'
 h_zeroCSAC = 7.2e-21;    %White frequency noise 
 h_minusoneCSAC =2.59e-23; % Flicker frequency noise
 h_minustwoCSAC = 2.70189e-27; %Random walk frequency noise
%  
%   %CSAC 'worst performance'

%replace these h_zero and h_minusone values with values from brown, see how
%it performs. as long as at tau = 3600, std dev = 1e-11 then it will be
%okay! May need to elevate it up a bit on the h_minustwo so that it passes
%through 1e-11 at tau =3600.


%  h_zeroCSAC = 7.2e-19;    %White frequency noise 
% % h_minusoneCSAC = 7e-24; % Flicker frequency noise
%  h_minustwoCSAC = 4.221715986e-27 ; %Random walk frequency noise
%  

%typical Rb Oscillator from Brown 
 


%Mems atomic clock model


 %Calculate the Sf and Sg terms from the Allan variance parameters given above
 
%  Sfa = 2*h_zero;
%  Sga = 8*pi*pi*h_minustwo;
 
 SfaTCXO = h_zeroTCXO/2;
 SgaTCXO = 2*pi*pi*h_minustwoTCXO;    %this should be on 3 or not? Note that im not using the KF model so i would need to divide by 3 here even though brown doesnt.

 
 SfaCSAC = h_zeroCSAC/2;
 SgaCSAC = 2*pi*pi*h_minustwoCSAC;  
 
%if you look at the matlab function wgn.m u can see if u use in dB 
%its expecting a value calculated with 10*log10

%  Sflog = 10*log10(Sfa);
%  Sglog = 10*log10(Sga);
%  
%  
%  %generate white noise driving functions
%  %note - sg term is for the ud!!
%  
%  ud = wgn(t,1,Sglog); %the drift 
%  ub = wgn(t,1,Sflog); %the bias  
%  
 
 
 %Generate White noise using central limit theorem
% NNoise = 10000;
% 
% 
% for j = 1:t
%     
% 
% 
% X = 0;
%    for i = 1:NNoise
%       U = rand;
%       X = X + U;
%    end
% 	
%    %/* for uniform randoms in [0,1], mu = 0.5 and var = 1/12 */
%    %/* adjust X so mu = 0 and var = 1 */
% 	
%    X = X - NNoise/2       ;         %/* set mean to 0 */
%    X = X * sqrt(12 / NNoise)  ;     %/* adjust variance to 1 */
%    
%    
%    
%    val1(j) = X;
% 
% %When the algorithm finishes, X will be our unit normal random. X can be further modified to have a particular mean and variance, e.g.:
% 
%   % X' = mean + sqrt(variance) * X
% 
% end
% 
% 
% 
% X = 0;
%    for i = 1:NNoise
%       U = rand;
%       X = X + U;
%    end
% 	
%    %/* for uniform randoms in [0,1], mu = 0.5 and var = 1/12 */
%    %/* adjust X so mu = 0 and var = 1 */
% 	
%    X = X - NNoise/2       ;         %/* set mean to 0 */
%    X = X * sqrt(12 / NNoise)  ;     %/* adjust variance to 1 */
%    
%    
%    
%    val2(j) = X;
% 
% %When the algorithm finishes, X will be our unit normal random. X can be further modified to have a particular mean and variance, e.g.:
% 
%   % X' = mean + sqrt(variance) * X
% 
% end

%====================================================================
%Generate white noise using polar method (variation on Box-Muller method)

% for j = 1:t
% 
% S = 2; %initialise it
%   while S >=1
%       U1=rand ;          % /* U1=[0,1] */
%       U2=rand ;       %    /* U2=[0,1] */
%       V1=2 * U1 -1  ;       %   /* V1=[-1,1] */
%       V2=2 * U2 - 1   ;      %  /* V2=[-1,1] */
%       S=V1 * V1 + V2 * V2;
%   end 
% 	
%    val1(j)=sqrt(-2 * log(S) / S) * V1;
%    val2(j)=sqrt(-2 * log(S) / S) * V2;
% 
% end
%====================================================================

 
%two independent white noises, uniformly distributed. zero mean, unit variance.
%randn is normal distributed while rand is uniformly distributed, thats why
%i dont use rand

%generate these outside 
%U1 = rand(t,1);  %generates between 0 and 1. 
U1 = U1 - mean(U1);  %set mean to 0 
U1 = U1*sqrt(12); %adjust variance to 1

%U2 = rand(t,1);
U2 = U2 - mean(U2);  %set mean to 0 
U2 = U2*sqrt(12); %adjust variance to 1

val1 = U1;
val2 = U2;

udTCXO = sqrt(SgaTCXO/Delta_t)*val1; %the drift  
ubTCXO = sqrt(SfaTCXO/Delta_t)*val2; %the bias  

udCSAC = sqrt(SgaCSAC/Delta_t)*val1; %the drift  
ubCSAC = sqrt(SfaCSAC/Delta_t)*val2; %the bias  


  
   %this kalman filter recursive loop is from page 235 figure 5.9 of Brown and Hwang
     
   %Enter Loop
     
  %Initial values - These can be set to an initial bias if desired, but just set to zero here.
  ClockBias_tTCXO(1) = 0; 
  ClockDrift_tTCXO(1) = 0;
  
    
  ClockBias_tCSAC(1) = 0; 
  ClockDrift_tCSAC(1) = 0;
  
  
  %generate measurements

  %with compromise? not considering flicker noise
  
  
  for i = 1:t
      
     ClockBias_tTCXO(i+1) = ClockBias_tTCXO(i) + ClockDrift_tTCXO(i)*Delta_t  + ubTCXO(i)*Delta_t; 
     ClockDrift_tTCXO(i+1) = ClockDrift_tTCXO(i) + udTCXO(i)*Delta_t;        
  
  end     
  
  
    
  for i = 1:t
      
     ClockBias_tCSAC(i+1) = ClockBias_tCSAC(i) + ClockDrift_tCSAC(i)*Delta_t  + ubCSAC(i)*Delta_t; 
     ClockDrift_tCSAC(i+1) = ClockDrift_tCSAC(i) + udCSAC(i)*Delta_t;        
  
  end   
  
  
   
  
  
  
  
  
  
  
  
  
  

  
%------------------------------------------------------
  %considering flicker noise (using the Q and C matrix)
%------------------------------------------------------



% 
% 
%  
%   %generate Q matrix (process noise covariance matrix) for the clock parameters
%     
%   Q(1,1) = (h_zero/2)*Delta_t + 2*h_minusone*Delta_t^2 +(2/3)*pi*pi*h_minustwo*Delta_t^3;
%   Q(2,2) = h_zero/(2*Delta_t) + 4*h_minusone +(8/3)*pi*pi*h_minustwo*Delta_t;
%   Q(1,2) = h_minusone*Delta_t +pi*pi*h_minustwo*Delta_t^2;
%   Q(2,1) = Q(1,2);
%   
%   
% 
%  
%    %generate C matrix (used if flicker noise is to be included in the equations)
%    
%   C(1,1) = sqrt(Q(1,1));
%   C(1,2) = 0;
%   C(2,1) = Q(2,1)/C(1,1);
%   C(2,2) = sqrt(Q(2,2) - C(2,1)^2);     
  
%   
%   Phi = [1 Delta_t;0 1];
%     
%   
%   z(:,1) = [ClockBias_t(1); ClockDrift_t(1)] ; %initial values
%   
%   for i = 1:t      
%        % Noise vector
%        W = [val2(i); val1(i)] ;       %this uses the same white noise as the non zero one
%        z(:,i+1) = Phi*z(:,i) + C*W;    
%   
%   end     
%   
%   
%   ClockBias_tflicker = z(1,:);
%   ClockDrift_tflicker = z(2,:);
%   
  
  


























 
%  %Convert from frequency to time. %Can vary fo according to temperature
% 
% %  ClockDrift_t = (ClockDrift_f/(2*pi*fo));
% %  ClockBias_t = (ClockBias_f/(2*pi*fo));   %should this be negative? 
% %  
% %  
% %  
% %  ClockDrift_t = (ClockDrift_f/fo);
% %  ClockBias_t = (ClockBias_f/fo);   %should this be negative? 
% 
% 
% ClockDrift_t = ClockDrift_f;
% ClockBias_t = ClockBias_f; 


%Calculate Temperature model



%Can relate this to doppler shift due to the clock at the L1.. simple 2ppm
%of 1575.42 MHz..






%Allan variance is a method of representing root mean square random drift error as a function of averaging times.
  





%calculate allan variance for averaging time of 0.01 sec

% time = ClockBias_t;
% difft = diff(time);
% squared = difft.^2;
% allanvar_sim = 0.5*mean(squared);
% %square root of allan variance
% allanvarsqrt_sim = sqrt(allanvar_sim);
% 



%calculate allan variance for averaging time of 0.1 sec






%calculate allan variance for averaging time of 1 sec

% %calculate allan var of the output
% time = dTpos;
% difft = diff(time);
% squared = difft.^2;
% allanvar_sim = 0.5*mean(squared);
% %square root of allan variance
% allanvarsqrt_sim = sqrt(allanvar_sim);




% %write to file for analysis in allanvar
% 
% 
% 
% fid = fopen('clockstufffreq.txt','w');
% fprintf(fid,'%e\n\n\n\n');
% fprintf(fid,'%e\n',ClockBias_t);
% 
% fclose(fid)



%------------- other things



% %Generate WGN using central limit theorem
% 
% 
% N = 100000;
% 
% 
% for j = 1:1000
%     
%     
% 
% X = 0;
%    for i = 1:N
%       U = rand;
%       X = X + U;
%    end
% 	
%    %/* for uniform randoms in [0,1], mu = 0.5 and var = 1/12 */
%    %/* adjust X so mu = 0 and var = 1 */
% 	
%    X = X - N/2       ;         %/* set mean to 0 */
%    X = X * sqrt(12 / N)  ;     %/* adjust variance to 1 */
%    
%    
%    
%    val(j) = X;
% 
% %When the algorithm finishes, X will be our unit normal random. X can be further modified to have a particular mean and variance, e.g.:
% 
%   % X' = mean + sqrt(variance) * X
% 
% 
% end
% 
% 
% 
% 
% Sfa = 100;
% 
% 
% 
% %want to generate white noise with power spectral amplitude of 100. 
% 
% del_t = 0.005;
% F_samp = 1/del_t;
% Ntot = 1000;
% Nseg = 1;
% NFFT = Ntot/Nseg;
% 
% 
% %find power spectral density of the function.
% %need to use this method for a random signal see p. 517,518,519 of signals and systems continous and discrete, ziemer, tranter, fannin
% %this is ia power spectrum estimation function using Welch method courtesy of matlab
% psdval = psd(val,NFFT);
% 
% 
% 
% 
% pda = fft(val);
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% Y = fft(val,1000);
% %The power spectrum, a measurement of the power at various frequencies, is
% 
% Pyy = Y.* conj(Y) / 1000;
% 
% 
% 
% %compare with wgn matlab function
% 
%  h_zero = 2e-18;
% 
%   Sfa = 2*h_zero;
% 
% Sflog = 10*log(Sfa);
%  
%  ud = wgn(10000,1,Sflog);
% 
% Z = fft(ud,1000);
% 
% 
% 
% %the mean of the power spectral density plot gives a value very close to the 
% %value of the input power of the wgn function. Therefore its correct to
% %input the Sf or Sg terms in 10*log10(Sf) in the wgn function. Because we want
% % a psd of constant Sf (approximate by the mean).
% 
% 
% %generate white gaussian noise
% 
%  h_zero = 2e-18;
%  Sfa = 2*h_zero;
%  
%  Sfa =  3.0292e-017;
% 
%  Input = 10*log10(Sfa);
%  ud = wgn(10000,1,Input);
% 
% 
%  %calculate the psd
% psdval = psd(ud,10000);
% 
% 
% %the 10*log10 of the mean of the the psdval should be the same as Sfalog
% 
% 
% meanpsd = mean(psdval); 



