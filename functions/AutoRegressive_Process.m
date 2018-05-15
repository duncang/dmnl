function [Xout] = AutoRegressive_Process(Xin,A,WN_Std_dev);


%This process approximates GPS pseudorange error based on a simple first
%order autoregressive model (AR(1)).
%see Integrity monitoring for multi-sensor integrated navigation systems
%(Igor Nikiforov)

%A is the auoregressive coefficient of AR(1). 

Ws = WN_Std_dev*randn(1);
Xout = A*Xin + sqrt(1-A^2)*Ws;

