function PR = DualFreqIono(PR1,PR2);
% function to calculate the dual-frequency based IONO correction
% Written by Duncan Greer 1 August 2005
%
% function PR = DualFreqIono(PR1,PR2);
%
% PR1 is the pseudorange measured on L1
% PR2 is the pseudorange measured on L2
% PR is the corrected pseudorange

% define gamma
L1_f = 1575.42e6; %Hz
L2_f = 1227.6e6; %Hz

gamma = (L1_f/L2_f)^2;  % unitless

% equation from ICD200 para 20.3.3.3.3.3
PR = (PR2 - gamma*PR1)/(1-gamma);
