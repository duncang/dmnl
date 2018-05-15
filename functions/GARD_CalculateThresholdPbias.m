function [a_out, lambda_out] = GARD_CalculateThresholdPbias(P_fa,P_md,DOFrange)
% function [a_out, lambda_out] = GARD_CalculateThresholdPbias(P_fa,P_md,DOFrange)
% 
%
% This code calculates the values for a values (used to calculate) threshold and pbias
% for different degrees of freedom. and PfalseAlarm and Pmiss 
% input the PFalseAlarm and Pmiss and the range of degrees' of freedoms you
% would like to calculate the a and pbias for.
% outputs are stored in the a vector and the lambdatrue vector
%BEWARE: Make sure the DOF are modified appropriately in a, for the purpose
%intended. The DOF are different between solution separation methods and
%the parity method.
% $Id: GARD_CalculateThresholdPbias.m 2792 2009-08-16 04:52:00Z greerd $

feature accel on;



numds = size(DOFrange,2);

% calculate thresholds
a = chi2inv(1-P_fa./DOFrange,DOFrange);   %added ./DOFrange as this is how the threshold is calculated in normalised SS method 29.1.08
TD = sqrt(a);

%% loop over degrees of freedom
for j = 1:numds
    % for the given threshold, find the non-centrality parameter which
    % meets the P_md

    disp(sprintf('calculating for DOF=%d',j));
    
    lambda_guess = 1;
    P_est = ncx2cdf(a(j),DOFrange(j),lambda_guess);
    while (P_est > P_md)
       lambda_guess = lambda_guess + 0.1;
       P_est = ncx2cdf(a(j),DOFrange(j),lambda_guess);
       
    end
    
    % save
    lambda_out(j) = lambda_guess;
    P_out(j) = P_est;
    
end

a_out = a;
