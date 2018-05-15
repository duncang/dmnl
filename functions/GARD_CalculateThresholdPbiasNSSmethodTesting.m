function [a_out_H,a_out_V, lambda_out_H, lambda_out_V] = GARD_CalculateThresholdPbiasNSSmethodTesting(P_fa,P_md_H,P_md_V,DOFrange)
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
% $Id: GARD_CalculateThresholdPbias.m 1198 2008-01-30 05:56:52Z bruggema $

%for testing against values given in "Derivation of the RAIM algorithm from
%First principles" by Robert Kelly

feature accel on;


numds = size(DOFrange,2);

% calculate thresholds
%DOF is 2 because of north and east directions
a_H = chi2inv(1-P_fa,DOFrange);   %added ./DOFrange as this is how the threshold is calculated in normalised SS method 29.1.08
TD_H = sqrt(a_H);

%DOF is 1 because of down direction
a_V = chi2inv(1-P_fa./DOFrange,1);   %added ./DOFrange as this is how the threshold is calculated in normalised SS method 29.1.08
TD_V = sqrt(a_V);


a_H = [5.1045, 5.4616,5.7391,5.9750,6.1858,6.3789,6.5576,6.7251,6.8834,7.0338];

a_H = a_H.^2;

%% loop over degrees of freedom
for j = 1:numds
    % for the given threshold, find the non-centrality parameter which
    % meets the P_md

    lambda_guess_H = 1;
    %P_est_H = ncx2cdf(a_H(j),DOFrange(j),lambda_guess_H);
    
    P_est_H = ncx2cdf(a_H(j),DOFrange(j),lambda_guess_H);
    while (P_est_H > P_md_H)
       lambda_guess_H = lambda_guess_H + 0.1;
       %P_est_H = ncx2cdf(a_H(j),DOFrange(j),lambda_guess_H);
       
        P_est_H = ncx2cdf(a_H(j),DOFrange(j),lambda_guess_H);
    end
    
    % save
    lambda_out_H(j) = lambda_guess_H;
    P_out_H(j) = P_est_H;
    
end

a_out_H = a_H;


%% loop over degrees of freedom
for j = 1:numds
    % for the given threshold, find the non-centrality parameter which
    % meets the P_md

    lambda_guess_V = 1;
    %P_est_V = ncx2cdf(a_V(j),DOFrange(j),lambda_guess_V);
    
     P_est_V = ncx2cdf(a_V(j),1,lambda_guess_V);
    while (P_est_V > P_md_V)
       lambda_guess_V = lambda_guess_V + 0.1;
       %P_est_V = ncx2cdf(a_V(j),DOFrange(j),lambda_guess_V);
       
        P_est_V = ncx2cdf(a_V(j),1,lambda_guess_V);
       
    end
    
    % save
    lambda_out_V(j) = lambda_guess_V;
    P_out_V(j) = P_est_V;
    
end

a_out_V = a_V;

