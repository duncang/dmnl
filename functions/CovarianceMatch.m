function [Pmd_Est_H, PrHPE, PrTestStat] = CovarianceMatch(lambda_k,TD_squared,N, DOF,sigma_re, PL, de)


%aim of this function is to calculate whatthe probabiltiy of missed
%detection is. 
%ie to make sure that HPL and VPL sati
%Troy Bruggemann 12.11.08

%note that to make a more sure assessemtn of HPL and VPL, many different
%test scases and scenarios have to be run. 
%the aim of this method is to reduce teh amount of time required to conduct
%simulations, the alternative is to do brute force montecarlo methods
%Note this is for the normalized SS method by Young and McGraw. Will need
%modification as necessary for Multiple solution separatino method by
%Brenner. 

% 
% INPUTS:
% 
% lambda_k  - noncentrality parameter of chi 2 distribtuion, is the maximum test statistic squared
% 
% TD_squared - is the detection threshold squared% 
% N - is the number of satellites
% 
% DOF - is the number of degrees of freedom in chi square distribution
% sigma_re - is the standard deviation of the position error
% 
% PL is the protection level (note , could be VPL)% 
% de  is the horizontal position error due to fault. 

% 
% OUTPUTS:
% Pmd_Est_H - probability taht HPE exceeds HPL  and test stat DOES NOT exceed threshold. 
% 
% PrHPE - probability that HPE exceeds HPL
% 
% PrTestStat - probability that test stat exceesd threshold.


PrTestStat = ncx2pdf(TD_squared,DOF,lambda_k);
PrTestStat = 1 - PrTestStat*N;   %NOTE BECAUSE threshold is calculated as pfa/N, need this here!!!



PrHPE = normcdf(PL-de,0,sigma_re);

PrHPE = 1-PrHPE ;  %do i need this? Because matlab goes from -inf to x, i want tfrom x to inf..i want the tail region on the main region

%PrTestStat(i) = 1-PrTestStat(i);


Pmd_Est_H = PrHPE*PrTestStat; 