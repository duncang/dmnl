function [xs,W] = GARD_GenerateSigmaPoints(UKF_x_hat_kminus,UKF_Px_kminus,UKF_Q,UKF_R,alpha,beta,kapa)
% function [xs,W] = GARD_GenerateSigmaPoints(UKF_x_hat_kminus,UKF_Px_kminus,UKF_Q,UKF_R,alpha,beta,kapa)
% Written by Duncan Greer 18 Dec 07
%
% $Id: GARD_GenerateSigmaPoints.m 1850 2008-07-14 04:52:47Z greerd $
%
%

NumberStates = length(UKF_x_hat_kminus);
ProcessNoiseStates = size(UKF_Q,1);
MeasurementNoiseStates = size(UKF_R,1);

Na = NumberStates+ProcessNoiseStates+MeasurementNoiseStates;

LAMBDA = alpha^2 * (Na+kapa) - Na;

UKF_xa_hat_kminus = [UKF_x_hat_kminus;zeros(ProcessNoiseStates,1);zeros(MeasurementNoiseStates,1)];

UKF_Pxa_kminus = [UKF_Px_kminus          zeros(NumberStates,ProcessNoiseStates) zeros(NumberStates,MeasurementNoiseStates);
              zeros(ProcessNoiseStates,NumberStates) UKF_Q           zeros(ProcessNoiseStates,MeasurementNoiseStates);
              zeros(MeasurementNoiseStates,NumberStates) zeros(MeasurementNoiseStates,ProcessNoiseStates)  UKF_R];



%blah = sqrt(Na+LAMBDA) * chol(UKF_Pxa_kminus);
blah = sqrt(Na+LAMBDA) * sqrtm(UKF_Pxa_kminus);
for i=0:Na
    if i==0
        xs_0 = UKF_xa_hat_kminus;
        W_0_m = LAMBDA / (Na + LAMBDA);
        W_0_c = W_0_m + (1 - alpha^2 + beta);
    else
        xs_i(:,i) = UKF_xa_hat_kminus + blah(i,:)';
        xs_i(:,i+Na) = UKF_xa_hat_kminus - blah(i,:)';
        W_i_m = 1 / (2 * (Na + LAMBDA));
        W_i_c = W_i_m;
    end
end

xs = [xs_0,xs_i];
W = [W_0_m,W_0_c,W_i_m,W_i_c];