function [UKF_x_hat_kminus,UKF_Px_kminus] = GARD_GenerateMeanFromSigmaPoints(xs,Ws,NumberStates)
% function [UKF_x_hat_kminus,UKF_Px_kminus] = GARD_GenerateMeanFromSigmaPoints(xs,Ws,NumberStates)
% Written by Duncan Greer 18 Dec 07
%
% $Id: GARD_GenerateMeanFromSigmaPoints.m 1850 2008-07-14 04:52:47Z greerd $
%
%        


Na = (size(xs,2)-1)/2;
%Na = size(xs,1);

W_0_m = Ws(1);
W_0_c = Ws(2);
W_i_m = Ws(3);
W_i_c = Ws(4);

            
UKF_W = eye(1,2*Na+1);
UKF_W(1,1) = W_0_m;
UKF_W(1,2:2*Na+1) = W_i_m;

UKF_X = xs(1:NumberStates,:);
UKF_x_hat_kminus = (UKF_W*UKF_X')';

nx = length(UKF_x_hat_kminus); nw = length(UKF_W);
UKF_Px_kminus=((ones(nx,1)*UKF_W).*(UKF_X-UKF_x_hat_kminus*ones(1,nw)))*(UKF_X-UKF_x_hat_kminus*ones(1,nw))';