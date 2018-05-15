function [UKF_y_hat_kminus,UKF_Py_kminus,UKF_Pxy_kminus] = GARD_GenerateMeasurementMeanFromSigmaPoints(xs,ys,Ws,NumberStates,NumberMeasurements)
% 
% Written by Duncan Greer 18 Dec 07
%
% $Id: GARD_GenerateMeasurementMeanFromSigmaPoints.m 1850 2008-07-14 04:52:47Z greerd $
%
%        
%Na = size(xs,1);
Na = (size(xs,2)-1)/2;
W_0_m = Ws(1);
W_0_c = Ws(2);
W_i_m = Ws(3);
W_i_c = Ws(4);



ys_kminus_0 = ys(:,1);
ys_kminus_i = ys(:,2:2*Na+1);

UKF_W = zeros(1,2*Na+1);
UKF_W(1,1) = W_0_m;
UKF_W(1,2:2*Na+1) = W_i_m;


UKF_X = xs(1:NumberStates,:);
UKF_x_hat_kminus = (UKF_W*UKF_X')';


UKF_Y = ys;
UKF_y_hat_kminus = (UKF_W*UKF_Y')';



nx = NumberStates; ny = length(UKF_y_hat_kminus); nw = length(UKF_W);
UKF_Py_kminus=((ones(ny,1)*UKF_W).*(UKF_Y-UKF_y_hat_kminus*ones(1,nw)))*(UKF_Y-UKF_y_hat_kminus*ones(1,nw))';
UKF_Pxy_kminus =((ones(nx,1)*UKF_W).*(UKF_X-UKF_x_hat_kminus*ones(1,nw)))*(UKF_Y-UKF_y_hat_kminus*ones(1,nw))';

