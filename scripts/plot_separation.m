% This script plots the solution separation vectors for a single epoch
% It is designed to work with the GARDSim_GPSINS_* series of scripts.
% Written by Duncan Greer
% $Id: plot_separation.m 1884 2008-07-15 05:54:33Z n2523710 $
%

% plot position

figure();hold on; grid on;
% plot separation vectors
plot(([UKF_x_hat_kplus(2), UKF_x_hat_kplus(2) - UKF_Beta_ss_H(2,1)]-pos_truth_llh(3,Epoch_lo*100+1))*RPh*cos(pos_truth_llh(2,Epoch_lo*100+1)), ...
     ([UKF_x_hat_kplus(1), UKF_x_hat_kplus(1) - UKF_Beta_ss_H(1,1)]-pos_truth_llh(2,Epoch_lo*100+1))*RMh,'g')

error_ellipse([UKF_B_ss_H(2,2,1),UKF_B_ss_H(2,1,1);UKF_B_ss_H(1,2,1),UKF_B_ss_H(1,1,1)]*RMh*RMh,...
              [(UKF_x_hat_kplus(2)-pos_truth_llh(3,Epoch_lo*100+1))*RPh*cos(pos_truth_llh(2,Epoch_lo*100+1)),...
               (UKF_x_hat_kplus(1)-pos_truth_llh(2,Epoch_lo*100+1))*RMh]);

for j=2:8 % for each subfilter
    plot(([UKF_x_hat_kplus(2), UKF_x_hat_kplus(2) - UKF_Beta_ss_H(2,j)]-pos_truth_llh(3,Epoch_lo*100+1))*RPh*cos(pos_truth_llh(2,Epoch_lo*100+1)), ...
         ([UKF_x_hat_kplus(1), UKF_x_hat_kplus(1)- UKF_Beta_ss_H(1,j)]-pos_truth_llh(2,Epoch_lo*100+1))*RMh,'r')
    
error_ellipse([UKF_B_ss_H(2,2,j),UKF_B_ss_H(2,1,j);UKF_B_ss_H(1,2,j),UKF_B_ss_H(1,1,j)]*RMh*RMh,...
              [(UKF_x_hat_kplus(2)-pos_truth_llh(3,Epoch_lo*100+1))*RPh*cos(pos_truth_llh(2,Epoch_lo*100+1)),...
               (UKF_x_hat_kplus(1)-pos_truth_llh(2,Epoch_lo*100+1))*RMh]);


end

% plot full filter
plot((UKF_x_hat_kplus(2)-pos_truth_llh(3,Epoch_lo*100+1))*RPh*cos(pos_truth_llh(2,Epoch_lo*100+1)),...
     (UKF_x_hat_kplus(1)-pos_truth_llh(2,Epoch_lo*100+1))*RMh,'b.') ;

% plot fault-free filter
plot((UKF_x_hat_kplus_Sub(2,1)-pos_truth_llh(3,Epoch_lo*100+1))*RPh*cos(pos_truth_llh(2,Epoch_lo*100+1)),...
     (UKF_x_hat_kplus_Sub(1,1)-pos_truth_llh(2,Epoch_lo*100+1))*RMh,'g.');

% plot faulted filters
plot((UKF_x_hat_kplus_Sub(2,2:8)-pos_truth_llh(3,Epoch_lo*100+1))*RPh*cos(pos_truth_llh(2,Epoch_lo*100+1)),...
    (UKF_x_hat_kplus_Sub(1,2:8)-pos_truth_llh(2,Epoch_lo*100+1))*RMh,'r.');

axis([-5 15 -15 5])
xlabel('East (m)');
ylabel('North (m)')
title(sprintf('Horizontal Solution Separation for the Epoch #%03d',Epoch_lo));
