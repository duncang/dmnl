% 
% chi_dist_plot.m
% Written by Duncan Greer
% This script plots a chi-squared distribution with various degrees of
% freedom.
%
% $Id: chi_dist_plot.m 905 2007-11-18 23:49:26Z greerd $
%
x = [0:0.01:50];

sigma = 1;

for DOF=1:10
    chipdf(DOF,:) = sigma * chi2pdf(x,DOF);
end


figure();
plot(x,chipdf(1:5,:)');
legend('1 DOF','2 DOF','3 DOF','4 DOF','5 DOF')
grid on;
axis([0 20 0 0.4]);
xlabel('x');
ylabel('P(x)');
title('\chi^{2} Distribution');

%% non-central
DELTA = 3;
for DOF=1:10
    ncchipdf(DOF,:) = sigma * ncx2pdf(x,DOF,DELTA);
end

figure();
plot(x,ncchipdf(1:5,:)');
legend('1 DOF','2 DOF','3 DOF','4 DOF','5 DOF')
grid on;
axis([0 20 0 0.4]);
xlabel('x');
ylabel('P(x)');
title('Non-Central \chi^{2} Distribution, \lambda = 3');