x = [-10:0.01:10];

Px = normpdf(x,0,1);

Px1 = normpdf(x-2,0,1);

aN = abs(norminv(0.01,0,1))+2;

%% create a fill region for pmd
x_md = [aN:0.01:10];
P_md = normpdf(x_md-2,0,1);

x_fa = [2:0.1:10];
P_fa = normpdf(x_fa,0,1);
% figure;
% plot(x,Px);
% hold on;
% plot(x,Px1,'r');
% plot([0 0], [0 0.5],'b','LineWidth',1.5);
% plot([2 2], [0 0.5],'b','LineWidth',1.5);
% plot([aN aN], [0 0.5],'r','LineWidth',1.5);
% axis([-5 7 0 0.5]);
% %fill(x_md,P_md,'r');
% grid on;
% xlabel('Test Statistic (m)');
% ylabel('P(Test Statistic)');

createfigure_pl;