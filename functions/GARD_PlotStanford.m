function fh = GARD_PlotStanford(error,protection_level,alert_limit,resolution,max_error,figure_title,use_colour)
% fh =
% GARD_PlotStanford(error,protection_level,alert_limit,resolution,max_error,figure_title,use_colour)
%
% plot stanford chart using histogram and intensity map
%

mhist = hist2d([error, protection_level],[0:resolution:max_error],[0:resolution:max_error]);

% convert to percentage
mhist = (mhist / sum(sum(mhist)))*100;

fh = pcolor(mhist');
shading flat;
caxis([0 2]);
%fh = image(mhist');
%fh = contourf(mhist',10);

if ~exist('use_colour')
    use_colour = 0;
end

if use_colour == 1
    jetmod = jet;
    jetmod(1,1:3) = 1;
    
    cmap = jetmod;
else
    graymod(:,1) = (log([1:64])/log(64));
    graymod(:,2) = (log([1:64])/log(64));
    graymod(:,3) = (log([1:64])/log(64));
    graymod(1,1:3) = 1;
    cmap = graymod;
    
    %cmap = 1-gray(max_error/resolution);

end


colormap(cmap);
colorbar
%grid off;
%set(fh,'Grid','off');
hold on;
%plot(error,protection_level,'b.');
plot([0 max_error],[alert_limit alert_limit],'k','LineWidth',2);
plot([alert_limit alert_limit],[0 alert_limit],'k','LineWidth',2);
axis([0 max_error 0 max_error]);
text(alert_limit/2 + 5,alert_limit/2-5,'MI');
text(alert_limit+(max_error-alert_limit)/2,alert_limit/2 - 5,'HMI');
text(alert_limit+(max_error-alert_limit)/2,alert_limit+(max_error-alert_limit)/2 - 10,'MI');
plot([0:max_error],[0:max_error],'k-','LineWidth',2);
xlabel('Position Error (m)');
ylabel('Protection Level (m)');
grid off;
title(figure_title);

