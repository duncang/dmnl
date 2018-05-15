


%% gps range data
clear time;

startindex = 60;

for i=1:size(gps.RangeData,2) - startindex
   time(i,1) = gps.RangeData(i+startindex).rtTimeStamp *1e-9;
   time(i,2) = gps.RangeData(i+startindex).GPSSec;
    
end

rtTimeStart = time(1,1);
gpsTimeStart = time(1,2);

time(:,3) = time(:,1) - rtTimeStart;
time(:,4) = time(:,2) - gpsTimeStart;

% figure(); hold on; grid on;
% plot(time(:,3))
% plot(time(:,4),'r')

drift = time(:,3) - time(:,4);

figure(); grid on; hold on;
plot(time(:,4)/60, drift * 1000);
xlabel('Test Time (mins)');
ylabel('Time Drift (msec)');
title('GPS Time to System Time (1Hz)');

% 
fitmodel = fit(time(:,4),drift,'poly1');

driftrate = fitmodel.p1 * 60 * 1000; % msec / minute

residual = drift - ((fitmodel.p1 * time(:,4)) + fitmodel.p2);

figure();  grid on; hold on;
plot(time(:,4) / 60, residual * 1000);



%% span data

clear spantime;

startindex = 1*50;
stopindex = 25*60*50;

spantime(:,1) = span(startindex:stopindex,2) * 1e-9;
spantime(:,2) = span(startindex:stopindex,4);

rtTimeStart = spantime(1,1);
gpsTimeStart = spantime(1,2);

spantime(:,3) = spantime(:,1) - rtTimeStart;
spantime(:,4) = spantime(:,2) - gpsTimeStart;

% figure(); hold on; grid on;
% plot(time(:,3))
% plot(time(:,4),'r')

spandrift = (spantime(:,3) - spantime(:,4));

figure(); grid on; hold on;
plot(spantime(:,4)/60, spandrift * 1000);
xlabel('Test Time (mins)');
ylabel('Time Drift (msec)');
title('SPAN Time to System Time (50Hz)');





%% all data
figure(); grid on; hold on;
plot(spantime(:,4)/60, (spantime(:,3) - spantime(:,4)) * 1000,'b');
plot(time(:,4)/60, (time(:,3) - time(:,4)) * 1000,'r');
xlabel('Test Time (mins)');
ylabel('Time Drift (msec)');
legend('SPAN (50Hz)', 'GPS (1Hz)');

