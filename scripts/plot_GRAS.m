% plot the australian coast, australian airports, GRS and GVS station coverages.


r2d =180/pi;


% generate the 95nm (5000') and 135 nm (10,000')  range contours
[x_40,y_40,z_40] = cylinder(((40*1.852)/6378)*180/pi);
[x_95,y_95,z_95] = cylinder(((95*1.852)/6378)*180/pi);
[x_135,y_135,z_135] = cylinder(((135*1.852)/6378)*180/pi);


GVS = AustAirports(1:100,5:6);
NumberGVS = length(GVS);


if ~exist('austcoast')
    load 'data/austcoast.dat';
end

if ~exist('AustAirports')
    load 'data\AustAirports.mat';
end


% plots the 5000' coverage
figure();
hold on;
for(GVS_i = 1:NumberGVS)
    fill(x_95(1,:)+GVS(GVS_i,2)*180/pi,y_95(1,:)+GVS(GVS_i,1)*180/pi,'g','EdgeColor','none'); 
end
% plot the australian airports
plot(AustAirports(:,6)*r2d,AustAirports(:,5)*r2d,'r.');
plot(austcoast(:,1),austcoast(:,2));
grid on;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('GRAS VHF Coverage at 5000 feet');
hold off;

% plots the 10,000' coverage
figure();
hold on;
for(GVS_i = 1:NumberGVS)
    fill(x_135(1,:)+GVS(GVS_i,2)*180/pi,y_135(1,:)+GVS(GVS_i,1)*180/pi,'m','EdgeColor','none'); 
end
for(GVS_i = 1:NumberGVS)
    fill(x_95(1,:)+GVS(GVS_i,2)*180/pi,y_95(1,:)+GVS(GVS_i,1)*180/pi,'g','EdgeColor','none'); 
end
for(GVS_i = 1:NumberGVS)
    fill(x_40(1,:)+GVS(GVS_i,2)*180/pi,y_40(1,:)+GVS(GVS_i,1)*180/pi,'c','EdgeColor','none'); 
end
%legend('10,000','5000','1000');

% plot the australian airports
plot(AustAirports(:,6)*r2d,AustAirports(:,5)*r2d,'r.');
plot(austcoast(:,1),austcoast(:,2));
grid on;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');
title('GRAS VHF Coverage at 10,000, 5000 and 1000 feet');
hold off;



