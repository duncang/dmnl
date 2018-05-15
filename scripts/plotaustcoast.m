% to plot teh australian coastline map
figure();
if ~exist('austcoast')
    load 'data/austcoast.dat';
end
plot(austcoast(:,1),austcoast(:,2));
grid on;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');