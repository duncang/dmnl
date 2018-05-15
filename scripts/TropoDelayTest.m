% test TropoDelay function

Elevation = [5:1:90];
UserHeight = [0:1000:7000];

for i=1:length(Elevation)
   for j=1:length(Height)
       TropoDelay(i,j) = GARD_TropoDelay(Elevation(i)*pi/180,UserHeight(j));
   end
end

figure();
mesh(UserHeight/1000,Elevation,TropoDelay)
xlabel('Height (km)');
ylabel('Satellite Elevation (deg)');
zlabel('Tropo Delay (m)');