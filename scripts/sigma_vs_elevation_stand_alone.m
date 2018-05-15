

sigma_multipath = 0.22;
sigma_UIRE = 4;
sigma_SNR = 0.22;
sigma_tropo = 0.7;
sigma_URA = sqrt(4); % corresponds to lowest URA index (0)

Elevation = 0:1:90;
for e=1:length(Elevation);
   
   sigma_i(e) = sqrt(sigma_URA^2 + sigma_UIRE^2 + sigma_SNR^2 + (sigma_tropo^2)/(sin(Elevation(e)*pi/180)^2) + (sigma_multipath^2)/(tan(Elevation(e)*pi/180)^2)); 
end

figure(); grid on; hold on;
plot(Elevation,sigma_i);
axis([0 90 0 10]);
xlabel('Elevation Angle (degrees)');
ylabel('\sigma_i (meters)');
