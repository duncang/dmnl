

sigma_multipath = 0.22;
sigma_UDRE = 0.5;
sigma_UIVE = 0.5;
sigma_SNR = 0.22;
sigma_tropo = 0.15;

Elevation = 0:1:90;
for e=1:length(Elevation);
   
   sigma_i(e) = sqrt(sigma_UDRE^2 + (GARD_IonoSlantFactor(Elevation(e)*pi/180)^2)*sigma_UIVE^2 + sigma_SNR^2 + (sigma_tropo^2)/(sin(Elevation(e)*pi/180)^2) + (sigma_multipath^2)/(tan(Elevation(e)*pi/180)^2)); 
end

figure(); grid on; hold on;
plot(Elevation,sigma_i);
axis([0 90 0 5]);
xlabel('Elevation Angle (degrees)');
ylabel('\sigma_i (meters)');
