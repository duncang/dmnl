% Ionex Reader Validation

% This function reads an IONEX file and plots hte TEC map for an epoch

Filename = 'data/igsg0050_mod.05i';

[IonoTECMap, GPSTime_Obs, Longitude_out, Latitude_out] = ReadIonex(Filename);


% create movie
% MovieFile = 'data/TEC_Movie.avi';
% SimulationMovie = avifile(MovieFile);

% % plot results
for Epoch = 1%:length(GPSTime_Obs)
    hold on;
    %pcolor(Longitude_out(Epoch,:),Latitude_out(Epoch,:),IonoTECMap(:,:,Epoch));
    contour(Longitude_out(Epoch,:),Latitude_out(Epoch,:),IonoTECMap(:,:,Epoch));
    grid on;
    xlabel('Longitude (deg');
    ylabel('Latitude (deg)');
    title('Ionosphere TECU 5 Jul 2005 from 00:00 UT');
    caxis([0 70]);
    colorbar;
    hold off;
    
%         % take an image frame and add to movie
%     Simulation(Epoch) = getframe(gcf);
%     
%     SimulationMovie = addframe(SimulationMovie,Simulation(Epoch));
end

% SimulationMovie = close(SimulationMovie);