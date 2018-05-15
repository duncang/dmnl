% plot IFR en-route Fixs on map

if ~exist('AustAirports','var')
    load data/AustAirports.mat;
end

if ~exist('AustFix','var')
    load data/AustFix.mat;
end

if ~exist('AustNav','var')
    load data/AustNav.mat;
end

if ~exist('austcoast')
    load 'data/austcoast.dat';
end

% check for mapping toolbox
if ~exist('geoshow','file')

    
    
    % plot coastline
    plotaustcoast;
    hold on;

    % plot airports
    plot(AustAirports(:,6)*180/pi,AustAirports(:,5)*180/pi,'k*');
    text(AustAirports(:,6)*180/pi+0.01,AustAirports(:,5)*180/pi,char(AustAirports(:,1:4)),'FontSize',6);

    % PLOT IFR Fixes
    plot(AustFix(:,7)*180/pi,AustFix(:,6)*180/pi,'r^')
    text(AustFix(:,7)*180/pi+0.01,AustFix(:,6)*180/pi,char(AustFix(:,1:5)),'FontSize',6);

    % plot NDB
    for i=1:size(AustNav.NDB,2)
        plot(AustNav.NDB(i).data(6)*180/pi,AustNav.NDB(i).data(5)*180/pi,'ro');
        text(AustNav.NDB(i).data(6)*180/pi+0.01,AustNav.NDB(i).data(5)*180/pi,strcat('NDB:',char(AustNav.NDB(i).data(1:4))),'FontSize',6);
    end

    % Plot VOR
    for i=1:size(AustNav.VOR,2)
        plot(AustNav.VOR(i).data(6)*180/pi,AustNav.VOR(i).data(5)*180/pi,'rd');
        text(AustNav.VOR(i).data(6)*180/pi+0.01,AustNav.VOR(i).data(5)*180/pi,strcat('VOR:',char(AustNav.VOR(i).data(1:4))),'FontSize',6);
    end

    hold off;grid on;
else
    % mapping toolbox available
    worldmap([-45 -10],[110 160]);
    % plot coastline
    geoshow(austcoast(:,2),austcoast(:,1));
    
    % plot airports
    geoshow(AustAirports(:,5)*180/pi,AustAirports(:,6)*180/pi,'DisplayType','point','MarkerEdgeColor','black','Marker','*');
    textm(AustAirports(:,5)*180/pi,AustAirports(:,6)*180/pi+0.01,char(AustAirports(:,1:4)),'FontSize',6);
    
    geoshow(AustFix(:,6)*180/pi,AustFix(:,7)*180/pi,'DisplayType','point','MarkerEdgeColor','red','Marker','^')
    textm(AustFix(:,6)*180/pi+0.01,AustFix(:,7)*180/pi,char(AustFix(:,1:5)),'FontSize',6);

    % plot NDB
    for i=1:size(AustNav.NDB,2)
        geoshow(AustNav.NDB(i).data(5)*180/pi,AustNav.NDB(i).data(6)*180/pi,'DisplayType','point','MarkerEdgeColor','red','Marker','o');
        textm(AustNav.NDB(i).data(5)*180/pi+0.01,AustNav.NDB(i).data(6)*180/pi,strcat('NDB:',char(AustNav.NDB(i).data(1:4))),'FontSize',6);
    end

    % Plot VOR
    for i=1:size(AustNav.VOR,2)
        geoshow(AustNav.VOR(i).data(5)*180/pi,AustNav.VOR(i).data(6)*180/pi,'DisplayType','point','MarkerEdgeColor','magenta','Marker','d');
        textm(AustNav.VOR(i).data(5)*180/pi+0.01,AustNav.VOR(i).data(6)*180/pi,strcat('VOR:',char(AustNav.VOR(i).data(1:4))),'FontSize',6);
    end
    
    %tightmap;
    %paperscale(1e7); % 1:10,000,000
    %previewmap;
    
    
    
end