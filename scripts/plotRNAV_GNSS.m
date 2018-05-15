
% the waypoint data is contained in the 'approach' data struct

if ~exist('approach')
    run data/YBBN_RNAV_GNSS_RWY14.m
end

if ~exist('austcoast')
    load 'data/austcoast.dat';
end
    
num_wpts = size(approach.WPTTable,1);

min_lat = min(approach.WPTTable(:,1))*180/pi;
min_lon = min(approach.WPTTable(:,2))*180/pi;
max_lat = max(approach.WPTTable(:,1))*180/pi;
max_lon = max(approach.WPTTable(:,2))*180/pi;

% check for mapping toolbox
if ~exist('geoshow','file')

    warning('Mapping toolbox not available');
    
    figure(); hold on; grid on;
    for i=1:num_wpts
       plot(approach.WPTTable(i,2)*180/pi, approach.WPTTable(i,1)*180/pi,'rh');
       text(approach.WPTTable(i,2)*180/pi+0.005, approach.WPTTable(i,1)*180/pi+0.005,approach.WPTName(i,:));

       if strcmp(approach.IAF,approach.WPTName(i,:)) 
           IAF = i; 
       end

       if strcmp(approach.IF,approach.WPTName(i,:))
           IF = i; 
       end

       if strcmp(approach.FAF,approach.WPTName(i,:))
           FAF = i; 
       end

       if strcmp(approach.MAPWP,approach.WPTName(i,:))
           MAWP = i; 
       end
       if strcmp(approach.MAPHP,approach.WPTName(i,:))
           MAHP = i; 
       end
    end

    plot([approach.WPTTable(IAF,2)*180/pi, approach.WPTTable(IF,2)*180/pi],...
         [approach.WPTTable(IAF,1)*180/pi, approach.WPTTable(IF,1)*180/pi],'k-','LineWidth',2);
    plot([approach.WPTTable(IF,2)*180/pi, approach.WPTTable(FAF,2)*180/pi],...
         [approach.WPTTable(IF,1)*180/pi, approach.WPTTable(FAF,1)*180/pi],'k-','LineWidth',2);
    plot([approach.WPTTable(FAF,2)*180/pi, approach.WPTTable(MAWP,2)*180/pi],...
         [approach.WPTTable(FAF,1)*180/pi, approach.WPTTable(MAWP,1)*180/pi],'k-','LineWidth',2);
    plot([approach.WPTTable(MAWP,2)*180/pi, approach.WPTTable(MAHP,2)*180/pi],...
         [approach.WPTTable(MAWP,1)*180/pi, approach.WPTTable(MAHP,1)*180/pi],'k--','LineWidth',2);


    % to calculate hte distance between wayponits
    %[brg,dist] = CalculateGC(approach.WPTTable(1,1),approach.WPTTable(1,2),approach.WPTTable(4,1),approach.WPTTable(4,2))



     % to plot the australian coastline map


    plot(austcoast(:,1),austcoast(:,2));
    grid on;
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');



    axis([min_lon-0.02 max_lon+0.06 min_lat-0.02 max_lat+0.06]);

    plot(pos_truth_llh(3,:)*180/pi,pos_truth_llh(2,:)*180/pi,'r-','LineWidth',1.5)

    daspect([1 1 1]);

    if exist('AustAirports')
        num_airps = size(AustAirports,1);
        for i=1:num_airps
            if strcmp(char(AustAirports(i,1:4)),approach.PLACE)
                airp_i = i;
                break;
            end
        end

        plot(AustAirports(airp_i,6)*180/pi,AustAirports(airp_i,5)*180/pi,'bd');
        text(AustAirports(airp_i,6)*180/pi+0.005,AustAirports(airp_i,5)*180/pi+0.005,char(AustAirports(airp_i,1:4)));


    end
else
    disp('Using mapping toolbox');
    
    % mapping toolbox available
    worldmap([min_lat-0.02 max_lat+0.06],[min_lon-0.02 max_lon+0.06]);
    gridm;
    
    % plot coastline
    geoshow(austcoast(:,2),austcoast(:,1));
    
    
    geoshow(pos_truth_llh(2,:)*180/pi,pos_truth_llh(3,:)*180/pi,'DisplayType','line','Color','black')
    %textm(AustFix(:,6)*180/pi+0.01,AustFix(:,7)*180/pi,char(AustFix(:,1:5)),'FontSize',6);

    for i=1:num_wpts
       geoshow(approach.WPTTable(i,1)*180/pi, approach.WPTTable(i,2)*180/pi,'DisplayType','point','Marker','d','MarkerEdgeColor','black');
       textm(approach.WPTTable(i,1)*180/pi+0.005, approach.WPTTable(i,2)*180/pi+0.005,approach.WPTName(i,:));

       if strcmp(approach.IAF,approach.WPTName(i,:)) 
           IAF = i; 
       end

       if strcmp(approach.IF,approach.WPTName(i,:))
           IF = i; 
       end

       if strcmp(approach.FAF,approach.WPTName(i,:))
           FAF = i; 
       end

       if strcmp(approach.MAPWP,approach.WPTName(i,:))
           MAWP = i; 
       end
       if strcmp(approach.MAPHP,approach.WPTName(i,:))
           MAHP = i; 
       end
    end

     geoshow([approach.WPTTable(IAF,1)*180/pi, approach.WPTTable(IF,1)*180/pi],...
             [approach.WPTTable(IAF,2)*180/pi, approach.WPTTable(IF,2)*180/pi],'LineWidth',2,'Color','black');
     geoshow([approach.WPTTable(IF,1)*180/pi, approach.WPTTable(FAF,1)*180/pi],...
          [approach.WPTTable(IF,2)*180/pi, approach.WPTTable(FAF,2)*180/pi],'LineWidth',2,'Color','black');
     geoshow([approach.WPTTable(FAF,1)*180/pi, approach.WPTTable(MAWP,1)*180/pi],...
          [approach.WPTTable(FAF,2)*180/pi, approach.WPTTable(MAWP,2)*180/pi],'LineWidth',2,'Color','black');
     geoshow([approach.WPTTable(MAWP,1)*180/pi, approach.WPTTable(MAHP,1)*180/pi],...
          [approach.WPTTable(MAWP,2)*180/pi, approach.WPTTable(MAHP,2)*180/pi],'LineWidth',2,'LineStyle','--','Color','black');

    
end
