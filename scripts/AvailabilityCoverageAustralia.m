
clear all
feature accel on
GPSConstants;


%--------------------------
%Optimised GPS Constellation - DO229C WAAS MOPS
%---------------------------

%dont forget to run GPSConstants above before running this code by itself.
FilenameOptimised = 'GPSConstellationOptimised.txt';

%read optimised constellation data from file
[NavData, TruePosVec] = freadOptimisedNav(FilenameOptimised);

%starting point DECEMBER 1, 1995 AT 00:00 UTC
GPSWeek = 829;
GPSSec = 432000;


TimeInterval = 3600; %seconds %for each hour.

%  set these to size 32 , otherwise if the highest satellite is only 21 (22
%drops out) then the size of this array is only 21 and it conflicts with
%previous sizes of 22, so set it fixed to 32
SatPos(1:4,1:32,1:4) = zeros(4,32,4);
SatVel(1:4,1:32,1:4) = zeros(4,32,4);
%this is the GPS processor which calculates simulated pseudoranges and
%satellite positions for the  GPSandModelEKF.m


%--------------------------
%Optimised GPS Constellation - DO229C WAAS MOPS
%---------------------------


for i = 1:8  %8 hour periods


    n = 1;

    Sat_PRN_Vec = zeros(1,32);  %initialise the vector with zeros. Columns of this vector is the PRN number of the satellite

    for SV = 1:32

        [SV_X_Data(SV) SV_Y_Data(SV) SV_Z_Data(SV) SV_T_Data(SV) ValidData_Satellite(SV)] = GPSOrbitPropagatorOptimal(GPSWeek, GPSSec, SV, NavData);

        [SV_Xvel_Data(SV) SV_Yvel_Data(SV) SV_Zvel_Data(SV) SV_Tvel_Data(SV) SV_Xacc_Data(SV) SV_Yacc_Data(SV) SV_Zacc_Data(SV) SV_Tacc_Data(SV) ValidDataSatVels(SV)] = GPSOrbitPropagatorOptimalVelocities(GPSWeek, GPSSec, SV, NavData);


        if  (ValidData_Satellite(SV) == 1 && ValidDataSatVels(SV) == 1)
            % if (ElevationBody(SV) > ElevationMask)  && (ValidData_Satellite(SV) == 1 && ValidDataSatVels(SV) == 1)
           
            SatPos(i,SV,1) = SV_X_Data(SV);
            SatPos(i,SV,2) = SV_Y_Data(SV);
            SatPos(i,SV,3) = SV_Z_Data(SV);
            SatPos(i,SV,4) = c*SV_T_Data(SV); %convert to metres

            SatVel(i,SV,1) = SV_Xvel_Data(SV);
            SatVel(i,SV,2) = SV_Yvel_Data(SV);
            SatVel(i,SV,3) = SV_Zvel_Data(SV);
            SatVel(i,SV,4) = c*SV_Tvel_Data(SV); %convert to metres

            Sat_PRN_Vec(i,SV) = 1;

        end
    end


    GPSSec = GPSSec+TimeInterval;

    %check for GPS week rollover

    if GPSSec > 604800
        GPSWeek = GPSWeek+1;
        GPSSec = 0;
    end

end



 Resolution = 1;  % degrees of latitude or longitude

    % define analysis area as australian FIR
    LatitudeMax = -27;
    LatitudeMin = -35;
    LongitudeMax = 140;
    LongitudeMin = 130;
    Height = 10000/3; %metres MSL  %3048 is 10,000 ft

    
    %put if statments in here to stop calculating positions for those lats
    %and lons which are in hte ocean below australia and other unwanted
    %regions.

    %Block out the following Regions
    
 

    % work out latitude and longitude bins
%     LongitudeBins = (LongitudeMax - LongitudeMin) / Resolution;
%     LatitudeBins = (LatitudeMax - LatitudeMin) / Resolution;

    %MaxIterations = floor(LongitudeBins * LatitudeBins);

    
    %want to block out this region, - southern region of australia.
%     LatitudeMinBlock = -40;
%     LatitudeMaxBlock = -36;
%     LongitudeMaxBlock = 135;
%     LongitudeMinBlock = 120;
           
    LatitudeMinBlock = 0;
    LatitudeMaxBlock = 0;
    LongitudeMaxBlock =0;
    LongitudeMinBlock = 0;
           
    %these contain the lat lon coordinates required
        
    Longitude = LongitudeMin;
    p = 1;
    canincrement = 0;
    while Longitude < LongitudeMax
        Latitude = LatitudeMin;
        q = 1;
        while Latitude < LatitudeMax

            %if (Longitude > LongitudeMaxBlock || Longitude < LongitudeMinBlock) && (Latitude > LatitudeMaxBlock || Latitude < LatitudeMinBlock)
                
             if (Longitude < LongitudeMaxBlock && Longitude > LongitudeMinBlock) && (Latitude < LatitudeMaxBlock && Latitude > LatitudeMinBlock)
                 
                             
               %don't want the ones in here, so skip them. 
               Latitude = Latitude + Resolution;
               %Longitude = Longitude + Resolution;       
              % skipp = 1;
              
               canincrement = 0;
             else
            
              canincrement = 1;
              LLHPositions(p,q,1) = Latitude;        
              LLHPositions(p,q,2) = Longitude;  
              LLHPositions(p,q,3) = Height;                 
              q = q+1;  
              Latitude = Latitude + Resolution;
              skipp = 0;
             end
        end       
    if canincrement == 1
         p = p+1;
    end     
         
         Longitude = Longitude + Resolution;
    end

    

    
    sizeposvec = size(LLHPositions);
    
    
    len = sizeposvec(1);
    wid = sizeposvec(2);
    
    
    ui = 1;
    for po = 1:len
        
        for wo = 1:wid
            
             LatitudePos(i,ui) = LLHPositions(po,wo,1)
              LongitudePos(i,ui) = LLHPositions(po,wo,1)
               HeightPos(i,ui) = LLHPositions(po,wo,1)
            
               
               
        [Position] = LLH2ECEF(LatitudePos(i,ui),LongitudePos(i,ui),HeightPos(i,ui));

        UserPosECEF(i,ui,1) = Position(1);
        UserPosECEF(i,ui,2) = Position(2);
        UserPosECEF(i,ui,3) = Position(3);
               
               
            ui = ui+1;
            
        end
    end
    
            



for i = 1:8  %this is the time

    for ui = 1:20

        n = 1;
        for SV = 1:32

            %calculate satellite elevations

            [AzimuthTemp(SV), ElevationTemp(SV)] = AzEl([UserPosECEF(i,ui,1),UserPosECEF(i,ui,2),UserPosECEF(i,ui,3)], [SatPos(i,SV,1), SatPos(i,SV,2), SatPos(i,SV,3)]);
            Azimuth(i,ui,SV) = rad2deg(AzimuthTemp(SV));
            Elevation(i,ui,SV) = rad2deg(ElevationTemp(SV));

            %
            %         [AzimuthBody(SV), ElevationBody(SV)] = AzElAntenna(UserPosECEF(1:3), [SV_X_Data(SV), SV_Y_Data(SV), SV_Z_Data(SV)],PHI,THETA,PSI);
            %
            %         AzimuthBody(SV) = rad2deg(AzimuthBody(SV));
            %         ElevationBody(SV)  = rad2deg(ElevationBody(SV));

            ElevationMask = 0;

            if (Elevation(i,ui,SV) > ElevationMask) && (ValidData_Satellite(SV) == 1 && ValidDataSatVels(SV) == 1)

                SV_AboveElevationMask(i,ui,n) = SV ;  %PRN numbers of the satellites above the elevation mask

                Sat_PRN_Vec(SV) = 1;

                n = n+1;

            end

        end %for SV = 1:32
                
            %number of satellites above elevation mask
            N(i,ui) = n-1;

    end %ui

end %i




for i = 1:8
    for ui = 1:20
        for k = 1:N(i,ui)


            SV =  SV_AboveElevationMask(i,ui,k);

            SVPos =  [SatPos(i,SV,1), SatPos(i,SV,2), SatPos(i,SV,3)];


            for m = 1:3
                ele(m) =  SVPos(m) - UserPosECEF(i,ui,m);
            end

            r_VecCalc(k) =  norm(ele);

            %Calculate DOPS

            M(k,1) =  -(SVPos(1) - UserPosECEF(i,ui,1))/r_VecCalc(k);
            M(k,2) =  -(SVPos(2) - UserPosECEF(i,ui,2))/r_VecCalc(k);
            M(k,3) =  -(SVPos(3) - UserPosECEF(i,ui,3))/r_VecCalc(k);
            M(k,4) = 1;

        end %for k

        Telev2 = T_ECEF2LTP(LongitudePos(i,ui),LatitudePos(i,ui));

        %make Telev2 into a 4 by 4

        Telev2(4,1:3) = [0 0 0];
        Telev2(1:4,4) = [0 0 0 1];

        H_LTP = M*Telev2';

        AA = (H_LTP'*H_LTP)^-1;

        var_x(i,ui) = AA(1,1);
        var_y(i,ui) = AA(2,2);
        var_z(i,ui) = AA(3,3);
        var_dt(i,ui) = AA(4,4);

        GDOP(i,ui) = sqrt(var_x(i,ui) + var_y(i,ui) + var_z(i,ui) + var_dt(i,ui));

        PDOP(i,ui) = sqrt(var_x(i,ui) + var_y(i,ui) + var_z(i,ui));

        HDOP(i,ui) = sqrt(var_x(i,ui) + var_y(i,ui));

        VDOP(i,ui) = sqrt(var_z(i,ui));

        TDOP(i,ui) = sqrt(var_dt(i,ui));



    end %ui


end %i



%Multiply by expected pseudorange noise to get position error



%Display results

%PLOT MAP OF COASTLINE


figure();
if ~exist('austcoast')
    load 'data/austcoast.dat';
end
plot(austcoast(260000:10:417056,1),austcoast(260000:10:417056,2),'.');


%surf(austcoast(260000:10:417056,1),austcoast(260000:10:417056,2),33*ones(size(austcoast(260000:10:417056,2))));


grid on;
xlabel('Longitude (deg)');
ylabel('Latitude (deg)');

hold;



%PLOT USER LOCATION



plot(LongitudePos(:,:),LatitudePos(:,:),'*r');





surf(LongitudePos(i,ui),LatitudePos(i,ui),HeightPos(i,ui));


plot3(LatitudePos(:,:),LongitudePos(:,:),HeightPos(:,:))


%PLOT GRAS ANTENNA LOCATION

%overlay color code for HDOP



contour3(LongitudePos(:,1:5),LatitudePos(:,1:5),HDOP(:,1:5));


contour(HDOP)


pcolor(LongitudePos(:,1:5),LatitudePos(:,1:5),HDOP(:,1:5));



pcolor(LLHPositions(:,:,2)',LLHPositions(:,:,1)',HPL');

    %shading facet;  % set shading to interpreted
grid on;
caxis([0 500]);  % set the color axis
colorbar;  % show a color bar
hold;

 plot(austcoast(:,1),austcoast(:,2),'k', 'LineWidth',2);
    grid on;
    xlabel('Longitude (deg)');
    ylabel('Latitude (deg)');        
    
shading interp;





















