function success = ExportGoogleEarth_KMLPathPoints(filename,name,description,altitudeMode,data,linecolor,polygoncolor,tessellate,extrude)
% Exports data as a path to be displayed in GoogleEarth. A path is either
% be displayed as a 2D path on the ground or a 3D path with your altitude
% data. (Note: Altitude is not mandatory)
% 
% More information availbale on:
% http://code.google.com/apis/kml/documentation/kml_tags_21.html
% http://code.google.com/apis/kml/documentation/kml_tut.html#paths
% 
% 
% input: 
%   filename:     name of .kml file
%   name:         Name of the Project-Data
%   description:  Description of the Project
%   altitudeMode: 'clampToGround' / 'relativeToGround' / 'absolute'
%   data:         A 2 or 3 coloumn matrix with the data in [DegLongitude DegLatitude] or [DegLongitude DegLatitude Altitude]
%   linecolor:    aabbggrr: (default=ffffffff)
%   polygoncolor: aabbggrr: (default=7f99ff99)
%   tessellate:   boolean: if lines follows terrain (default=1)
%   extrude:      boolean: if line is extruded to the ground (default=1)
%
% GoogleEarth can easely handle 10'000+ datapoints, but I couldn't make it
% displaying 100'000. I guess there is a limit somewhere. 
% 
% Example usage:
% ExportGoogleEarth_KMLPath(...
%     'flight.kml',...
%     'Flight Alex Gurtner 20070614',...
%     'This is the flight with the fisheye lens. From Archerfield to Cleveland, Redland Bay, Goldcoast, Logan, Archerfield',...
%     'absolute',...
%     [Pos_LLH(:,2)*180/pi, Pos_LLH(:,1)*180/pi, round(Pos_LLH(:,3))]...
%     )
% 
% Initial: Alex Gurtner, 12/10/2007 v1.0
% Mods: Alex, 12/10/2007 added support for no altitude data v1.1

%Modified by Troy Bruggemann 10.11.09 
%converts csv and placename data for gps pole coordinates to kml.


display(['Exporting data to ' filename '...']);

success = 0;

try
    
    % setting not yet set values
%     if nargin < 6
%         linecolor = 'ffffffff';
%         display('setting linecolor (default=ffffffff)');
%     end
%     if nargin < 7
%         polygoncolor = '7f99ff99';
%         display('setting polygoncolor (default=7f99ff99)');
%     end
%     if nargin < 8
%         tessellate = 1;
%         display('setting tessellate (default=1)');
%     end
%     if nargin < 9
%         extrude = 1;
%         display('setting extrude (default=1)');
%     end
%     
%     
%     % in case of no altitude data, we draw the line to the ground
%     if size(data,2) < 3
%         display('MESSAGE: No altitude given, setting altitudeMode to ''clampToGround''');
%         altitudeMode = 'clampToGround';
%     end
    
    
    % Creating KML-File
    [fid,message] = fopen(filename,'w+');
    if fid == -1
        display('ERROR: Failed opening file');
        display(['System message: ' message]);
        return;
    end
    
    
%     fprintf(fid, '<?xml version="1.0" encoding="UTF-8"?>\n');
%     fprintf(fid, '<kml xmlns="http://earth.google.com/kml/2.2">\n');
%     fprintf(fid, '  <Document>\n');
%     fprintf(fid, '    <name>%s</name>\n',name);
%     fprintf(fid, '    <description>%s</description>\n',description);
%     fprintf(fid, '    <Style id="TheStyle">\n');
%     fprintf(fid, '      <LineStyle>\n');
%     fprintf(fid, '        <color>%s</color>\n',linecolor);
%     fprintf(fid, '        <width>2</width>\n');
%     fprintf(fid, '      </LineStyle>\n');
%     fprintf(fid, '      <PolyStyle>\n');
%     fprintf(fid, '        <color>%s</color>\n',polygoncolor);
%     fprintf(fid, '      </PolyStyle>\n');
%     fprintf(fid, '    </Style>\n');
%     fprintf(fid, '    <Placemark>\n');
%     fprintf(fid, '      <name>Absolute Extruded</name>\n');
%     fprintf(fid, '      <description>none</description>\n');
%     fprintf(fid, '      <styleUrl>#TheStyle</styleUrl>\n');
%     fprintf(fid, '      <LineString>\n');
%     fprintf(fid, '        <extrude>%d</extrude>\n',extrude);
%     fprintf(fid, '        <tessellate>%d</tessellate>\n',tessellate);
%     fprintf(fid, '        <altitudeMode>%s</altitudeMode>\n',altitudeMode);
%     fprintf(fid, '        <coordinates> ');
%     
    
    % Writing data to KML-file
    % -112.2550785337791,36.07954952145647,2357
    if size(data,1) > 20000
        disp('WARNING: GoogleEarth may fail to display so many path-points')
    end
    
%     if size(data,2) < 3
%         disp('MESSAGE: Could not find altitude data')
%         for i=1:length(data(:,1))
%             fprintf(fid, '%f,%f\n', data(i,1), data(i,2) );
%         end
%     else
%         disp('MESSAGE: Found altitude data' )
%         for i=1:length(data(:,1))
%             fprintf(fid, '%f,%f,%f\n', data(i,1), data(i,2), round(data(i,3)) );
%         end
%     end
    
    
    % Closing coords tags
%     fprintf(fid, '        </coordinates>\n');
%     fprintf(fid, '      </LineString>\n');
%     fprintf(fid, '    </Placemark>\n');
%     
%     
%     
%     % add waypoints    
%     fprintf(fid, '     <Folder>\n');
%     fprintf(fid, '      <name>Waypoints</name>\n');
    
   % inc_percent = 5;
    
   % for iwp=1:floor(length(data(:,1))*inc_percent/100):length(data(:,1))   
      
          
   
   
 fprintf(fid, '<?xml version="1.0" encoding="utf-8" ?>\n');
 fprintf(fid, '   <kml xmlns="http://earth.google.com/kml/2.1">\n');
 fprintf(fid, '   <Folder>\n');
 fprintf(fid, '    <description>\n');
 fprintf(fid, '    <![CDATA[ GPS Pole Sites  ]]>\n');
 fprintf(fid, '     </description>\n');
   
   
    for iwp=1:length(data(:,1))       
        
        fprintf(fid, '<Placemark>\n');
        fprintf(fid, '   <name>\n');
         fprintf(fid, '    <![CDATA[%d ]]>\n', data(iwp,3));
         fprintf(fid, '  </name>\n');
        fprintf(fid, '  <description>\n');
        fprintf(fid, '  <![CDATA[ <br><br><br>\n <table border="1" padding="0">\n');
        fprintf(fid, '      <tr><td>ID</td><td>%d</td></tr>\n',data(iwp,3));
        fprintf(fid, '      <tr><td>Lat</td><td>%f</td></tr>\n', data(iwp,2));
        fprintf(fid, '      <tr><td>Long</td><td>%f</td></tr>\n', data(iwp,1));       
        fprintf(fid, '              ]]>\n'); 
        fprintf(fid, '</description>\n');
        
        fprintf(fid, '<visibility>1</visibility>\n'); 
         fprintf(fid, '<open>0</open>\n');  
         fprintf(fid, '<Point>\n');          
         fprintf(fid, '  <extrude>1</extrude>\n');
          fprintf(fid, '  <altitudeMode>relativeToGround</altitudeMode>\n');            
          fprintf(fid, '  <coordinates>%f,%f,0</coordinates>\n',data(iwp,1),data(iwp,2));
          fprintf(fid, '</Point>\n');
        fprintf(fid, '</Placemark>\n');

          
          
        
%         fprintf(fid, '      <styleUrl>root://styles#default</styleUrl>\n');
%         fprintf(fid, '      <Style>\n');
%         fprintf(fid, '        <IconStyle>\n');
%         fprintf(fid, '          <Icon>\n');
%         fprintf(fid, '            <href>root://icons/palette-4.png</href>\n');
%         fprintf(fid, '            <y>128</y>\n');
%         fprintf(fid, '            <w>32</w>\n');
%         fprintf(fid, '            <h>32</h>\n');
%         fprintf(fid, '          </Icon>\n');
%         fprintf(fid, '        </IconStyle>\n');
%         fprintf(fid, '      </Style>\n');
%         fprintf(fid, '      <Point>\n');
%         fprintf(fid, '        <extrude>%d</extrude>\n', extrude);
%         fprintf(fid, '        <altitudeMode>%s</altitudeMode>\n',altitudeMode);
%         fprintf(fid, '        <coordinates>\n');
        
        %if size(data,2) < 3
        %    fprintf(fid, '%f,%f\n', data(iwp,1), data(iwp,2) );
        %else
         %   fprintf(fid, '%f,%f,%f\n', data(iwp,1), data(iwp,2), round(data(iwp,3)) );
        %end
        
        %fprintf(fid, '       </coordinates>\n');
        %fprintf(fid, '      </Point>\n');
        %fprintf(fid, '    </Placemark>\n');
        
    end
    
    
%     
%     <Placemark>
% - <name>
% - <![CDATA[ 910204 ]]> 
%   </name>
% - <description>
% - <![CDATA[ <br><br><br>
%     <table border="1" padding="0">
%     <tr><td>ID</td><td>910204</td></tr>
%     <tr><td>Lat</td><td>-27.6731</td></tr>
%     <tr><td>Long</td><td>147.945</td></tr>
%         
%   ]]> 
%   </description>
%   <visibility>1</visibility> 
%   <open>0</open> 
% - <Point>
%   <extrude>1</extrude> 
%   <altitudeMode>relativeToGround</altitudeMode> 
%   <coordinates>147.945197,-27.673085,0</coordinates> 
%   </Point>
%   </Placemark>

    
    
    
    
    
    
    
    
    
    
    fprintf(fid, ' </Folder>\n');
    
    
    
    
    
    
    
    
    
    
    
    
    
    % closing KML file
    %fprintf(fid, '  </Document>\n');
    fprintf(fid, '</kml>\n');
    
    fclose(fid);

    success = 1;
    display(['... export completed']);
    
catch
    
    display(['ERROR: Export was not succesfull']);
    
    % try to close file, if it was opened
    if exist('fid')
        fclose(fid);
    end
    
    % no success
    success = 0;
    
end

