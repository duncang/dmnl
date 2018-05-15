function nav = GARD_LoadNav(filename)


fid = fopen(filename);

linedata = fgetl(fid);
linedata = fgetl(fid);

ndbindex = 0;
vorindex = 0;

while ~feof(fid)
   linedata = fgetl(fid);
   
   if str2num(linedata(1:2)) == 99
       break;
   end
   
   linelength = length(linedata);
   
   NavType = str2num(linedata(1:2));
   NavLat = str2num(linedata(4:13));
   NavLon = str2num(linedata(15:25));
   NavEle = str2num(linedata(27:32));
   
   if(NavLat > -10 || NavLat < -50)
        continue;
    end
    
    if(NavLon < 110 || NavLon > 160)
        continue;
    end
   
   switch NavType
       case 1

       case 2
           % ndb
           ndbindex = ndbindex + 1;
           NavFreq = str2num(linedata(36:38));
           NavCode = linedata(52:55);
           NavName = linedata(57:linelength);
           
           nav.NDB(ndbindex).data(1:4) = double(NavCode);
           nav.NDB(ndbindex).data(5) = NavLat*pi/180;
           nav.NDB(ndbindex).data(6) = NavLon*pi/180;
           nav.NDB(ndbindex).data(7) = NavEle;
           nav.NDB(ndbindex).name = NavName;
           nav.NDB(ndbindex).frequency = NavFreq;
           
       case 3
           % vor
           vorindex = vorindex+1; 
           NavFreq = str2num(linedata(34:38))/100;
           NavCode = linedata(52:55);
           NavName = linedata(57:linelength);
           nav.VOR(vorindex).data(1:4) = double(NavCode);
           nav.VOR(vorindex).data(5) = NavLat*pi/180;
           nav.VOR(vorindex).data(6) = NavLon*pi/180;
           nav.VOR(vorindex).data(7) = NavEle;
           nav.VOR(vorindex).name = NavName;
           nav.VOR(vorindex).frequency = NavFreq;
       case 4
           % ILS

       case 5
            % localizer
       case 6
            % Glide slope
       case 7
            % outer marker
       case 8
            % middle marker
       case 9
            % inner marker
       case 10

       case 11

       case 12
           %  	DME (including the DME element of an ILS, VORTAC or
           %  	VOR-DME). 

       case 13
           %  	DME (including the DME element of an ILS, VORTAC or
           %  	VOR-DME).
       
   end
   
end

fclose(fid);