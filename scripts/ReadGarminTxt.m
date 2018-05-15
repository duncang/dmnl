%%% this function reads Garmin Text-out data stream
%% Ref: http://www.gpsinformation.org/dale/interface.htm
%% Written By Duncan Greer 4 June 2006
%% $Id: ReadGarminTxt.m 1884 2008-07-15 05:54:33Z n2523710 $
%%

Filename = 'c:\garmin_4jun06_1501.txt';
DataFile = fopen(Filename);

% read line by line until EOF
Epoch=0;
while(~feof(DataFile))
    
    DataLine = fgetl(DataFile);
    Epoch = Epoch + 1;
    
    %%%% Valid Data
    
%     'd' if current 2D differential GPS position
%     'D' if current 3D differential GPS position
%     'g' if current 2D GPS position
%     'G' if current 3D GPS position
%     'S' if simulated position
%     '_' if invalid position

    if(DataLine(31) == '_')
        ValidData(Epoch) = -1;
        break;
    else
        ValidData(Epoch) = 1;
    end
    
    %%%% TIME 
    
    Hour = str2num(DataLine(8:9));
    Minute = str2num(DataLine(10:11));
    Second = str2num(DataLine(12:13));
    
    Time(Epoch) = Hour * 3600 + Minute * 60 + Second;
    
    %%% POSITION
    
    if(DataLine(14) == 'N')
        LatSign = 1;
    else
        LatSign = -1;
    end
    
    if(DataLine(22) == 'E')
        LongSign = 1;
    else
        LongSign = -1;
    end
    
    if(DataLine(35) == '+')
        AltSign = 1;
    else
        AltSign = -1;
    end
    
    Latitude(Epoch) = LatSign*(str2num(DataLine(15:16)) + str2num(DataLine(17:21))/(60000));
    Longitude(Epoch) = LongSign*(str2num(DataLine(23:25)) + str2num(DataLine(26:30))/(60000));
    Altitude(Epoch) = AltSign*(str2num(DataLine(36:40)));
    
    % EPH - estimated horizontal position error
    EPH(Epoch) = str2num(DataLine(32:34));
    
    %%%%% Velocity data
    if(DataLine(41) == 'E')
        VelESign = 1;
    else
        VelESign = -1;
    end
    
    if(DataLine(46) == 'N')
        VelNSign = 1;
    else
        VelNSign = -1;
    end
    
    if(DataLine(51) == 'U')
        VelUSign = 1;
    else
        VelUSign = -1;
    end
    
    VelE(Epoch) = VelESign*(str2num(DataLine(42:45))/10);
    VelN(Epoch) = VelNSign*(str2num(DataLine(47:50))/10);
    VelU(Epoch) = VelUSign*(str2num(DataLine(52:55))/10);
    
    
    
end

