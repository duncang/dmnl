

DataPath = 'data/Flight_Data/quas_ins_yred-ykry/';
FilePrefix = 'log_20100804-quas_ins_';
NavDataPath = strcat(DataPath,FilePrefix,'nav.out');
GPSDataPath = strcat(DataPath,FilePrefix,'gps.out');
IMUDataPath = strcat(DataPath,FilePrefix,'imu.out');
SPANDataPath = strcat(DataPath,FilePrefix,'span.out');


imudata = dlmread(IMUDataPath);
navdata = dlmread(NavDataPath);
spandata = dlmread(SPANDataPath);


