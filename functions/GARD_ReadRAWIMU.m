function rawimu = ReadRAMIMU(file)
% rawimu = ReadRAMIMU(file)
%
%
%

% file = 'log2/log_span_25jun09_ybaf_2.ASC';

% %RAWIMUSA,1537,361705.645;1537,361705.644826900,00000000,30387,415,-5398,
% -23039,-18100736,-5117*c01b540c
f = fopen(file,'r');

%preallocate rawimu
rawimu.log_week = zeros(1,130000);
rawimu.log_sec = zeros(1,130000);
rawimu.obs_week = zeros(1,130000);
rawimu.obs_sec = zeros(1,130000);
rawimu.status = zeros(1,130000);
rawimu.XAccel = zeros(1,130000);
rawimu.YAccel = zeros(1,130000);
rawimu.ZAccel = zeros(1,130000);
rawimu.XRate = zeros(1,130000);
rawimu.YRate = zeros(1,130000);
rawimu.ZRate = zeros(1,130000);

gyro_sf = 0.1e-8/3600; % deg/LSB
accel_sf = 0.05e-15;   % m/sec/LSB

count = 0;
while ~feof(f)
    line = fgets(f);
    
    if length(line) < 9
        continue;
    end
    
    %% short header
    if strcmp(line(1:9), '%RAWIMUSA')
        count = count+1;
        
        if mod(count,100) == 0
            disp(count);
        end
        
        %data = textscan(line,'%s,%n,%n;%n,%n,%n,%n,%n,%n,%n,%n,%n%s');
        
        % get the start
        [tok, rem] = strtok(line,',');
        % tok should equal '%RAWIMUSA'
        
        
        % get the week
        [tok, rem] = strtok(rem,',');
        rawimu.log_week(count) = str2num(tok);
        
        % get the second
        [tok, rem] = strtok(rem,';');
        rawimu.log_sec(count) = str2num(tok);
        
        [tok, rem] = strtok(rem,',');
        rawimu.obs_week(count) = str2num(tok);
        
        [tok, rem] = strtok(rem,',');
        rawimu.obs_sec(count) = str2num(tok);
        
        [tok, rem] = strtok(rem,',');
        rawimu.status(count) = str2num(tok);
        
        [tok, rem] = strtok(rem,',');
        rawimu.ZAccel(count) =  str2num(tok) * accel_sf;
        
        [tok, rem] = strtok(rem,',');
        rawimu.YAccel(count) =  str2num(tok) * -accel_sf;
        
        [tok, rem] = strtok(rem,',');
        rawimu.XAccel(count) =  str2num(tok) * accel_sf;
        
        [tok, rem] = strtok(rem,',');
        rawimu.ZRate(count) =  str2num(tok) * gyro_sf;
        
        [tok, rem] = strtok(rem,',');
        rawimu.YRate(count) =  str2num(tok) * -gyro_sf;
        
        [tok, rem] = strtok(rem,'*');
        rawimu.XRate(count) =  str2num(tok) * gyro_sf;
    end
    
    %% long header
    if strcmp(line(1:8), '#RAWIMUA')
        count = count+1;
        
        if mod(count,100) == 0
            disp(count);
        end
        
        %data = textscan(line,'%s,%n,%n;%n,%n,%n,%n,%n,%n,%n,%n,%n%s');
        
        % get the start
        [tok, rem] = strtok(line,',');
        % tok should equal '%RAWIMUSA'
        
        % bit of garbage at the start
        [tok, rem] = strtok(rem,','); %port
        [tok, rem] = strtok(rem,','); %??
        [tok, rem] = strtok(rem,','); %cpu free
        [tok, rem] = strtok(rem,','); %clock state
        
        % get the week
        [tok, rem] = strtok(rem,',');
        rawimu.log_week(count) = str2num(tok);
        
        % get the second
        [tok, rem] = strtok(rem,',');
        rawimu.log_sec(count) = str2num(tok);
        
        % bit more garbage
        [tok, rem] = strtok(rem,','); % imu status
        [tok, rem] = strtok(rem,','); % ??
        [tok, rem] = strtok(rem,';'); % ??
        
        
        [tok, rem] = strtok(rem,',');
        rawimu.obs_week(count) = str2num(tok);
        
        [tok, rem] = strtok(rem,',');
        rawimu.obs_sec(count) = str2num(tok);
        
        [tok, rem] = strtok(rem,',');
        rawimu.status(count) = str2num(tok);
        
        [tok, rem] = strtok(rem,',');
        rawimu.ZAccel(count) =  str2num(tok) * accel_sf;
        
        [tok, rem] = strtok(rem,',');
        rawimu.YAccel(count) =  str2num(tok) * -accel_sf;
        
        [tok, rem] = strtok(rem,',');
        rawimu.XAccel(count) =  str2num(tok) * accel_sf;
        
        [tok, rem] = strtok(rem,',');
        rawimu.ZRate(count) =  str2num(tok) * gyro_sf;
        
        [tok, rem] = strtok(rem,',');
        rawimu.YRate(count) =  str2num(tok) * -gyro_sf;
        
        [tok, rem] = strtok(rem,'*');
        rawimu.XRate(count) =  str2num(tok) * gyro_sf;
    end
end 


fclose(f);

