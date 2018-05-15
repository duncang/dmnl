function GARD_WriteSimGenMotionFile(OutputFileName,pos_truth_ecef,vel_truth,att_truth,sensors)


NumberSamples = length(pos_truth_ecef);
dt = 0.01;
vel_truth_ecef = zeros(NumberSamples,4);
acc_truth_ecef = zeros(NumberSamples,4);
jrk_truth_ecef = zeros(NumberSamples,4);

disp('Converting to ECEF');
for i = 1:NumberSamples
   C_BN = GARD_EulerToDCM(att_truth(2,i),att_truth(3,i),att_truth(4,i));
   [lat,long,hgt] = ECEF2LLH(pos_truth_ecef(2:4,i));
   C_NE = T_ECEF2NED(lat,long)';
   
   vel_truth_ecef(i,2:4) = (C_NE*vel_truth(2:4,i))';
   acc_truth_ecef(i,2:4) = (C_NE*(C_BN*(sensors(2:4,i)+[0;0;9.79])))';
   
end


pos_truth_ecef = pos_truth_ecef';
att_truth = att_truth';
sensors = sensors';

% first differentiate accel data to obtain jerk
disp('Calculating Jerks..');
jrk_truth_ecef(1:NumberSamples-1,2) = diff(acc_truth_ecef(:,2))/dt;
jrk_truth_ecef(1:NumberSamples-1,3) = diff(acc_truth_ecef(:,3))/dt;
jrk_truth_ecef(1:NumberSamples-1,4) = diff(acc_truth_ecef(:,4))/dt;


disp('Writing File');
Outfile = fopen(OutputFileName,'w');

    % write header
    i = 1;
    fprintf(Outfile,'%f,mot,v1_m1,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,0,0,0,0,0,0\n', ...
                pos_truth_ecef(i,1),...
                pos_truth_ecef(i,2),pos_truth_ecef(i,3),pos_truth_ecef(i,4),...
                vel_truth_ecef(i,2),vel_truth_ecef(i,3),vel_truth_ecef(i,4),...
                acc_truth_ecef(i,2),acc_truth_ecef(i,3),acc_truth_ecef(i,4),...
                jrk_truth_ecef(i,2),jrk_truth_ecef(i,3),jrk_truth_ecef(i,4),...
                att_truth(i,2),att_truth(i,3),att_truth(i,4),...
                sensors(i,5),sensors(i,6),sensors(i,7));
        
    fprintf(Outfile,'RU\n');
    
    for i=2:NumberSamples
            fprintf(Outfile,'%f,mot,v1_m1,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,0,0,0,0,0,0\n', ...
                pos_truth_ecef(i,1),...
                pos_truth_ecef(i,2),pos_truth_ecef(i,3),pos_truth_ecef(i,4),...
                vel_truth_ecef(i,2),vel_truth_ecef(i,3),vel_truth_ecef(i,4),...
                acc_truth_ecef(i,2),acc_truth_ecef(i,3),acc_truth_ecef(i,4),...
                jrk_truth_ecef(i,2),jrk_truth_ecef(i,3),jrk_truth_ecef(i,4),...
                att_truth(i,4),att_truth(i,3),att_truth(i,2),...
                sensors(i,5),sensors(i,6),sensors(i,7));
    end


fclose(Outfile);

