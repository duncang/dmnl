function output = Convert_NED2ECEF(input);
%Troy Bruggemann 27 Feb 06 not verified yet
%this functino is for simulink.



Latitude = input(1);
Longitude = input(2);
North = input(3);
East = input(4);
Down = input(5);



NED = [North East Down]';



%split vector




TMatrix = T_ECEF2NED(Latitude,Longitude);

Tned2ecef = TMatrix';

ECEF = Tned2ecef*NED;



%output vector

output = ECEF;

