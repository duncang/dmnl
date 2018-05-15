function [DOP] = GARD_CalculateDOPS2(UserPos,SVPos)
%DOPS
%put the DOPS in a vector, the order is [GDOP PDOP HDOP VDOP TDOP]
%
% $Id: GARD_LSQ.m 895 2007-11-09 06:25:40Z greerd $
%

Xb = UserPos(1);   
Yb = UserPos(2);
Zb = UserPos(3);
dtb = UserPos(4);


TempBase_B_Pos = [Xb,Yb,Zb, dtb];

N = size(SVPos,1);

%Calculated slant ranges

for k = 1:N

    for m = 1:3
        ele(m) =  SVPos(k,m) - TempBase_B_Pos(m);
    end

    r_VecCalc(k) =  norm(ele);

end




%Design Matrix
%generate elements of M matrix

for k = 1:N
    M(k,1) =  -(SVPos(k,1) - TempBase_B_Pos(1))/r_VecCalc(k);
    M(k,2) =  -(SVPos(k,2) - TempBase_B_Pos(2))/r_VecCalc(k);
    M(k,3) =  -(SVPos(k,3) - TempBase_B_Pos(3))/r_VecCalc(k);
    M(k,4) = 1;
end
    

%Transform M matrix into LTP coordinates (for HDOP/VDOP calculation)
        
Position = [Xb, Yb, Zb];


[Latitude,Longitude,Height] = ECEF2LLH(Position);



Telev2 = T_ECEF2LTP(Longitude,Latitude);

%make Telev2 into a 4 by 4

Telev2(4,1:3) = [0 0 0];
Telev2(1:4,4) = [0 0 0 1];
           
    
    
H_LTP = M*Telev2';   



AA = (H_LTP'*H_LTP)^-1;

var_x = AA(1,1);
var_y = AA(2,2);
var_z = AA(3,3);
var_dt = AA(4,4);

GDOP = sqrt(var_x + var_y + var_z + var_dt);

PDOP = sqrt(var_x + var_y + var_z);

HDOP = sqrt(var_x + var_y);

VDOP = sqrt(var_z);

TDOP = sqrt(var_dt);

DOP = [GDOP PDOP HDOP VDOP TDOP];