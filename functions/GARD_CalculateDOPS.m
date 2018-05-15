function [GDOP PDOP HDOP VDOP TDOP] = GARD_CalculateDOPS(H)
% function [GDOP PDOP HDOP VDOP TDOP] = GARD_CalculateDOPS(H)
% Returns the DOP values for the given Nx4 H-matrix
% Written by Duncan Greer and Troy Bruggemann
% $Id: GARD_CalculateDOPS.m 1850 2008-07-14 04:52:47Z greerd $
%

AA = inv(H' * H);

var_x = AA(1,1);
var_y = AA(2,2);
var_z = AA(3,3);
var_dt = AA(4,4);

GDOP = sqrt(var_x + var_y + var_z + var_dt);
PDOP = sqrt(var_x + var_y + var_z);
HDOP = sqrt(var_x + var_y);
VDOP = sqrt(var_z);
TDOP = sqrt(var_dt);

