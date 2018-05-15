function DeltaTropo = GARD_TropoDelay(Elevation,UserHeight)
% function DeltaTropo = GARD_TropoDelay(Elevation,UserHeight)
% Implements a simplified tropospheric delay model from
% the Blue Book, Page 541
% UserHeight given in metres
% Elevation given in radians

DeltaTropo = 2.47 * exp(-0.133*UserHeight/1000) / (sin(Elevation) + 0.0121);

