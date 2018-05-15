function TMatrix = T_Body2NED(PHI,THETA, PSI)
% generate  body frame to Fixed frame (NED) transformation matrix based on Euler
% angles. 

%by Troy Bruggemann 14 June. Taken from p 102 Nelson.  
%this has been verified with other sources to be correct (seeing Nelson has many mistakes in it)
%note: the transpose of this matrix can be used for NED2Body
%transformation.



Cpsi = cos(PSI);
Ctheta = cos(THETA);
Cphi =  cos(PHI);

Spsi = sin(PSI);
Stheta = sin(THETA);
Sphi = sin(PHI);





TMatrix(1,1) = Ctheta*Cpsi;
TMatrix(1,2) =  Sphi*Stheta*Cpsi - Cphi*Spsi;
TMatrix(1,3) =  Cphi*Stheta*Cpsi + Sphi*Spsi;

TMatrix(2,1) = Ctheta*Spsi;
TMatrix(2,2) = Sphi*Stheta*Spsi + Cphi*Cpsi;
TMatrix(2,3) = Cphi*Stheta*Spsi - Sphi*Cpsi;

TMatrix(3,1) = -Stheta;
TMatrix(3,2) = Sphi*Ctheta;
TMatrix(3,3) = Cphi*Ctheta;


