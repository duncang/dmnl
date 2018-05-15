function TMatrix = T_pqrbody2EulerRate(PHI,THETA)
% generates angular velocities in body frame (pqr) to euler rates (phi
% dot, theta dot, psi dot
%transformation matrix based on euler angles
%NOTE: for pqr body to EulerRate the transpose of this matrix CANNOT be
%used. Use the T_EulerRate2pqrbody.m function

%by Troy Bruggemann 27 July. Taken from p 103 Nelson.  
%this has been verified to be correct


Ctheta = cos(THETA);
Cphi =  cos(PHI);


Stheta = sin(THETA);
Sphi = sin(PHI);



TMatrix(1,1) = 1;
TMatrix(1,2) =  Sphi*tan(THETA);
TMatrix(1,3) =  Cphi*tan(THETA);

TMatrix(2,1) = 0;
TMatrix(2,2) =Cphi;
TMatrix(2,3) = -Sphi;

TMatrix(3,1) = 0;
TMatrix(3,2) = Sphi*sec(THETA);
TMatrix(3,3) = Cphi*sec(THETA);
