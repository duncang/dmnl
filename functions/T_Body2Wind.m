function TMatrix = T_Body2Wind(alpha, beta)
% generate  body to wind transformation matrix  

%by Troy Bruggemann 22 july 2006.  





Calpha = cos(alpha);
Cbeta = cos(beta);
Salpha = sin(alpha);
Sbeta = sin(beta);




TMatrix(1,1) = Calpha*Cbeta;
TMatrix(1,2) =  Sbeta;
TMatrix(1,3) =  Salpha*Cbeta;

TMatrix(2,1) = -Calpha*Sbeta;
TMatrix(2,2) = Cbeta;
TMatrix(2,3) = -Salpha*Sbeta;

TMatrix(3,1) = -Salpha;
TMatrix(3,2) = 0;
TMatrix(3,3) = Calpha;


