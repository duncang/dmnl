%pinv example with svd

A = [9 12 15; 6 12 6; 3 3 3;]




[U, S, V] = svd(A)



%find inverse of A


Ainv = V*inv(S)*U'


AinvPinv = pinv(A)

Ainvinv = inv(A)
