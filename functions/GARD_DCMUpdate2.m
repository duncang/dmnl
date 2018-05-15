function C_BN_new = GARD_DCMUpdate2(C_BN,Sigma)
% function C_BN_new = GARD_DCMUpdate2(C_BN,Sigma)
% C_BN is the current direction cosine matrix
% Sigma is an angle vector per Titterton 11.5
%
% References: Titterton Ch. 11
%

% form skew form of sigma (sigma_x)
Sigma_x = [0,-Sigma(3),Sigma(2); ...
           Sigma(3),0,-Sigma(1); ...
           -Sigma(2),Sigma(1),0];

%Sigma_x_dash = -Sigma_x';

Sigma_x_2 = [-(Sigma(2)^2 + Sigma(3)^2), Sigma(1)*Sigma(2), Sigma(1)*Sigma(3);  ...
             Sigma(1)*Sigma(2),-(Sigma(1)^2 + Sigma(3)^2),Sigma(2)*Sigma(3);  ...
             Sigma(1)*Sigma(3),Sigma(2)*Sigma(3),-(Sigma(1)^2 + Sigma(2)^2)];
       
% calculate sigma, the magnitude of the Sigma vector
s = sqrt(Sigma(1)^2 + Sigma(2)^2 + Sigma(3)^2);

% calculate the coefficients
a1 = 1 - s^2/factorial(3) + s^4/factorial(5);
a2 = 1/factorial(2) - s^2/factorial(4) + s^4/factorial(6);

% Direction Cosine matrix attitude update
Ak = eye(3) + a1 * Sigma_x + a2 * Sigma_x_2;

C_BN_new = C_BN * Ak;


C_BN_new = GARD_OrthogonaliseDCM(C_BN_new);