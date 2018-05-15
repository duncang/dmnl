function [g] = Earth_Gravity(P,e2,ye,k)
s2 = sin(P(1))^2;
g = ye*(1+k*s2)/((1-e2*s2)^(1/2))- ...
    (3.0877e-6-0.0044e-6*s2)*P(3)+0.072e-12*P(3)*P(3);