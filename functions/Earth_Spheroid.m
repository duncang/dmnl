function [M,N] = Earth_Spheroid(P,a,e2)
s2 = sin(P(1))^2;
M = (a*(1-e2))/((1-e2*s2)^(3/2));
N = a/((1-e2*s2)^(1/2));