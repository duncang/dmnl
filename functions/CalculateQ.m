function y = CalculateQ(x)
%Troy Bruggemann 2006
%This function yet to be verified! NOTE: Using the matlab chi-square
%functions instead of this one
%forms Q matrix used in "Integrated GPS/Inertial Fault Detection
%Availability, Mats Brenner

%Adapted from  y = ncchisq2(x) function from 
%"GPS RAIM:Calculation of thresholds and protection radius using chi-square methods - a geometric approach"  (R. Grover Brown).
%RTCA paper No.491-94/SC159-584

%N = number of terms in series minus one
%lambda is the noncentral parameter
N = 49;
lambda = 68.413;
prevterm = 1;
sum = 1;
for j = 1:N
    term = prevterm.*lambda.*x./(4.*j.*j);
    sum = sum+term;
    prevterm = term;
end
y = 0.5.*sum.*exp(-0.5.*(x+lambda))