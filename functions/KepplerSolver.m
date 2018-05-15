function E_final = KepplerSolver(M,e)

% Keppler Solver iteraly solves for the eccentric anomoly using iterative
% method described on page 164 of GPS theory and applciations volume 1

% written by Duncan Greer (c) CRCSS 2005

% M error threshold
Threshold = 1e-13;
% initial guess
E0 = M + e*sin(M) / (1 - sin(M+e) + sin(M));

% do iterations until the error is less than the thresshold
n = 1;
E(n) = E0;

while 1
    E(n+1) = E(n) - (E(n) - e*sin(E(n)) - M) / (1 - e*cos(E(n)));
    
    % break loop when M can be solved from E
    M_check = E(n+1) - e * sin(E(n+1));
    
    M_error = M - M_check;
    if(abs(M_error) < Threshold)
        break;
    end
    n = n + 1;
end

E_final = E(length(E));