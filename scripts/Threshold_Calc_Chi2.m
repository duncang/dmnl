

% calculate Td for 2 or more DOF using chi-square statistics
clear;

Pfa = 1e-6;
    
for N = 6:14
    DOF = N-4;



    a_try = 1;
    error = 1;
    iters = 0;
    step = 10;
    while (abs(error) > 1e-10)

        iters = iters + 1;
        if (error > 0)
            a_try = a_try - step/iters ;
        else
            a_try = a_try + step/iters ;
        end

        Pcalc = 1-chi2cdf(a_try,DOF);
        error = Pfa - Pcalc;

        if iters > 2000
            disp('Iterations exceeded 20 - stopping!');
            break;
        end

    end

    iters;

    a_try;

    Td(N) = sqrt(a_try)
    
end

