

%clear;

Pmd = 0.001;
%Pfa = 1e-6;


for N = 6:14;
    DOF = N-4;

    %Td = 5.2565; % note - this is dependant on the DOF.

    a = Td(N)^2;

    lambda_try = 100;
    iters = 0;
    error = 1;
    step = 10;
    while (abs(error) > 1e-9)
        iters = iters+1;

        if error > 0
            lambda_try = lambda_try - step/iters;
        else
            lambda_try = lambda_try + step/iters;
        end

        Pcalc = ncx2cdf(a,DOF,lambda_try);

        error = Pmd - Pcalc;

        if iters > 2000
            disp('Iterations exceeded 20 - stopping!');
            break;
        end

    end

    iters;
    lambda_try;
    pbias(N) = sqrt(lambda_try)
end

