
clear

% calculates detection threshold Td for the 5-in-view case

Pfa = 1.0e-6/2;
T_try = 1;
iters = 0;
error = 1;
step = 2;
overshoot = 0;

while (abs(error) > 1e-13)
    iters = iters + 1;
    if (error > 0)
        T_try = T_try + step/iters ;
    else
        T_try = T_try - step/iters ;
    end

    
        
    Pcalc = quad(@normpdf,T_try,100,1e-15);
    T_save(iters) = T_try;
    P_save(iters) = Pcalc;
    
    error = Pcalc - Pfa;
    
    e_save(iters) = error;
    

    

    if iters > 600
        disp('Iterations exceeded 20 - stopping!');
        break;
    end
end


