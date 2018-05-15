% pbias calculation for the 5-in-view case


clear

Pmd = 0.001;
Td = 4.8916; % calculated using Threshold_Calc.m

iters = 0;
error = 1;
step = 1;

ll_try = 0;

while (abs(error) > 1e-6)
    iters = iters + 1;
    if (error > 0)
        ll_try = ll_try + step/iters ;
    else
        ll_try = ll_try - step/iters ;
    end

    pbias = ll_try + Td;
        
    Pcalc = quad(@normpdf,pbias - Td,pbias + Td,1e-15);
    
    ll_save(iters) = ll_try;
    P_save(iters) = Pcalc;
    
    error = Pcalc - Pmd;
    
    e_save(iters) = error;
    

    

    if iters > 20
        disp('Iterations exceeded 20 - stopping!');
        break;
    end
end


