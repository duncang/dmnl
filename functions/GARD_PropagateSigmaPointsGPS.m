function [xs_0_k,xs_i_k] = GARD_PropagateSigmaPointsGPS(Na,xs_0,xs_i,gps_dt,NumberStates,ProcessNoiseStates)




    % propagate system dynamics
    for i=0:2*Na
        if i==0
            %% zeroth point
            
           
            
           
            
            
            
            % position
            xs_0_k(1) = xs_0(1) + xs_0(4)*gps_dt  + (xs_0(9)*gps_dt^2)/2 + (xs_0(12)*gps_dt^3)/3;
            xs_0_k(2) = xs_0(2) + xs_0(5)*gps_dt  + (xs_0(10)*gps_dt^2)/2 + (xs_0(13)*gps_dt^3)/3;
            xs_0_k(3) = xs_0(3) + xs_0(6)*gps_dt  + (xs_0(11)*gps_dt^2)/2 + (xs_0(14)*gps_dt^3)/3;
            xs_0_k(7) = xs_0(7) + xs_0(8)*gps_dt;
            
            % velocity
            xs_0_k(4) = xs_0(4)  + xs_0(9)*gps_dt + (xs_0(12)*gps_dt^2)/2;
            xs_0_k(5) = xs_0(5)  + xs_0(10)*gps_dt + (xs_0(13)*gps_dt^2)/2;
            xs_0_k(6) = xs_0(6)  + xs_0(11)*gps_dt + (xs_0(14)*gps_dt^2)/2;
            xs_0_k(8) = xs_0(8);
            
            % acceleration
            xs_0_k(9) = xs_0(9) + xs_0(12)*gps_dt;
            xs_0_k(10) = xs_0(10) + xs_0(13)*gps_dt;
            xs_0_k(11) = xs_0(11) + xs_0(14)*gps_dt;
            
             % jerk
            xs_0_k(12) = xs_0(12);
            xs_0_k(13) = xs_0(13);
            xs_0_k(14) = xs_0(14);
            
            %% add process noise
            V = xs_0(NumberStates+1:NumberStates+ProcessNoiseStates)';
            xs_0_k(1:NumberStates) = xs_0_k(1:NumberStates) + V;
            
            
        else
            
            % i-points
           
            
           
            
            % position
            xs_i_k(1,i) = xs_i(1,i) + xs_i(4,i)*gps_dt + (xs_i(9,i)*gps_dt^2)/2 + (xs_i(12,i)*gps_dt^3)/3;
            xs_i_k(2,i) = xs_i(2,i) + xs_i(5,i)*gps_dt + (xs_i(10,i)*gps_dt^2)/2 + (xs_i(13,i)*gps_dt^3)/3;
            xs_i_k(3,i) = xs_i(3,i) + xs_i(6,i)*gps_dt + (xs_i(11,i)*gps_dt^2)/2 + (xs_i(14,i)*gps_dt^3)/3;
            xs_i_k(7,i) = xs_i(7,i) + xs_i(8,i)*gps_dt;

            % velocity
            xs_i_k(4,i) = xs_i(4,i) + xs_i(9,i)*gps_dt + (xs_i(12,i)*gps_dt^2)/2;
            xs_i_k(5,i) = xs_i(5,i) + xs_i(10,i)*gps_dt + (xs_i(13,i)*gps_dt^2)/2;
            xs_i_k(6,i) = xs_i(6,i) + xs_i(11,i)*gps_dt + (xs_i(14,i)*gps_dt^2)/2;
            xs_i_k(8,i) = xs_i(8,i);
            
             % acceleration
            xs_i_k(9,i) = xs_i(9,i) + xs_i(12,i)*gps_dt;
            xs_i_k(10,i) = xs_i(10,i) + xs_i(13,i)*gps_dt;
            xs_i_k(11,i) = xs_i(11,i) + xs_i(14,i)*gps_dt;
            
             % jerk
            xs_i_k(12,i) = xs_i(12,i);
            xs_i_k(13,i) = xs_i(13,i);
            xs_i_k(14,i) = xs_i(14,i);
             
            
            %% add process noise
            V = xs_i(NumberStates+1:NumberStates+ProcessNoiseStates,i);
            xs_i_k(1:NumberStates,i) = xs_i_k(1:NumberStates,i) + V;
            
        end
    end
    
    
    
    