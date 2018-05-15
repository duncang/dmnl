function [phi, theta] = GARD_GravityVector(A_xb,A_yb,A_zb)

phi =  atan2(-A_yb,sqrt(A_xb^2 + A_zb^2));% roll
theta = atan2(A_xb,-A_zb); % pitch


        
        