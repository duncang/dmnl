function [Q] = skew_symmetric(omega1,omega2,omega3)
% function [Q] = skewbody(q)
% Returns the skew symmetric form of vector
%
% Written by Troy B. 5.2.09
%




 Q = [  0, -omega3, omega2 ; 
     omega3, 0, -omega1;
      -omega2, omega1, 0;];
            