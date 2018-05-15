function [euler] = QuatToEuler(quat)
% function [euler] = QuatToEuler(quat)
%
% $Id: QuatToEuler.m 1883 2008-07-15 05:53:55Z n2523710 $
%


   q0 = quat(1);
   q1 = quat(2);
   q2 = quat(3);
   q3 = quat(4);
   
   
   
   phi = atan2( 2*(q0*q1 + q2*q3), q0^2 + q3^2 - q1^2 -q2^2);
   theta = asin(2*(q0*q2 - q1*q3));
   psi = atan2( 2*(q0*q3 + q1*q2), q0^2 + q1^2 - q2^2 - q3^2);
   
   
   euler(1) = phi;
   euler(2) = theta;
   euler(3) = psi;