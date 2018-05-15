function [quat] = EulerToQuat(euler)
% [quat] = EulerToQuat(euler)
% $Id: EulerToQuat.m 1874 2008-07-15 04:42:16Z n2523710 $
%


phi = euler(1);  %roll
theta = euler(2);  %pitch
psi = euler(3);   %yaw


%note that for this function to give the right +/- signs for the
%quaternions, the yaw angle must be restricted between +/- pi. 

%normalise psi if over 360 deg
if psi > 2*pi    
    
    psi = psi - 2*pi*(floor(psi/(2*pi)));
    
end

%now restrict between +/- pi
if psi >= pi
    
    psi = psi - 2*pi;
    
end



Cphi2 = cos(phi/2);
Ctheta2 = cos(theta/2);
Cpsi2 = cos(psi/2);

Sphi2 = sin(phi/2);
Stheta2 = sin(theta/2);
Spsi2 = sin(psi/2);


q0 = (Cphi2*Ctheta2*Cpsi2 + Sphi2*Stheta2*Spsi2);
q1 = (Sphi2*Ctheta2*Cpsi2 - Cphi2*Stheta2*Spsi2);
q2 = (Cphi2*Stheta2*Cpsi2 + Sphi2*Ctheta2*Spsi2);
q3 = (Cphi2*Ctheta2*Spsi2 - Sphi2*Stheta2*Cpsi2);



quat(1) = q0;
quat(2) = q1;
quat(3) = q2;
quat(4) = q3;


% %test code for the restriction algorithm
% psiin(i) = 0;
% for i = 1:4000
%    
%     
%     psiinstore(i) = psiin(i);
%      
%     if psiin(i) > 2*pi
% 
%         psiin(i) = psiin(i) - 2*pi*(floor(psiin(i) /(2*pi)));
% 
%     end
% 
%     %now restrict between +/- pi
%     if  psiin(i) >= pi
%          psiin(i) =  psiin(i) - 2*pi;
%     end
%     
%   psiout(i) = psiin(i) ; 
%   
%    psiin(i+1) = psiinstore(i) + 1*pi/180;
% 
% end
% 
% plot(psiinstore*180/pi)
% hold;
% plot(psiout*180/pi)
%     




