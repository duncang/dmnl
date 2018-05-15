function C_BN_new = GARD_DCMUpdate(C_BN,del_phi,del_theta,del_psi,dt)

OMEGAedot = 7.2921151467e-005;

OMEGA_b = [ 0         -del_psi     del_theta; ...
        del_psi    0          -del_phi; ...
       -del_theta    del_phi     0        ];

% find the skew symmetric form of the navigation frame inertial
% rotation rate
OMEGA_in = [0, -OMEGAedot * sin(-27*pi/180), 0; ...
            OMEGAedot * sin(-27*pi/180), 0, -OMEGAedot * cos(-27*pi/180);
            0, OMEGAedot * cos(-27*pi/180), 0];


% Direction Cosine matrix attitude update
C_BN_new = C_BN + (C_BN * OMEGA_b - OMEGA_in * C_BN * dt);

C_BN_new = GARD_OrthogonaliseDCM(C_BN_new);