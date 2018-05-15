function P_dual_out = GARD_DualCovariancePropagator(P_dual_in,PHI_k,Q_k,K0,Kn,H0,R,i)

s0 = size(K0,1);
sn = size(Kn,1);

if s0 ~= sn
    %Kn not padded with zeros
    Kn = [Kn(1:17,1:i-1), zeros(17,1),Kn(1:17,i:7),Kn(1:17,8:8+i-2),zeros(17,1),Kn(1:17,8+i-1:14),Kn(1:17,15)];

end

PHI_kk = [PHI_k,        zeros(s0,s0);
          zeros(s0,s0), PHI_k       ];
Q_kk = [Q_k,          zeros(s0,s0);
        zeros(s0,s0), zeros(s0,s0)];

P_dual_up = PHI_kk * P_dual_in * PHI_kk' + Q_kk;


PHI = [eye(s0) - K0*H0', zeros(s0,s0);
       (Kn-K0)*H0'     , eye(s0) - Kn*H0'];

GAMMA = [K0;
         K0 - Kn];
         
P_dual_out = PHI * P_dual_up * PHI' + GAMMA * R * GAMMA';
         