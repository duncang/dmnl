function psi_out = pibound(psi)

if(psi > pi)
    psi_out = psi - 2*pi;
elseif(psi < -pi)
    psi_out = psi + 2*pi;
else
    psi_out = psi;
end

