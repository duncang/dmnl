function q_out = quatupdate(q_in,del_phi,del_theta,del_psi)

OMEGA = skewbody(del_phi,del_theta,del_psi);

%expm(-0.5*OMEGA)

s = 0.5*sqrt(del_phi^2+del_theta^2+del_psi^2);
del_q = eye(4,4)*cos(s) - 0.5 * OMEGA * sin(s)/s;

q_out = del_q * q_in;

%q_out = q_out/norm(q_out);

