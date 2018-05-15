

N = 5;
sigma = 33;
pbias = sigma * 7.9818;
td = sigma * 4.8916;

G =  [-0.8607445743 -0.3446039300 0.3746557209 1;
       0.2109370676  0.3502943374 0.9125784518 1;
      -0.0619331310 -0.4967359072 0.8656891623 1;
      -0.7248969588  0.4759681238 0.4979746422 1;
      -0.4009266538  0.1274180997 0.9072058455 1];

A = inv(G'*G)*G';
S = eye(N,N) - G*A;

for i=1:N
    slope(i) = sqrt(A(1,i)^2 + A(2,i)^2) / sqrt(S(i,i));
end

[slope_MAX,slope_i] = max(slope);

[Q,R] = qr(G);
P = Q(:,5:N)'; % lower n-4 rows of Q';

rbias = pbias / norm(P(:,slope_i));

e = zeros(N,1);
e(slope_i) = rbias;

E = A * e;

HPL = norm(E(1:2));

HPL_2 = pbias * slope_MAX;

   
 