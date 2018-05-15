

clear;
interval = 0.01;
sigma_x = 1.0;
sigma_y = 1.0;
mu_x = 0.0;
mu_y = 0.0;
rho = 0.0;



x = [-10:interval:10];
y = [-10:interval:10];

MU = [mu_x mu_y];
SIGMA = [sigma_x^2 rho*sigma_x*sigma_y;rho*sigma_x*sigma_y sigma_y^2];

for i=1:length(x)
    for j = 1:length(y)
        %z(i,j) = normpdf2(mu_x,sigma_x,mu_y,sigma_y,rho,x(i),y(j));
        z(i,j) = mvnpdf([x(i) y(j)],MU,SIGMA);
    end
end

%surf(x,y,z);




for i=1:length(x)
    for j = 1:length(y)
        cdf_z(i,j) = mvncdf([x(i) y(j)],MU,SIGMA);
    end
end

%surf(x,y,cdf_z);


[V,D] = eigs(SIGMA);

figure(); grid on; hold on;


% now, to find the values of x and y which correspond to a probability of
% Pfa - for uncorrelated x and y, the magnitude (sqrt (x^2 + y^2) should be
% constant
Pfa = 0.001;

P_tail = 1-Pfa/2;

[I,J] = find(cdf_z > P_tail);

% print results
P_tail
x(I)
y(J)
sqrt(x(I)^2 + y(J)^2)



norminv(1-Pfa/2,MU,SIGMA)

% integrate along dK then revolve through 2pi

K = 3.3;
dK = interval;
sum = 0;
for k=K:dK:10
    sum = sum + mvnpdf([k 0],MU,SIGMA);
end


angular_interval = 0.1*pi/180;
intervals = 2*pi / angular_interval;
dphi = K * sin(angular_interval);
   
sum2 = sum * dphi * intervals ;

1-sum2


