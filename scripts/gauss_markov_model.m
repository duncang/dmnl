

gm_var = 0.05;
tau = 300;
sequence_length = 10000;

noise = randn(1,sequence_length)*sqrt(gm_var);
x(1) = noise(1);
y(1) = 1;
for t=2:sequence_length
    y(t) = t;
%    x(t) = x(t-1) + noise(t) - 1/tau * x(t-1);
%    variance(t) = sqrt(t) * sqrt(gm_var);
    [x(t),variance(t)] = GaussMarkov_Process(x(t-1), 1/tau,sqrt(gm_var),1);
end

figure();hold on;
plot(noise,'r');
plot(x);
plot(sqrt(variance),'g');
plot(-sqrt(variance),'g');

%plot(filter(ones(1,500)/500,1,x),'g');
hold off;


figure;
plot(abs(xcorr(x,x))/sequence_length);