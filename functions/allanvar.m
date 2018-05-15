function avar = allanvar(data,dt,tau);
% function avar = allanvar(data,tau);
% Calculates the Allan Variance of data over integration time tau
% written by Duncan Greer Febuary 2006
%
%
% data - a column vector containing the data set
% dt - the time step of the data
% tau - the averaging time
%


% allanvar(dTpos(1,:),0.02,0.1);
% tau = 1;
% dt = 0.02;
% data = dTpos(1,:);
% get the number of available bins
binlength = floor(tau/dt);
n = floor(length(data) / binlength);
if(n < 9)
    disp('Warning - Data set is too short.  At least 9 bins should be available');
end



for index = 1:n
    % get the bin data set
    binset = data((index-1)*binlength+1:index*binlength);
    binmean(index) = mean(binset);
end

sum = 0;

for index = 1:n-1
    sum = sum + (binmean(index+1) - binmean(index))^2;
end

avar_2 = sum / (2 * (n-1));

avar = sqrt(avar_2);