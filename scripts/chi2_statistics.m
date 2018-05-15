

clear;

x = [0:0.1:50];

% central chi-squared distribution
for v=1:6
    y(:,v) = chi2pdf(x,v); 
end


% non-central chi-squared distribution
lambda = 10*ones(1,length(x));
for v=1:6
    yn(:,v) = ncx2pdf(x,v,lambda); 
end


