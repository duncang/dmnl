function y = normpdf2(mu_x,sigma_x,mu_y,sigma_y,rho,X,Y)
% function y = normpdf2(mu_x,sigma_x,mu_y,sigma_y,rho,X,Y)
% Implements 2-dimentional (bivariate) normal distribution
% http://en.wikipedia.org/wiki/Multivariate_normal_distribution


a = 1/(2*pi*sigma_x*sigma_y*sqrt(1-rho^2));

b = (-1/(2 * 1-rho^2)) * ( ((X - mu_x)^2)/sigma_x^2 + (((Y - mu_y)^2)/sigma_y^2) - (2 * rho * (X - mu_x) * (Y - mu_y) / (sigma_x*sigma_y))  );

y = a * exp(b);

