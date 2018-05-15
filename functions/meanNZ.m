function [MEAN] = meanNZ(signal);

%By Troy Bruggemann (c) 2005


%This function calculates the mean of a set of data excluding any zero elements.
%This is useful for GARDSIm when a set of position data has 0's in it due to lack of a solution.
%Signal - the data, must be a row vector 
%MEAN - the mean of the data

k = 1;
h = size(signal);
for i = 1:(h(2)-1)
         
    if signal(i) ~= 0           
          newsignal(k) = signal(i);
       k = k+1;
    end
    
end

hh = size(newsignal);
MEAN = sum(newsignal)/hh(2);
