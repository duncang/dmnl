function [K_out] = CEP_TableLookup_PmdH0_H(value);

%for rare normal case, (fault free) 
% Pmd_H0_H = 1e-10;
%gets the relevant K value for c value (ratio) ie sigma_x/sigma_y between 0.0 and 1.0 where
%sigma_y is the largest
%only relevant for Pmd value of 1- Pmd = 1-  1e-10 = 0.9999999999


%note that value is chosen conservatively, meaning values might be 1 or 2
%metres higher in the HPL (HPL = K*sigma) than they could be otherwise

if value >= 0 & value < 0.1

K_out = 5.5;

elseif value >=0.1 & value < 0.2
K_out = 5.5;

elseif value >= 0.2 & value < 0.3
K_out = 5.5;

elseif value >= 0.3 & value < 0.4
K_out = 5.5;

elseif value >= 0.4 & value < 0.5
K_out = 5.5;

elseif value >= 0.5 & value < 0.6
K_out = 5.5;

elseif value >= 0.6 & value < 0.7
K_out = 5.6;

elseif value >= 0.7 & value < 0.8
K_out = 5.6;

elseif value >= 0.8 & value < 0.9
K_out = 5.7;

elseif value >= 0.9 & value < 1.0
K_out = 5.7;
    
    
elseif value == 1.0
K_out = 5.8;    
    

end


