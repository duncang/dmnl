function [K_out] = CEP_TableLookup_Pmd(value);

%gets the relevant K value for c value (ratio) ie sigma_x/sigma_y between 0.0 and 1.0 where
%sigma_y is the largest
%only relevant for Pmd value of 1- Pmd = 1- 0.001 = 0.99900
%note that value is chosen conservatively, meaning values might be 1 or 2
%metres higher in the HPL (HPL = K*sigma) than they could be

if value >= 0 & value < 0.1

K_out = 3.3;

elseif value >=0.1 & value < 0.2
K_out = 3.3;

elseif value >= 0.2 & value < 0.3
K_out = 3.4;

elseif value >= 0.3 & value < 0.4
K_out = 3.4;

elseif value >= 0.4 & value < 0.5
K_out = 3.4;

elseif value >= 0.5 & value < 0.6
K_out = 3.4;

elseif value >= 0.6 & value < 0.7
K_out = 3.4;

elseif value >= 0.7 & value < 0.8
K_out = 3.5;

elseif value >= 0.8 & value < 0.9
K_out = 3.6;

elseif value >= 0.9 & value < 1.0
K_out = 3.8;
    
    
elseif value == 1.0
K_out = 3.8;    
    

end


