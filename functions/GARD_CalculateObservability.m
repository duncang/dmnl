function [ObsRank,ObsMatrix] = GARD_CalculateObservability(F,H);


HSize = size(H,1);

clear ObsMatrix;
for i=1:NumberStates
   ObsMatrix((i-1)*HSize+1:i*HSize,1:NumberStates) = H *  F^(i-1);
end

ObsRank = rank(ObsMatrix);   