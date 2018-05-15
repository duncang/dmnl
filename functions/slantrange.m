function range = slantrange(vector1,vector2)

% calculate the slant range between 2 position vectors

pointingvector = vector2 - vector1;

range = sqrt(pointingvector(1)^2 + pointingvector(2)^2 + pointingvector(3)^2);