function range = pseudorange(SVpos,U)

% calculate the pseudorange range between the satellite position SV and
% estimated user position U

c = 2.99792458e8;  % speed of light m/s

pointingvector = SVpos - U(1:3);


%range = magnitude(pointingvector) + c * U(4);
range = magnitude(pointingvector) + U(4);