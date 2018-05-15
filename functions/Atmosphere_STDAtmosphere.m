


function [P T p a] = Atmosphere_STDAtmosphere(alt,R,S,g)
% Set the reference altitudes
refAlt = [0 11000 20000 32000 47000 51000 71000 84900 86000]; Index = 1;
% Check if data is within limits
if (alt<refAlt(1))||(alt>=refAlt(end))
    warning('Altitude is outside standard atmosphere data');
    refSlope = -6.5;    refTemp = 288.15;   refPres = 101325.0;
else
    % Search for the Reference altitude
    while (((refAlt(Index)<= alt)&&(alt<refAlt(Index+1)))~=1)
        Index = Index+(sign(alt-refAlt(Index)));
    end
    % Get reference values from the standard atmosphere data
    switch (Index)
        case 1
            refSlope = -6.5;    refTemp = 288.15;   refPres = 101325.0;
        case 2
            refSlope = 0;       refTemp = 216.65;   refPres = 22632.00;
        case 3
            refSlope = 1;       refTemp = 216.65;   refPres = 5474.900;
        case 4
            refSlope = 2.8;     refTemp = 228.65;   refPres = 868.0187;
        case 5
            refSlope = 0;       refTemp = 270.65;   refPres = 110.9063;
        case 6
            refSlope = -2.8;    refTemp = 270.65;   refPres = 66.93890;
        case 7
            refSlope = -2;      refTemp = 214.65;   refPres = 3.9564;    
        case 8
            refSlope = 0;       refTemp = 186.946;  refPres = 0.3769;
        otherwise
            refSlope = 0;       refTemp = 0;   refPres = 0;
    end
end
% Calculate atmosphere data
if (refSlope == 0)
    T = refTemp;
    P = refPres*exp((-g*(alt-refAlt(Index))/(R*refTemp)));
else
    T = refTemp+(refSlope/1000)*((alt-refAlt(Index)));
    P = refPres*(T/refTemp)^(-g/(R*(refSlope/1000)));
end
p = P/(R*T); a = sqrt(S*R*T);

