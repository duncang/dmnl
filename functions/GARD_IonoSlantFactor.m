function slant_factor = GARD_IonoSlantFactor(Elevation)

% Convert Elevation in radians to semi-circles
E = Elevation / pi;

% calculate the slant factor, F
slant_factor = 1 + 16 * (0.53 - E) ^3;


