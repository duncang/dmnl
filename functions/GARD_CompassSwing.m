function [Xsf,Ysf,Xoff,Yoff] = GARD_CompassSwing(swing)
% function [Xsf,Ysf,Xoff,Yoff] = GARD_CompassSwing(swing)
% swing = n x 3 vector of compass measurements swung by 360 deg
%
%

Xmin = min(swing(:,1));
Xmax = max(swing(:,1));

Ymin = min(swing(:,2));
Ymax = max(swing(:,2));

Xsf = max(1.0, (Ymax - Ymin) / (Xmax - Xmin) );
Ysf = max(1.0, (Xmax - Xmin) / (Ymax - Ymin) );

Xoff = ((Xmax - Xmin)/2 - Xmax) * Xsf; 
Yoff = ((Ymax - Ymin)/2 - Ymax) * Ysf;

disp('Compass Calibration Results:');
disp(sprintf('   Xsf = %f',Xsf));
disp(sprintf('   Ysf = %f',Ysf));
disp(sprintf('   Xoff = %f',Xoff));
disp(sprintf('   Yoff = %f',Yoff));

