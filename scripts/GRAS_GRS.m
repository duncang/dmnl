% GRAS reference station (GRS) locations
% Written by Dundan Greer 24 January 2006
%
% Source: Proposed GRAS Architecture - AsA website
%
% $Id: GRAS_GRS.m 1884 2008-07-15 05:54:33Z n2523710 $
%
f2m = 1/3.2;
d2r = pi/180;

Weipa = [degmin2deg(-12, 40.7), degmin2deg(141,55.5),63*f2m];
Mackay = [degmin2deg(-21, 10.3),degmin2deg(149,10.8),19*f2m];
Brisbane = [degmin2deg(-27,23.0),degmin2deg(153,07.1),13*f2m];
Sydney = [degmin2deg(-33, 56.8),degmin2deg(151, 10.6),21*f2m];
Melbourne = [degmin2deg(-37, 40.4), degmin2deg(144, 50.6),434*f2m];
Launceston = [degmin2deg(-41, 32.7), degmin2deg(147, 12.9),562*f2m];
Adelaide = [degmin2deg(-34, 56.7), degmin2deg(138, 31.8),20*f2m];
AliceSprings = [degmin2deg(-23, 48.4), degmin2deg(133, 54.1),1789*f2m];
Darwin = [degmin2deg(-12, 24.9), degmin2deg(130, 52.6),103*f2m];
Perth = [degmin2deg(-31, 56.4), degmin2deg(115, 58.8),67*f2m];
Carnarvon = [degmin2deg(-24, 52.8), degmin2deg(113, 40.3),13*f2m];
Derby = [degmin2deg(-17, 22.2), degmin2deg(123,39.6),24*f2m];

GRS_LLH=[Weipa;Mackay;Brisbane;Sydney;Melbourne;Launceston; ...
    Adelaide;AliceSprings;Darwin;Perth;Carnarvon;Derby];

for Index=1:length(GRS_LLH)
    GRS_ECEF(Index,:)=LLH2ECEF(GRS_LLH(Index,1)*d2r, ...
                               GRS_LLH(Index,2)*d2r, ...
                               GRS_LLH(Index,3));
end

