% This function requires the Google Earth Toolbox from
% http://www.mathworks.com/matlabcentral/fileexchange/loadFile.do?objectId=12954
%
% to plot a 3d path in google earth from flight data - NOTE order is
% LONGITUDE, LATITUDE, HEIGHT in decimal degrees and metres

kmlstr3d = ge_plot3(pos_truth_llh(3,:)*180/pi,pos_truth_llh(2,:)*180/pi,pos_truth_llh(4,:),'altitudeMode','absolute');
ge_output('demo3d.kml',kmlstr3d);
