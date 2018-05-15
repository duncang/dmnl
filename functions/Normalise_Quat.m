function [q0norm,q1norm,q2norm,q3norm] = Normalise_Quat(q0,q1,q2,q3);
%this function normalises the quaternions to a unit quaternion


normquat = norm([q0,q1,q2,q3]);


q0norm = q0/normquat;
q1norm = q1/normquat;
q2norm = q2/normquat;
q3norm = q3/normquat;




