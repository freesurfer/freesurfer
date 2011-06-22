function angles = MRIrotMat2angles(R)
% angles = MRIrotMat2angles(R)
% Convert 4x4 rotation matrix to 3 euler angles
%
% angles is a 3x1 vector in radians. 
% angles(1) - pitch - rotation about x or LR axis (gamma)
% angles(2) - yaw   - rotation about y or AP axis (beta)
% angles(3) - roll  - rotation about z or SI axis (alpha)
%
% Ref: Craig, Intro to Robotics, first Fixed Angle Set:
%   Rxyz(g,b,a) = Rz(a)*Ry(b)*Rx(g);
% Only consistent withing +/- pi/2

beta = asin(-R(3,1));
alpha = atan2(R(2,1)/cos(beta),R(1,1)/cos(beta));
gamma = atan2(R(3,2)/cos(beta),R(3,3)/cos(beta));

angles = [gamma beta alpha]';

%c2beta = R(1,1).^2 + R(2,1).^2;
%beta2 = acos(sqrt(c2beta))
%keyboard

return;

