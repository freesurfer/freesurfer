function dR = dMRIangles2RotMat(angles,k)
% dR = dMRIangles2RotMat(angles)
% Derivative of rotation matrix wrt the kth angle
% Convert 3 euler angles into a 4x4 rotation matrix
% angles is a 3x1 vector in radians. 
% angles(1) - pitch - rotation about x or LR axis (gamma)
% angles(2) - yaw   - rotation about y or AP axis (beta)
% angles(3) - roll  - rotation about z or SI axis (alpha)
%  Ref: Craig, Intro to Robotics
% This function matches MRIangles2RotMat() in transform.c
% R = Rz*Ry*Rx;

gamma = angles(1);
beta  = angles(2);
alpha = angles(3);

Rx = zeros(3,3);
Rx(1,1) = +1;
Rx(2,2) = +cos(gamma);
Rx(2,3) = -sin(gamma);
Rx(3,2) = +sin(gamma);
Rx(3,3) = +cos(gamma);

Ry = zeros(3,3);
Ry(1,1) = +cos(beta);
Ry(1,3) = +sin(beta);
Ry(2,2) = 1;
Ry(3,1) = -sin(beta);
Ry(3,3) = +cos(beta);

Rz = zeros(3,3);
Rz(1,1) = +cos(alpha);
Rz(1,2) = -sin(alpha);
Rz(2,1) = +sin(alpha);
Rz(2,2) = +cos(alpha);
Rz(3,3) = +1;

dR = zeros(3,3);
if(k==1)
  dR(1,1) = +0;
  dR(2,2) = -sin(gamma);
  dR(2,3) = -cos(gamma);
  dR(3,2) = +cos(gamma);
  dR(3,3) = -sin(gamma);
  dR = Rz*Ry*dR;
end
if(k==2)
  dR(1,1) = -sin(beta);
  dR(1,3) = +cos(beta);
  dR(2,2) = 0;
  dR(3,1) = -cos(beta);
  dR(3,3) = -sin(beta);
  dR = Rz*dR*Rx;
end
if(k==3)
  dR(1,1) = -sin(alpha);
  dR(1,2) = -cos(alpha);
  dR(2,1) = +cos(alpha);
  dR(2,2) = -sin(alpha);
  dR(3,3) = 0;
  dR = dR*Ry*Rx;  
end

M = zeros(4,4);
M(1:3,1:3) = dR;
dR = M;

return;

