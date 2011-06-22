function R = MRIangles2RotMat(angles)
% R = MRIangles2RotMat(angles)
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

% This will be a 3x3 matrix
%R3 = MatrixMultiply(Rz,Ry,NULL);
%R3 = MatrixMultiply(R3,Rx,R3);
R3 = Rz*Ry*Rx;

R = zeros(4,4);
R(1:3,1:3) = R3;
R(4,4) = 1;

return;

% This is just a check
R2 = zeros(4,4);
R2(4,4) = 1;

R2(1,1) = cos(alpha)*cos(beta);
R2(1,2) = -sin(alpha)*cos(gamma) + cos(alpha)*sin(beta)*sin(gamma);
R2(1,3) = sin(alpha)*sin(gamma) + cos(alpha)*sin(beta)*cos(gamma);

R2(2,1) = sin(alpha)*cos(beta);
R2(2,2) = cos(alpha)*cos(gamma) + sin(alpha)*sin(beta)*sin(gamma);
R2(2,3) = -cos(alpha)*sin(gamma) + sin(alpha)*sin(beta)*cos(gamma);

R2(3,1) = -sin(beta);
R2(3,2) = cos(beta)*sin(gamma);
R2(3,3) = cos(beta)*cos(gamma);




