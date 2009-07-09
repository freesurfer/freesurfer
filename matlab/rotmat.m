function M = rotmat(deg)
% M = rotmat(deg)

rads = radians(deg) ;

M = eye(3) ;
M(1,1) = cos(rads) ;
M(1,2) = sin(rads) ;
M(2,1) = -sin(rads) ;
M(2,2) = cos(rads) ;

