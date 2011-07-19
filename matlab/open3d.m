function vout = open3d(vin)
% vout = function open3d(vin)

[width, height, depth] = size(vin);
vout = erode3d(vin) ;
vout = dilate3d(vout) ;
