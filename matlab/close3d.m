function vout =  close3d(vin)
% vout = function close3d(vin)

[width, height, depth] = size(vin);
vout = dilate3d(vin) ;
vout = erode3d(vout) ;
