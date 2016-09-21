% This functions crops a 3d volume using the boundaries provided by
% CropLabelVol
% Icropped=applyCropping(I,cropping)
function Icropped=applyCropping(I,cropping)


i1=cropping(1);
j1=cropping(2);
k1=cropping(3);

i2=cropping(4);
j2=cropping(5);
k2=cropping(6);

Icropped=I(i1:i2,j1:j2,k1:k2);