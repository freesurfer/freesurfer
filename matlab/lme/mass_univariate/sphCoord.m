function [phi,theta] = sphCoord(Coord3D)
% [phi,theta] = sphCoord(Coord3D)
% 
% This function transforms spherical 3D XYZ spherical coordinates to 2D 
% (theta,phi) spherical coordinates.
%
% Input
% Coord3D: Data Matrix with the spherical XYZ coordinates at the first
% three columns respectively.
%
% Output
% theta: Latitude  
% phi: Longitud.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:11 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:11 $
%    $Revision: 1.1.2.2 $

% [theta,phi,r] = cart2sph(overlayData(:,1),overlayData(:,2),overlayData(:,3));
% //It was confusion

X=Coord3D(:,1); Y=Coord3D(:,2); Z=Coord3D(:,3);
phi=atan2(Y,X);
theta=atan2(sqrt(X.^2+Y.^2),Z);

