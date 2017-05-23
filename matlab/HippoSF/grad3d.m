% [Gmodule GX GY GZ]=grad3d(X)
% Gradient for 3D volumes
function [Gmodule GX GY GZ]=grad3d(X)

v=[-1 0 1];

GX=imfilter(X,reshape(v,[3 1 1]));
GY=imfilter(X,reshape(v,[1 3 1]));
GZ=imfilter(X,reshape(v,[1 1 3]));

GX(1,:,:)=X(2,:,:)-X(1,:,:); GX(end,:,:)=X(end,:,:)-X(end-1,:,:);
GY(:,1,:)=X(:,2,:)-X(:,1,:); GY(:,end,:)=X(:,end,:)-X(:,end-1,:);
GZ(:,:,1)=X(:,:,2)-X(:,:,1); GZ(:,:,end)=X(:,:,end)-X(:,:,end-1);

Gmodule=sqrt(GX.*GX+GY.*GY+GZ.*GZ);
