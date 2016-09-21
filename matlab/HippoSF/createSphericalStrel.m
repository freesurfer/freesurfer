%
% strel=createSphericalStrel(r,pixsize)
%
% Creates a spherical structuring element or radius r
% If pixsize is provided, the voxels dimensions are taken into
% consideration
%
function strel=createSphericalStrel(r,pixsize)

if exist('pixsize','var')==0
    pixsize=[1 1 1];
end
if numel(pixsize)==1
    pixsize=[pixsize pixsize pixsize];
end

d=ceil(2*r./pixsize+1);
d=d+mod(d+1,2);
strel=zeros(d);
for ii=1:d(1)
    for jj=1:d(2)
        for kk=1:d(3)
            if (((ii-(d(1)+1)/2)*pixsize(1))^2+((jj-(d(2)+1)/2)*pixsize(2))^2+((kk-(d(3)+1)/2)*pixsize(3))^2)<=r^2
                strel(ii,jj,kk)=1;
            end
        end
    end
end
strel=logical(strel);

