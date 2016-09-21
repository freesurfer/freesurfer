% This functions crops a volume around the voxels greater than 'threshold' 
% It allows a margin given by 'margin'.
% Defaults are margin=10, threshold=0
% [cropped cropping]=cropLabelVol(M,margin,threshold)
function [cropped,cropping]=cropLabelVol(V,margin,threshold)

if exist('margin','var')==0 % default to 10 voxels
    margin=10;
end
if exist('threshold','var')==0 
    threshold=0;
end
if numel(margin)==1
    margin=repmat(margin,[1 3]);
end


f=find(V>threshold);
[i,j,k]=ind2sub(size(V),f);

i1=max(1,min(i)-margin(1));
j1=max(1,min(j)-margin(2));
k1=max(1,min(k)-margin(3));

i2=min(size(V,1),max(i)+margin(1));
j2=min(size(V,2),max(j)+margin(2));
k2=min(size(V,3),max(k)+margin(3));

cropping=[i1 j1 k1 i2 j2 k2];
cropped=V(i1:i2,j1:j2,k1:k2);

