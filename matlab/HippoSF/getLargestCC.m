% extracts the biggest connected component from a binary image (2D or 3D)
% Y=getLargestCC(X)
%
function Y=getLargestCC(X)

X=X>0;  % in case...

if length(size(X))==2
    [LAB n_lab]=bwlabel(X,4);
else
    [LAB n_lab]=bwlabeln(X,6);
end

volumes=zeros(1,n_lab);
for v=1:n_lab
   volumes(v)=sum(LAB(:)==v); 
end
pos_max=find(volumes==max(volumes));
if ~isempty(pos_max)
    Y=LAB==pos_max(1);
else
    Y=logical(zeros(size(X)));
end
