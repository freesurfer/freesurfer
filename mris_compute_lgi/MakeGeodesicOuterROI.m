function ListFacesInHullROI = MakeGeodesicOuterROI (mesh_outer, iV, facesInROI);

FacesByLines=[];
FacesByLines=mesh_outer.faces(facesInROI,:);

[index,value]=find(FacesByLines==iV);
index=index';
t=FacesByLines(index,:);

more off
l=unique(t(:));
l2=l;
oldsize=1;
while size(l,1)>oldsize
        for r=1:size(l2,1)
        [index2,value2]=find(FacesByLines==l2(r));
        index2=index2';
        index=[index index2];
        t2=FacesByLines(index2,:);
        t=[t; t2];
        end 
     t=unique(t,'rows');
     oldsize=size(l);
     oldl=l;
     l=unique(t(:));
     l2=setdiff(l,oldl);
     indexListHull=unique(index);
end

ListFacesInHullROI=t;

if size(facesInROI,2) ~= size(indexListHull,2)
  disp([' non geodesic regions removed from the patch '])
end
