function [ListVerticesInPialROI,ListFacesInPialROI] = PropagateGeodesic(mesh_total,m,iV,wholepath)



[centerIndexMT,centerValuesMT]=mesh_vertex_nearest(mesh_total.vertices,m.vertices(iV,:));
index=mesh_total.facesOfVertex(centerIndexMT).faceList;
clear centerIndexMT;
clear centerValuesMT;

t=mesh_total.faces(index,:); 

% pour trouver la liste des faces qui contiennent le vertex le plus proche,
% ensuite il faut enlever successivement toutes les lignes qui contiennent
% au moins un des vertex

more off
l=unique(t(:));
l2=l;
oldsize=1;
while size(l,1)>oldsize
    fprintf('.')
    index2=[];
    for r = l2' 
        index2 = [index2, mesh_total.facesOfVertex(r).faceList];
    end
    index2 = unique (index2);
    for u = 1: size(index2,2)
        if isInGeodesicROI(mesh_total,wholepath,index2(u))
            t2=mesh_total.faces(index2(u),:);
            t=[t; t2];
            t=unique(t,'rows');
        %else
        %index2(u) 
        %fprintf('one vertex not added .. ')
         end               
    end
    oldsize=size(l);
    oldl=l;
    l=unique(t(:));
    l2=setdiff(l,oldl);
    %size(l2)
    l=[oldl ; l2];
end

disp('region delimited!')

ListVerticesInPialROI=unique(l);

ListFacesInPialROI = t;
