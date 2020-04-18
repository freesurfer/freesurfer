function [AdjM,cn] = AdjMtx(Surf,maskvtx)
% AdjM = AdjMtx(Surf,maskvtx)
% 
% This function finds the adjacent vertices for all vertices along the 
% surface.
%
% Input
% Surf: Surface. It is a structure with Surf.tri = t x 3 matrix of triangle 
% indices, 1-based, t=#triangles and Surf.coord = 3 x nv matrix of 
% coordinates, nv=#vertices.
% maskvtx: Mask's vertices (1-based). Default [] (all vertices included).
%
% Output
% AdjM: Matrix of nv x cn whos rows contain the adjacent vertices of the
% vertex indicated by the row's number. If a given vertex has less adjacent
% vertices than cn then the remaining entries in its corresponding row in
% AdjM are zero.
% cn: Connectivity number (maximum number of adjacent vertices to any 
% given vertice along the surface). It is 6 for Freesurfer's meshes.
%
% Original Author: Jorge Luis Bernal Rusiel 
%
if nargin < 2
    maskvtx = [];
end;
surftri = Surf.tri;
[~,cn] = mode(double(surftri(:)));
if ~isempty(maskvtx)
    logictri = surftri;
    for i=1:3
        logictri(:,i) = ismember(surftri(:,i),maskvtx);
    end;
    surftri = surftri(sum(logictri,2) == 3,:);
end;
nv = size(Surf.coord,2);
AdjM = zeros(nv,cn);
last = ones(nv,1);
ntri = size(surftri,1);
for i=1:ntri
    tri = surftri(i,:);
    [AdjM(tri(1),:), last(tri(1))] = add(AdjM(tri(1),:),last(tri(1)),tri(2));
    [AdjM(tri(1),:), last(tri(1))] = add(AdjM(tri(1),:),last(tri(1)),tri(3));
    [AdjM(tri(2),:), last(tri(2))] = add(AdjM(tri(2),:),last(tri(2)),tri(1));
    [AdjM(tri(2),:), last(tri(2))] = add(AdjM(tri(2),:),last(tri(2)),tri(3));
    [AdjM(tri(3),:), last(tri(3))] = add(AdjM(tri(3),:),last(tri(3)),tri(1));
    [AdjM(tri(3),:), last(tri(3))] = add(AdjM(tri(3),:),last(tri(3)),tri(2));
end;
end

 

function [v,last] = add(v,last,el)
   if sum(v == el) == 0
       v(last) = el;
       last = last + 1;
   end;  
end
    
    
    
    
    
    
    
    
