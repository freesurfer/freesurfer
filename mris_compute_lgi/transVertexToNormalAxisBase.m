function result = transVertexToNormalAxisBase(mesh_outer, iV, coord)

% Originally written by Lucas Tamarit, tamarit@cisa.unige.ch
% Adapted to the local GI whole computation on 2007/15/11
% 
% create a new referential based on the normal of the outer surface at the
% center of the outer region of interest. In order to give the distance of
% a given point (coord) the main axis of this referential. 


B.origin = mesh_outer.vertices(iV,:); % set the origin of the referential 
n = mesh_outer.vertexNormal(iV,:); % use the normal to the outer smoothed surface at that point 

B.axe = n ./ norm(n);
B.u = getOrthogonalVector(B.axe);
B.v = cross(B.axe,B.u);
B.v = B.v ./ norm(B.v);

vt = coord - B.origin;
vt = [B.u; B.v; B.axe]*vt'; 

result = vt;
