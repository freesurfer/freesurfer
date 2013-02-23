function Dist = lme_mass_RgDist(Surf,maskvtx,vtxInd,Dtype)
% Dist = lme_mass_RgDist(Surf,maskvtx,vtxInd,Dtype)
% 
% This function computes the distances among all vertices inside a region. 
%
% Input
% Surf: Triangular surface. It is a structure with Surf.tri = t x 3 matrix 
% of triangle indices, 1-based, t=#triangles and Surf.coord = 3 x nv matrix 
% of coordinates, nv=#vertices.
% maskvtx: Mask's vertices (1-based). 
% vtxInd: Vector whose entries are the indices in maskvtx of the vertices 
% that belong to the region.
% Dtype: Type of distances to be computed among the surface nodes. It is
% 'euc' for Euclidean or 'surf' for computing distances along the surface 
% using the triangulation squeme. Default 'euc'.
%
% Output
% Dist: Square symmetric matrix with distances among all locations inside
% the region.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:10 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:10 $
%    $Revision: 1.1.2.2 $
%
if (nargin < 3)
    error('Too few inputs');
elseif nargin < 4
    Dtype = 'euc';
end;
nv = length(vtxInd);
if nv == 1
    Dist = 0;
elseif strcmpi(Dtype,'euc')
    Coord = Surf.coord(:,maskvtx(vtxInd));
    Dist = zeros(nv,nv);
    aux = ones(1,nv);
    for i=1:nv
        c = Coord(:,i);
        Dist(i,:) = sqrt(sum((Coord - kron(aux,c)).^2));
    end;
elseif strcmpi(Dtype,'surf')
        global geodesic_library;                
        geodesic_library = 'libgeodesic'; 
        logictri = ismember(Surf.tri,maskvtx(vtxInd));
        patch_tri = Surf.tri(logictri(:,1) | logictri(:,2) | logictri(:,3),:);
        [patch_vtx,~,loc] = unique(patch_tri);
        patch_tri = reshape(loc,size(patch_tri,1),3);
        [~,Rgvtx_loc] = ismember(maskvtx(vtxInd),patch_vtx);
        patch_coord = Surf.coord(:,patch_vtx)';
        n = size(patch_coord,1);    
        mesh = geodesic_new_mesh(patch_coord,patch_tri);         %initilize new mesh
        algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm
        Dist = zeros(n,nv);
        for i=1:nv
            source_points = {geodesic_create_surface_point('vertex',Rgvtx_loc(i),patch_coord(Rgvtx_loc(i),:))};
            geodesic_propagate(algorithm,source_points);   %propagation stage of the algorithm (the most time-consuming)
           [~,Dist(:,i)] = geodesic_distance_and_source(algorithm);     %find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1
        end;
        Dist = Dist(Rgvtx_loc,:);
        geodesic_delete;
else 
    error('Valid values for Dtype are ''euc'' or ''surf''');
end;


