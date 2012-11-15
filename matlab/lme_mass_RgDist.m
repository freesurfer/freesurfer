function Dist = lme_mass_RgDist(Coord,Rgvtx,AdjM)
% Dist = lme_mass_RgDist(Coord,Rgvtx,AdjM)
% 
% This function computes the distances among all vertices inside a region. 
% If AdjM is empty the Euclidean distance is computed, otherwise the
% distances are computed along the surface using the triangulation scheme.
%
% Input
% AdjM: Matrix whos rows contain the adjacent vertices along the surface of
% any vertex indicated by its row number (see AdjMtx).
% Coord: 3 x nv matrix of coordinates, nv=#vertices in the surface.
% Rgvtx: Vector whose entries are the indices of the vertices inside the
% region.
%
% Output
% Dist: Square symmetric matrix with distances among all locations inside
% the region.
%
% $Revision: 1.1 $  $Date: 2012/11/15 15:17:52 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: vinke $
%    $Date: 2012/11/15 15:17:52 $
%    $Revision: 1.1 $
%
if (nargin < 3)
    error('Too few inputs');
end;
nv = length(Rgvtx);
Coord = Coord(:,Rgvtx);
if nv == 1
    Dist = 0;
else
    Dist = zeros(nv,nv);
    aux = ones(1,nv);
    for i=1:nv
        c = Coord(:,i);
        Dist(i,:) = sqrt(sum((Coord - kron(aux,c)).^2));
    end;
    if ~isempty(AdjM)
        RgAdj = RgAdjMtx(AdjM,Rgvtx);
        Dist1 = dijkstra(RgAdj,Coord');
        Dist = (Dist + Dist1)/2;
    end;
end;