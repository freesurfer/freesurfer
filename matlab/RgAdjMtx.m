function RgAdj = RgAdjMtx(AdjM,Rgvtx)
% RgAdj = RgAdjMtx(AdjM,Rgvtx)
% 
% This function finds a zero-one adjacency matrix for the vertices inside a
% given ROI.
%
% Input
% AdjM: Matrix whos rows contain the adjacent vertices along the original 
% surface of any vertex indicated by its row number (see AdjMtx).
% Rgvtx: Vector whose entries are the indices of the vertices in AdjM.
%
% Output
% RgAdj a nvROIxnvROI (nvROI=#vertices in the ROI) adjacency matrix, where 
% RgAdj(i,j) is nonzero (=1) if and only if an edge connects point i to 
% point j.
%
% Original Author: Jorge Luis Bernal Rusiel 
%
if nargin < 2 
    error('Too few inputs');   
end;
nv = length(Rgvtx);
RgAdj = zeros(nv,nv,'int8');
AdjM = AdjM(Rgvtx,:);
for i=1:nv
    [tf,loc] = ismember(AdjM(i,:),Rgvtx);   
    RgAdj(i,loc(tf == 1)) = 1;
end;
end


