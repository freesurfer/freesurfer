function Ce = fmri_ecovar(e,VoxId)
%
% Ce = fmri_ecovar(e)
% Ce = fmri_ecovar(e,VoxId)
%
% Computes covariance matrix of e.  e is typically the
% residual temporal error at each voxel.
%
% e has dimension nRows x nCols x nTP x nRuns.
% VoxId - 1D list of subset of voxels to use.
% Ce will have dimension nTP x nTP x nRuns.
%
% $Id: fmri_ecovar.m,v 1.1 2003/03/04 20:47:39 greve Exp $

[nRows nCols nTP nRuns] = size(e);

e = reshape(e, [nRows*nCols nTP nRuns]);
e = permute(e, [2 1 3]);

Ce = zeros(nTP,nTP,nRuns);

for r = 1:nRuns,

  if(nargin ==2) Ce(:,:,r) = e(:,VoxId,r) * e(:,VoxId,r)'; %'
  else           Ce(:,:,r) = e(:,:,r) * e(:,:,r)';
  end

end


if(nargin==2) nV = length(VoxId); 
else          nV = nRows * nCols;
end

Ce = Ce/nV;

return;
