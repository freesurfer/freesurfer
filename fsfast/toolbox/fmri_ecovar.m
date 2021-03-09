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
%


%
% fmri_ecovar.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

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
