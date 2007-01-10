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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
