function xyzlab = MRIseg2labelxyz(segmri,segid)
% xyzlab = MRIseg2labelxyz(segmri,<segid>)
%
% Computes the xyz for use in a label of the voxels in
% the segmentation volume.
%
% segmri - MRI struct. eg, segmri = MRIread('segvol');
% segid - id to find in the seg (default is 1)
%
% xyzlab - xyz in segvol's RAS space (defined by segmri.vox2ras0)
%


%
% MRIseg2labelxyz.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
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


xyzlab = [];
if(nargin < 2 | nargin > 3)
  fprintf('xyzlab = MRIseg2labelxyz(segmri,segid)\n');
  return;
end

if(~exist('segid','var')) segid = 1; end

% Get list of voxels with seg id
indlab = find(segmri.vol == segid);
nlab = length(indlab);

% Convert indices to row, col, slice
[r c s] = ind2sub(segmri.volsize,indlab);
crs = [c r s]' - 1 ; % 0-based
crs = [crs; ones(1,nlab)];

% Convert row, col, slice to XYZ
xyz1 = segmri.vox2ras0 * crs;
xyz = xyz1(1:3,:)';

return;
