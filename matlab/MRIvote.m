function [v nocc] = MRIvote(vol)
% [v nocc] = MRIvote(vol)
%
% For each voxel, selects most freqeuntly occurring value across 
% frames. nocc is the number of times it occurs.
% See also mri2.c::MRIvote()
% 
% $Id: MRIvote.m,v 1.1 2007/08/01 21:27:43 greve Exp $

%
% MRIvote.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/08/01 21:27:43 $
%    $Revision: 1.1 $
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

volsize = size(vol);

mat = fast_vol2mat(vol);
[nf nv] = size(mat);
u = unique(mat);
nu = length(u);

vu = zeros(nu,nv);
for nthu = 1:nu
  vu(nthu,:) = sum(mat==u(nthu));
end

[nocc uind] = max(vu);
v = u(uind)';

v    = fast_mat2vol(v,volsize);
nocc = fast_mat2vol(nocc,volsize);

return;


