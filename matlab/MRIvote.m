function [v nocc] = MRIvote(vol)
% [v nocc] = MRIvote(vol)
%
% For each voxel, selects most freqeuntly occurring value across 
% frames. nocc is the number of times it occurs.
% See also mri2.c::MRIvote()
% 
% $Id: MRIvote.m,v 1.2 2011/03/02 00:04:12 nicks Exp $

%
% MRIvote.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.2 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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


