function mcvect = bmm_mcvect(v1,v2)
% mcvect = bmm_mcvect(v1,v2)
%
% Not tested yet!
%
% Construct a modified count vector (mcv). Used in the dependence
% model in binomial mixture modeling.
%
% v1 - binarization from method 1
% v2 - binarization from method 2
% In both v1 and v2, the size of the last dimension is the number
% of trials/replicants. The earler dimensions are spatial.
%
% mcvect will have four frames at each voxel indicating the number
% of times the following conditions were met
%   1 : neigther v1 nor v2 : vii = ~v1 & ~v2;
%   2 : not v1 but v2      : via = ~v1 &  v2;
%   3 : v1 but not v2      : vai =  v1 & ~v2; 
%   4 : v1 and v2          : vaa =  v1 &  v2;
% mcvect will have the same number of spatial dimensions as the
% input. It will have one more dimension with 4 elements as
% described above. The sum over this dimension will be constant 
% across all voxels and equal to number of replicants.
%
% Based on appendix of Genovese, et al, 1997. Estimating Test-Retest
% Reliability in Functional MR Imaging I: Statistical Methodology.
% MRM 38:497-507.
%


%
% bmm_mcvect.m
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


mcvect = [];
if(nargin ~= 2)
  fprintf('mcvect = bmm_mcvect(v1,v2)');
  return;
end

sz = size(v1);
sz2 = size(v2);

% Should check that they have the same size
% sz =? sz2

ndim = length(sz);
if(ndim == 2 & sz(2)==1) ndim = 1; end
nv = prod(sz);

% Number of trials/replicants (informational, not used)
M = sz(ndim);

vii = ~v1 & ~v2;
via = ~v1 &  v2;
vai =  v1 & ~v2;
vaa =  v1 &  v2;

mcvect = cat(ndim+1, vii, via, vai, vaa); % tmp

% Now sum over the replicant dimension
mcvect = squeeze(sum(mcvect,ndim));

return;











