function fn = fast_fnorm(f,fdim,demean)
% fn = fast_fnorm(f,<fdim>,<demean>)
%
% functional normalization - divides each voxel by sqrt 
% sum of the squares across frames. If demean flag is non-zero, 
% the mean is removed first.
%
% Assumes functional dimension is 1 unless set by fdim
%
%


%
% fast_fnorm.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

fn = [];

if(nargin < 1 | nargin > 3)
  fprintf('fn = fast_fnorm(f,<fdim>,<demean>)\n');
  return;
end

if(~exist('fdim','var')) fdim = []; end
if(isempty('fdim')) fdim = 1; end

if(~exist('demean','var')) demean = []; end
if(isempty('demean')) demean = 1; end

szf = size(f);
if(length(szf) < fdim)
  fprintf('ERROR: fdim = %d, but f only has %d dims\n',length(szf))
  return;
end

nf = szf(fdim);
repmatsz = ones(size(szf));
repmatsz(fdim) = nf;

fn = f;

if(demean)
  fmn = mean(f,fdim);
  fn = fn - repmat(fmn,repmatsz);
end

fss2 = sqrt(sum(fn.^2,fdim));
ind = find(fss2==0);
fss2(ind) = 1e10;

fn = fn./repmat(fss2,repmatsz);

return;
