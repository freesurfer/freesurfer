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

fn = [];

if(nargin < 1 | nargin > 3)
  fprintf('fn = fast_fnorm(f,<fdim>,<demean>)\n');
  return;
end

if(~exist('fdim','var')) fdim = []; end
if(isempty(fdim)) fdim = 1; end

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
