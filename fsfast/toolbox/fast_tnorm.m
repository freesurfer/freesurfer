function ftnorm = fast_tnorm(f)
% imgtnorm = fast_tnorm(f)
% 
% Temporal Normalization f: N1xN2x...Nt
%


%
% fast_tnorm.m
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

if(nargin ~= 1)
  msg = 'imgtnorm = fast_tnorm(f)';
  qoe(msg);error(msg);
end

szf = size(f);
ndimf = length(szf);

aszf = szf(1:ndimf-1);
nvox = prod(aszf);
nt   = szf(ndimf);

f = reshape(f, [nvox nt]);
fmn = mean(f,2);
fstd = std(f,[],2);

ftnorm = (f - repmat(fmn,[1 nt])) ./ repmat(fstd,[1 nt]) ;

ftnorm = reshape(ftnorm,szf);

return;
