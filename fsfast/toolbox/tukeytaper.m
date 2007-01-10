function w = tukeytaper(nf,M)
% w = tukeytaper(nf,<M>)
%
% half-cosine taper from 1 to M.
% If M is not specified, M=nf
%
% Truncation is forced a n >= M.
%
%


%
% tukeytaper.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
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

w = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('w = tukeytaper(nf,<M>)\n');
  return;
end

if(exist('M')~=1) M = []; end
if(isempty(M))    M = nf; end

if(M > nf)
  fprintf('ERROR: M = %d must be less than nf = %d\n',M,nf);
  return;
end

tau = [0:nf-1]';

w = .5*(1+cos(pi*tau/M));
w(M:end) = 0;

return;
