function [n, F] = fast_synthnoise(nf,nc,acf)
% [n F] = fast_synthnoise(nf,nc,<acf>)
%
%


%
% fast_synthnoise.m
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


n = [];
F = [];

if(nargin < 2 | nargin > 3)
  fprintf('[n F] = fast_synthnoise(nf,nc,<acf>)\n');
  return;
end

if(exist('acf')~=1) acf = []; end
if(isempty(acf))
  n = randn(nf,nc);
  F = eye(nf);
  return;
end

nacf = length(acf);
if(nacf ~= nf)
  fprintf('ERROR: acf len = %d, nf = %d\n',nacf,nf);
  return;
end

S = toeplitz(acf);
mineig = min(eig(S));
if(mineig < 0)
  fprintf('ERROR: acf is not pos def\n');
  return;
end

F = chol(S);


n = F * randn(nf,nc);

% Correct for edge effects %
if(0)
  c = (F*ones(nf,1));
  n = n ./ repmat(c,[1 nc]);
end

return;
