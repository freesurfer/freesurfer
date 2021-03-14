function [n, F] = fast_synthnoise(nf,nc,acf)
% [n F] = fast_synthnoise(nf,nc,<acf>)
%
%


%
% fast_synthnoise.m
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

F = chol(S)';


n = F * randn(nf,nc);

% Correct for edge effects %
if(0)
  c = (F*ones(nf,1));
  n = n ./ repmat(c,[1 nc]);
end

return;
