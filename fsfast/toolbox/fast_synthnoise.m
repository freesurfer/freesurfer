function [n, F] = fast_synthnoise(nf,nc,acf)
% [n F] = fast_synthnoise(nf,nc,<acf>)
%
% $Id: fast_synthnoise.m,v 1.1 2003/04/11 23:08:56 greve Exp $


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
