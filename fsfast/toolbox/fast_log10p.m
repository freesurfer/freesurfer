function [log10p, nnz] = fast_log10p(p)
% [log10p nnz] = fast_log10p(p)
% Computes the log10 of p, keeping the sign
% and handling zeros appropriately
%
% log10p(indnz) = -sign(pnz) .* log10(abs(pnz));
%
% nnz is the number of p values not equal to zero.
%
% $Id: fast_log10p.m,v 1.1 2003/05/12 18:32:27 greve Exp $

log10p = [];
if(nargin ~= 1)
  fprintf('[log10p nnz] = fast_log10p(p)\n');
  return;
end

% Get the voxels that are zero 
indnz = find(p ~= 0);

log10p = zeros(size(p));
pnz = p(indnz);
log10p(indnz) = -sign(pnz) .* log10(abs(pnz));

return;