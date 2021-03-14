function [log10p, nnz] = fast_log10p(p)
% [log10p nnz] = fast_log10p(p)
% Computes the log10 of p, keeping the sign
% and handling zeros appropriately
%
% log10p(indnz) = -sign(pnz) .* log10(abs(pnz));
%
% nnz is the number of p values not equal to zero.
%
%


%
% fast_log10p.m
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
