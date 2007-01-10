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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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
