function [M, indnz, indz] = fast_rmzerocols(M0)
% [M indnz indz] = fast_rmzerocols(M0)
%
% Removes columns from M0 that are all zero
% indnz are the indicies of the non-zero columsn
% indz  are the indicies of the zero columsn


%
% fast_rmzerocols.m
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

if(nargin ~= 1)
  fprintf('[M indnz indz] = fast_rmzerocols(M0)\n');
  return;
end

ncols = size(M0,2);

nth = 1;
M = [];
indnz = [];
indz = [];
for c = 1:ncols
  indnzc = find(M0(:,c) ~= 0);
  if(~isempty(indnzc))
    M = [M M0(:,c)];
    indnz = [indnz c];
    nth = nth + 1;
  else
    indz = [indz c];
  end
end

return;
