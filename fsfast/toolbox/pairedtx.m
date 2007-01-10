function [X, C] = pairedtx(npairs)
% [X C] = pairedtx(npairs)
%
% Creates a design matrix for a paired test for npairs of data. It
% is assumed that each pair appears consecutively. The design
% matrix will have size 2*npairs by npairs+1. C can be used to test
% whether the paired difference is equal to zero. C is a
% vector of length npairs+1 with all the components equaling 0
% except the last which is 1.
%
%


%
% pairedtx.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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
  fprintf('[X C] = pairedtx(npairs)\n');
  return;
end

S = repmat([1 1]',[1 1 npairs]);
S = fast_blockdiag(S);
D = repmat([1 -1]',[npairs 1]);
X = [S D];

C = zeros(1,npairs+1);
C(end) = 1;

return;

