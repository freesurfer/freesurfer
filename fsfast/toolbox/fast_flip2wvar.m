function [wvar, ivar] = fast_flip2wvar(mn1,mn2,var1,var2)
% [wvar, ivar] = fast_flip2wvar(mn1,mn2,var1,var2)
%
% Compute white and instability noise given two measures
% with different proportions.
%
% $Id: fast_flip2wvar.m,v 1.1 2007/05/03 00:01:45 greve Exp $

%
% fast_flip2var
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/05/03 00:01:45 $
%    $Revision: 1.1 $
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

wvar = [];
ivar = [];

if(nargin ~= 4)
  fprintf('[wvar, ivar] = fast_flip2wvar(mn1,mn2,var1,var2)\n');
  return;
end

% This is the factor by which the instability should change due to
% changing the acq parameters:
v = (mn1./mn2).^2;

wvar = (v*var2 - var1)./(v-1);
ivar = var1 - wvar;

return;





