function [gdi, tgdi] = tdr_graddumpinterp(gd,dt)
% [gdi tgdi] = tdr_graddumpinterp(gd,<dt>)
% sample gradient dump file at 2x rate
% dt defaults to 10 usec


%
% tdr_graddumpinterp.m
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

gdi = [];

if(nargin < 1 || nargin > 2)
  fprintf('[gdi tgdi] = tdr_graddumpinterp(gd,<dt>)\n');
  return;
end

if(~exist('dt','var')) dt = []; end
if(isempty(dt)) dt = 10; end % 10 usec

ngd = size(gd,1);
tgd = dt*[0:ngd-1]';

ngdi = 2*ngd;
dti = dt/2;
tgdi = dti*[0:ngdi-1]';

gdi = interp1(tgd,gd,tgdi);

return









