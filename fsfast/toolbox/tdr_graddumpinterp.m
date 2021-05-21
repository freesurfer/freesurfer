function [gdi, tgdi] = tdr_graddumpinterp(gd,dt)
% [gdi tgdi] = tdr_graddumpinterp(gd,<dt>)
% sample gradient dump file at 2x rate
% dt defaults to 10 usec


%
% tdr_graddumpinterp.m
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









