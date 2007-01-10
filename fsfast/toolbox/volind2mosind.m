function [im, tszmos] = volind2mosind(iv, szvol, tszmos)
% [im tszmos] = volind2mosind(iv, szvol, tszmos)


%
% volind2mosind.m
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

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: [im tszmos] = volind2mosind(iv, szvol, <tszmos>)';
  error(msg);
end

szvol = szvol(1:3);
Nvr = szvol(1);
Nvc = szvol(2);
Nvs = szvol(3);

% Size of Mosaic measured in Tiles %
if(nargin == 2) tszmos = []; end
tszmos = defmossize(Nvs, tszmos);
Ntr = tszmos(1);
Ntc = tszmos(2);

[rv cv sv] = ind2sub(szvol,iv);
[rm cm] = volsub2mossub(rv,cv,sv,szvol,tszmos);

Nmr = Ntr*Nvr;
Nmc = Ntc*Nvc;
szmos = [Nmr Nmc];

im = sub2ind(szmos, rm, cm);

return;
