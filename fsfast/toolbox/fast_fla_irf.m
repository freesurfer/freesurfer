function [irf, irftaxis] = fast_fla_irf(fla,beta)
% [irf irftaxis] = fast_fla_irf(fla,beta)
%
% Returns the Impulse Response Function (IRF) for the current
% fx (ie, nthfx) and the time axis for the IRF. Uses beta
% to weight the components of the assumption matrix. The
% nthfx must be an ERM.
%
%


%
% fast_fla_irf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

irf = [];
irftaxis = [];

if(nargin ~= 2)
  fprintf('[irf irftaxis] = fast_fla_irf(fla,beta)\n');
  return;
end

if(fla.nthfx > length(fla.fxlist))
  fprintf('fla_irf: nthfx=%d too big\n',fla.nthfx);
  return;
end

if(~fast_fxcfg('iserm',fla))
  fprintf('fla_irf: nthfx=%d is not an ERM\n',fla.nthfx);
  return;
end

fx = fla.fxlist(fla.nthfx).fx;
A = fast_fxcfg('irfmatrix',fla);

if(strcmp(fx.fxtype,'fixed'))
  ind = fx.regind{1};
else
  ind = fx.regind{fla.nthrun};
end

irf = A*beta(ind,:);
irftaxis = fast_fxcfg('irftaxis',fla);

return;
