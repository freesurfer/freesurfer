function [p,p2] = fast_trf_sup_surf(t,dof,fwhm,area)
% [p p2] = fast_trf_sup_surf(t,dof,fwhm,area)
%
% Based on Moo Chung, et al, Cortical thickness analysis in autism
% with heat kernel smoothing. NeuroImage, 25, 1256-1265. 2005.
%
%


%
% fast_trf_sup_surf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

p = [];
if(nargin ~= 4)
  fprintf('[p p2] = fast_trf_sup_surf(t,dof,fwhm,area)\n');
  return;
end

p0 = tTest(dof,t)/2;

gstd = fwhm/sqrt(log(256));
lambda = 1/(2*gstd*gstd);

tmp = t.*(1 + (t.^2)/dof).^(-(dof-1)/2);

f = (lambda*gamma((dof+1)/2)) /...
    ( ((2*pi)^(3/2)) * sqrt(dof/2) * gamma(dof/2) );
p2 = f*tmp;

p = 2*p0 + area*p2/2;
%p = p0;

return;















