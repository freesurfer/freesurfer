function newpar = remap_par(par,condidmap)
%
% newpar = remap_par(par,condidmap)
%
%


%
% remap_par.m
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

if(nargin ~= 2)
  msg = 'USAGE: newpar = remap_par(par,condidmap)'
  qoe(msg);error(msg);
end

nconds = length(condidmap);

%tpar = reshape1d(par(:,1,:,:));
cpar = reshape1d(par(:,2,:,:));
cnewpar = zeros(size(cpar));

for cond = 0:nconds-1
  condid = condidmap(cond+1);
  ind = find(cpar == condid);
  cnewpar(ind) = cond;
end

npars = prod(size(par))/(2*size(par,1));

%tpar2 = reshape(tpar, [size(par,1) 1 npars]);
cnewpar = reshape(cnewpar, [size(par,1) 1 npars]);
newpar = reshape(par, [size(par,1) 2 npars]);
newpar(:,2,:) = cnewpar;
newpar = reshape(newpar, [size(par)]);

return;
