function newpar = remap_par(par,condidmap)
%
% newpar = remap_par(par,condidmap)
%
%


%
% remap_par.m
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
