function newpar = remap_par(par,condidmap)
%
% newpar = remap_par(par,condidmap)
%
% $Id: remap_par.m,v 1.1 2003/03/04 20:47:41 greve Exp $

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
