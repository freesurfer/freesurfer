function condidmap = par2condidmap(par,nullid)
%
% condidmap = par2condidmap(par,nullid)
%
% Extracts and sorts all unique condition numbers in paradigm.
%
% $Id: par2condidmap.m,v 1.1 2003/03/04 20:47:41 greve Exp $

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: condidmap = par2condidmap(par,<nullid>)'
  qoe(msg); error(msg);
end

condid = reshape1d(par(:,2,:,:));
condidmap = unique(condid);

if(nargin == 2)
  n = find(condidmap == nullid);
  if(isempty(n))
    condidmap = [nullid; reshape1d(condidmap)];
  else
    tmp = condidmap(1);
    condidmap(1) = nullid;
    condidmap(n) = tmp;
  end
end


return;
