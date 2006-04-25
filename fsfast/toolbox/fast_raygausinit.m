function params = fast_raygausinit(y)
% params = fast_raygausinit(y)
% params = [alpha Rmu Gmu Gstd];
%
% $Id: fast_raygausinit.m,v 1.3 2006/04/25 05:19:30 greve Exp $

if(nargin ~= 1)
  fprintf('params = fast_raygausinit(y)\n');
  return;
end

if(0)
%md = mode(y);
ysort = sort(y);
ny = length(y);
md = ysort(round(ny/2));
indray = find(y < md);
indgaus = find(y >= md);
end

ymn = mean(y);
indray  = find(y < ymn);
indgaus = find(y > ymn);

alpha = .25;
Rmu = mean(y(indray));
Gmu = mean(y(indgaus));
Gstd = std(y(indgaus));

params = [alpha Rmu Gmu Gstd];

return;
