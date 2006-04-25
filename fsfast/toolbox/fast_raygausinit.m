function params = fast_raygausinit(y)
% params = fast_raygausinit(y)
% params = [alpha Rmu Gmu Gstd];
%
% $Id: fast_raygausinit.m,v 1.2 2006/04/25 05:04:32 greve Exp $

if(nargin ~= 1)
  fprintf('params = fast_raygausinit(y)\n');
  return;
end

%md = mode(y);
ysort = sort(y);
ny = length(y);
md = ysort(round(ny/2));


indray = find(y < md);
indgaus = find(y >= md);
alpha = .5;
Rmu = mean(y(indray));
Gmu = mean(y(indgaus));
Gstd = std(y(indgaus));

params = [alpha Rmu Gmu Gstd];

return;
