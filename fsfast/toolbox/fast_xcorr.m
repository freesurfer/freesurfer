function rho = fast_xcorr(x,y,demean)
% rho = fast_xcorr(x,y,demean)
% Default is to demean. Set demean=0 to turn off.
%   rho = sum(x.*y)./sqrt(sum(x.^2) .* sum(y.^2));

rho = [];
if(nargin < 2 | nargin > 3)
  fprintf('rho = fast_xcorr(x,y,demean)\n');
  return;
end

if(nargin == 2) demean = 1; end

[nt nv] = size(x);

if(demean)
  x = x - repmat(mean(x),[nt 1]);
  y = y - repmat(mean(y),[nt 1]);
end

rho = sum(x.*y)./sqrt(sum(x.^2) .* sum(y.^2));

return;




