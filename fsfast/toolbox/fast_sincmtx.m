function M = fast_sincmtx(len,delta,rownorm)
% M = fast_sincmtx(len,delta,<rownorm>)
% Create a sinc interpolation matrix. Given a vector y,
% yhat = M*y will be the sinc interpolation of y shifted
% backwards by delta, where delta is a fration of the 
% distance between y. If rownorm is set to non-zero, then
% each row of the matrix M is normalized so that the
% mean is 1.

M = [];

if( nargin ~= 2 & nargin ~= 3 )
  fprintf('M = fast_sincmtx(len,delta,<rownorm>)\n');
  return;
end
if(exist('rownorm') ~= 1) rownorm = []; end
if(isempty(rownorm)) rownorm = 0; end

M = zeros(len);
i = 1:len;
for n = 1:len
  r = sinc(i - n - delta);
  if(rownorm) r = r/sum(r); end
  M(n,:) = r;
end

return;

%-----------------------------------------------------------%
function y = sinc(x)

y      = ones(size(x));
inz    = find(x ~= 0);
y(inz) = sin(pi*x(inz)) ./ (pi*x(inz));

return;