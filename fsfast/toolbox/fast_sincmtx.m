function M = fast_sincmtx(len,delta,rownorm,hw)
% M = fast_sincmtx(len,delta,<rownorm>,<hw>)
% Create a sinc interpolation matrix. Given a vector y,
% yhat = M*y will be the sinc interpolation of y shifted
% backwards by delta, where delta is a fration of the 
% distance between y. If rownorm is set to non-zero, then
% each row of the matrix M is normalized so that the
% mean is 1. hw is the width of the hann window.

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
  x = i - n - delta;
  if(exist('hw') == 1)
    r = sinc(x) .* hann(x,hw);
  else
    r = sinc(x);
  end
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

%-----------------------------------------------------------%
function h = hann(x,tau)

alpha = 0.5; % 0.54 for Hamming
h = zeros(size(x));
ind = find(abs(x) < tau);
h(ind) = alpha + (1-alpha)*cos(pi*x(ind)/tau);

return;