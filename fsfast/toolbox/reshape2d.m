function m = reshape2d(x)
% m = reshape2d(x)
% Reshape into a 2d array where
% size(m,1) = size(x,1)

if(nargin ~= 1)
  msg = 'm = reshape2d(x)'
  qoe(msg);error(msg);
end

szx = size(x);
xdim = length(szx);

if(xdim < 2)
  msg = sprintf('Dimension of x is only %d',xdim);
  qoe(msg);error(msg);
end

m = reshape(x, [szx(1) prod(szx(2:xdim))]);
return
