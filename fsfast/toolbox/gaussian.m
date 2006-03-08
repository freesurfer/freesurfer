function y = gaussian(x,mean,stddev)
%
% Gaussian Distribution Function
%
% y = gaussian(x,mean,stddev)
%
% y = exp( -((x-mean).^2)/(2*(stddev.^2))) / (stddev * sqrt(2*pi));
% y = y/sum(y);

if(nargin ~= 3) 
  msg = 'USAGE: y = gaussian(x,mean,stddev)'
  error(msg);
end

y = exp( -((x-mean).^2)./(2*(stddev.^2))) ./ (stddev * sqrt(2*pi));

return
