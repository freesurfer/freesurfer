function [r, rmean, wstd] = randr(rsize,rmean)
%
% [r rmean wstd] = randr(rsize,rmean)
%
% Generates matrix of given dimension whose elements have a rayleigh
% distribution. The expected value is rmean. If rmean is not
% specified, then rmean = sqrt(pi/2), which makes wstd=1. wstd
% is the standard devaition of the white noise used to generate
% the data:
%   r = abs(wstd*(randn(rsize) + i*randn(rsize)));
%
% Note:
%  rmean = rstd/sqrt(4/pi-1);
%  rmean = wstd*sqrt(pi/2);
%  rstd  = rmean*sqrt(4/pi-1);
%  wstd  = rmean/sqrt(pi/2);
%  wstd  = rstd/sqrt(2-pi/2);
%
% $Id: randr.m,v 1.1 2005/03/19 00:23:31 greve Exp $

if(exist('rsize') ~= 1) rsize  = 1; end
if(exist('rmean') ~= 1) 
  wstd = 1;
  rmean = wstd*sqrt(pi/2);
else
  wstd = rmean/sqrt(pi/2);
end

r = abs(wstd*(randn(rsize) + i*randn(rsize)));

return;
