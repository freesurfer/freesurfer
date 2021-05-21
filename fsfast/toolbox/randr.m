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
%  rmean = wstd*sqrt(pi/2);
%  rstd  = wstd*sqrt(2-pi/2);
%  rstd  = rmean*sqrt(4/pi-1);
%  rmean = rstd/sqrt(4/pi-1);
%  wstd  = rmean/sqrt(pi/2);
%  wstd  = rstd/sqrt(2-pi/2);
%  fsnr  = rmean/rstd = 1/sqrt(4/pi-1) = 1.9137
%


%
% randr.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if(exist('rsize') ~= 1) rsize  = 1; end
if(exist('rmean') ~= 1) 
  wstd = 1;
  rmean = wstd*sqrt(pi/2);
else
  wstd = rmean/sqrt(pi/2);
end

r = abs(wstd*(randn(rsize) + i*randn(rsize)));

return;
