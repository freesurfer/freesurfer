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
%


%
% randr.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
