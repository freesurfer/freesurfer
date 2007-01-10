function [mrgl, r, s, W] = fast_cvm_normposdef(m,nmax)
% [mrgl r] = fast_cvm_normposdef(m,nmax)
%
% Normalization and force to be positive by contructing a
% cholesky decomposition matrix from the auotcorrelation
% function


%
% fast_cvm_normposdef.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

mrgl = [];
if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: [mrgl r] = fast_cvm_normposdef(m, nmax)';
  qoe(msg);error(msg);
end

ntrs = size(m,1);

if(nargin ~= 2) nmax = ntrs; end

if(nmax > ntrs)
  msg = sprintf('nmax (%d) > ntrs (%d)\n',nmax,ntrs);
  qoe(msg);error(msg);
end

% estimate auto-correlation function
r = zeros(ntrs,1);
for n = 1:ntrs
  r(n) = mean(diag(m,n-1));
end
r = r/r(1);

R      = fft(r(1:nmax));
Rmag   = abs(R);
Rphase = angle(R);

Smag = 1./sqrt(Rmag);
Sphase = -Rphase;

S = Smag .* cos(Sphase) + sqrt(-1)*Smag .* sin(Sphase);
s0 = real(ifft(S));
s0 = s0/s0(1);

s = zeros(ntrs,1);
s(1:nmax) = s0;

tmp = zeros(ntrs,1);
tmp(1) = 1;

W = toeplitz(tmp,s);
mrgl = inv(W*W'); %'


return;
