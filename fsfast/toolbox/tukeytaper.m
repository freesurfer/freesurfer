function w = tukeytaper(nf,M)
% w = tukeytaper(nf,<M>)
%
% half-cosine taper from 1 to M.
% If M is not specified, M=nf
%
% Truncation is forced a n >= M.
%
%


%
% tukeytaper.m
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

w = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('w = tukeytaper(nf,<M>)\n');
  return;
end

if(exist('M')~=1) M = []; end
if(isempty(M))    M = nf; end

if(M > nf)
  fprintf('ERROR: M = %d must be less than nf = %d\n',M,nf);
  return;
end

tau = [0:nf-1]';

w = .5*(1+cos(pi*tau/M));
w(M:end) = 0;

return;
