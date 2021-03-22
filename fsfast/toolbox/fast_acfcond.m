function [cnd, mineig, S] = fast_acfcond(acf,taper)
% [cnd mineig S] = fast_acfcond(acf,<taper>)


%
% fast_acfcond.m
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

cnd = [];
mineig = [];
S = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('[cnd mineig S] = fast_acfcond(acf,<taper>)\n');
  return;
end

if(exist('taper') ~= 1) taper = []; end
if(~isempty(taper))
  acf = acf .* taper;
end

S = toeplitz(acf);
mineig = min(eig(S));
cnd = cond(S);

return;







