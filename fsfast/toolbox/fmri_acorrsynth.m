function R = fmri_acorrsynth(alpha, rho, nMaxDelay)
%
% R = fmri_acorrsynth(alpha, rho, nMaxDelay)
%
% Synthesizes an autocorrelation function based on the model
%    R(0) = 1
%    R(n) = (1-alpha)*rho^n,  0 < abs(n) <= nMaxDelay
%
% R(n) has 2*nMaxDelay elements.
%
% If alpha and rho have multiple elements, then multiple Rs
% computed and occupy a seperate column in R.
%
% See also fast_ar1w_acf.m
%
%


%
% fmri_acorrsynth.m
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

if(nargin ~= 3)
  msg = 'Usage: R = fmri_acorrsynth(alpha, rho, nMaxDelay)';
  qoe(msg);error(msg);
end

Sided = 2;
nR = max(nMaxDelay) + 1; 
nRuns = size(alpha,2);

if(Sided == 1) R=zeros(nR,nRuns);
else           R=zeros(2*(nR-1)+1,nRuns);
end

for r = 1:nRuns,

  if(Sided == 1)
    R(1,r)=1;
    R(2:nR,r) = (1-alpha(r))*(rho(r).^([1:nR-1]));
  else
    R(nR,r)=1;
    x = [1:nR-1]';
    z = (1-alpha(r))*(rho(r).^(x));
    R(x,r) = z(nR-x);
    R(nR+x,r) = z(x);
  end

end

return;

