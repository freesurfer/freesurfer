function [alpha, rho] = fmri_acorrfit(R,nFit)
%
% [alpha, rho] = fmri_acorrfit(R)
% [alpha, rho] = fmri_acorrfit(R,nFit)
%
% Fits the autocorrelation function R to the model:
%    R(n) = 1,                n = 0
%    R(n) = (1-alpha)*rho^n,  n != 0
% and returns alpha and rho.  Components of R <= 0 are excluded
% from the fit.
%
% R should be two-sided column vector, normalized with zero lag = 1.
% If R has multiple columns, an alpha and rho are fit for each
% column.  If nFit is supplied, then only the first nFit components
% of R are fit.  nFit must be a scalar.
%
%


%
% fmri_acorrfit.m
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE [alpha rho] = fmri_acorrfit(R,nFit)';
  qoe(msg); error(msg);
end

if(length(nFit) ~= 1)
  msg = sprintf('nFit must be a scalar (%d)',length(nFit));
  qoe(msg); error(msg);
end

[nR nRuns] = size(R);

nR1 = (nR+1)/2;
R1 = R(nR1:nR,:);


for r = 1:nRuns,

  if(nargin == 1) nFit = nR1; end

  if(nFit > nR1) nFit = nR1; end
  X = [ ones(nFit-1,1) [1:nFit-1]' ]; %'
  R2 = R1([2:nFit],r);

  %% Exclude points where R <= 0 %%
  RMask = (R2 > 0);
  iRMask = find(R2>0);
  if(isempty(~iRMask) | R2(1)<0)
    alpha(r) = 1;
    rho(r)   = 0;
  else
    iRMask = iRMask(2:length(iRMask));
    NN = [1:length(R2)]';
    rho(r) = mean(exp(log(R2(iRMask)./R2(1))./NN(iRMask)));
    alpha(r) = 1.0 - R2(1)/rho(r);
  end

end

return;

