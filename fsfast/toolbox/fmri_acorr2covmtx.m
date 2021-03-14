function CVM = fmri_acorr2covmtx(R,nCVM,Sided)
%
% CVM = fmri_acorr2covmtx(R)
% CVM = fmri_acorr2covmtx(R,nCVM)
% CVM = fmri_acorr2covmtx(R,nCVM,Sided)
%
% Arguments:
%   R - correlation function (normalized)
%   nCVM - number of rows/cols of the final covariance matrix.
%     if unspecified, nCVM = length(R);
%   Sided - 1 or 2. 
%      1 indicates one-sided (ie, zero lag to maxlag) 
%      2 indicates two-sided (ie, symetric, -maxlag to maxlag) (default)
%
%


%
% fmri_acorr2covmtx.m
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

if(nargin ~= 1 & nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: CVM = fmri_acorr2covmtx(R,<nCVM,<Sided>>)';
  qoe(msg);error(msg);
end

Sided = 1;
if(nargin == 1 | nargin == 2)  Sided = 2; end

[nR nRuns] = size(R);

if(Sided == 1)
  R1 = R;
else
  R1 = R([(nR-1)/2 + 1:nR],:);
end

[nR nRuns] = size(R1);

if(nargin == 1) 
  nCVM  = nR;
end

for r = 1:nRuns,
  if(nR < nCVM)
    R2 = cat(1,R1(:,r),zeros(nCVM-nR,1));
  else 
    R2 = R1(1:nCVM,r);
  end

  CVM(:,:,r) = toeplitz(R2);
end


return;

