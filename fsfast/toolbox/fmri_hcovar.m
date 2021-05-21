function Ch = fmri_hcovar(X,Cn)
%
% Ch = fmri_hcovar(X)
% Ch = fmri_hcovar(X,Cn)
%
% Computes the covariance matrix of the hemodynamic estimates.
% X is the stimulus convolution matrix. Cn is the normalized
% (ie, diag(Cn)=1) noise covariance matrix.  If unspecified,
% it is assumed that Cn=I (ie, white noise).
%
%


%
% fmri_hcovar.m
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
  msg = 'USAGE: Ch = fmri_hcovar(X, <Cn>)';
  qoe(msg);
  error(msg);
end

[nTP nTotEst nRuns] = size(X);

Ch = 0;

if(nargin == 1) % Assume white noise, ie Cn = I
  for r = 1:nRuns,
     Ch = Ch + X(:,:,r)' * X(:,:,r); %'
  end

else % Use Cn as noise covariance matrix

  if(size(Cn,1) ~= nTP | size(Cn,2) ~= nTP | ...
     (size(Cn,3) ~=1 & size(Cn,3) ~= nRuns))
    msg = 'X and Cn dimensions are inconsistent';
    qoe(msg);
    error(msg);
  end

  if(size(Cn,3) == 1) % only one matrix specified, use for all runs
    iCn = inv(Cn);
    for r = 1:nRuns,
       Ch = Ch + X(:,:,r)' * iCn * X(:,:,r); %'
    end
  else  % use different Cn for each run
    for r = 1:nRuns,
       Ch = Ch + X(:,:,r)' * inv(Cn(:,:,r)) * X(:,:,r); %'
    end  
  end

end

Ch = inv(Ch);

return
