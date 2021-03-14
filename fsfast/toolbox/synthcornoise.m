function [n, CVM, W2, R1] = synthcornoise(ntrs,ncols,alpha,rho)
% [n CVM W2 R1] = synthcornoise(ntrs,ncols,alpha,rho)
%
% Synthesizes a functional slice of correlated noise using the
% alpha-rho model (see fmri_acorrsynth()) to construct the 
% autocorrelation function.  Synthesizes the noise using 2*ntrs
% points and then takes the middle half so as to avoid edge effects.
% The output data matrix n will have dimension ntrs X ncols.
%
% The asymptotoic covariance matrix is CVM. This can be compared to
% the sampled covariance matrix CVMhat = n*n'/nvoxs; %'
%
% W2 is the square of the whitening matrix that can be used to remove 
% the correlation from n, W2 = W * W'; W = chol(W2);           %' 
% 
% alpha = 1 or rho = 0 will create white noise


%
% synthcornoise.m
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

if(nargin ~= 4)
  msg = 'USAGE: [n, CVM, W2] = synthcornoise(ntrs,ncols,alpha,rho)';
  qoe(msg); error(msg);
end

ntrs2 = 2*ntrs;

% Two-sided autocorrelation function %
R = fmri_acorrsynth(alpha, rho, ntrs2);

% One-sided autocorrelation function %
R1 = R(ntrs2+1:2*ntrs2); % size is changed below

% Assymptotic Covariance Matrix %
CVM = toeplitz(R1);

% Filter to achieve CVM %
CVMFilter = chol(CVM);

% Generate White Noise %
y = randn(ntrs2,ncols);

% Filter to create Correlated Noise %
n = CVMFilter*y;

% Square of the Whitening Matrix
W2 = inv(CVM);

% Only take the middle ntrs %
nthird = floor(ntrs/3);
nkeep = nthird:nthird+ntrs-1;
CVM = CVM(nkeep,nkeep);
W2  = W2(nkeep,nkeep);
n   = n(nkeep,:);
R1 = R1(1:ntrs);

return
