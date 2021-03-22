function [err, racfexp] = fast_ar1_fiterr(rho,racf,R,w)
% err = fast_ar1_fiterr(rho,racf,R,<w>)
%
% error function for fitting an AR1 noise model while also taking
% into account the bias introduced by GLM fitting. This function
% can be used with fminsearch, which uses Nelder-Mead 
% unconstrained non-linear search. 'Unconstrained' may cause
% problems. 
%
% rho - AR1 parameter
% racf is actual residual ACF (nf-by-1). Should be created with
%       racf = fast_acorr(r,'unbiasedcoeff');
% R is residual forming matrix (nf-by-nf)
% w is the weight of each delay (nf-by-1). If unspecified,
%   no weighting is performed. In simulations, works best
%   without weighting. But those are simulations. Consider
%   weighting with w = 1./[1:nf]';
% 
% err = sum( abs(racf-racfexp).*w ); % L1 norm
% 
% Autocorrelation function of AR1+White Noise
%   acf(n) = rho^n, n=0:nf-1
%
% Example:
%  rho = racf(2); % Init
%  rhoopt = fminsearch('fast_ar1_fiterr',rho,[],racf,R);
%  nacfopt = rho.^[0:nf-1]';
%  [e racfexp] = fast_ar1_fiterr(rho,racfmn,R);
%  plot(1:nf,racf,1:nf,racfexp)
% 
%
%
%
% (c) Douglas N. Greve, 2004
%


%
% fast_ar1_fiterr.m
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

% Number of frames
nf = length(racf);

% Create ideal noise ACF based on AR1+W parameters
nacf = rho.^[0:nf-1]';

% Create the noise covariaance matrix
Sn = toeplitz(nacf);

% Create the expected residual covariance matrix
%Srexp = R*Sn*R;
% Extract the expected residual ACF
%racfexp = fast_cvm2acor(Srexp);

% But this is about 4X faster
racfexp = (R(1,:)*Sn*R)';
racfexp = racfexp/racfexp(1);

% Error is L1 difference between actual and expected
if(exist('w','var'))
  err = sum( abs(racf-racfexp).*w );
else
  err = sum( abs(racf-racfexp) );
end

return;

