function [err, racfexp] = fast_arnw_fiterr(p,racf,R,w,modeltype)
% [err, racfexp] = fast_arnw_fiterr(p,racf,<R>,<w>,<modeltype>)
%
% error function for fitting an ARN+White noise model while also
% taking into account the bias introduced by GLM fitting. This
% function can be used with fminsearch, which uses Nelder-Mead
% unconstrained non-linear search. 
%
% p = [phi1 phi2 ... phiN alpha]; 
%   phi1-phiN are the ARN parameters.
%     Eg, for an AR1, phi1 = 0.5 corresponds
%     to an acf = (0.5).^[1:nf]
%   alpha is the white noise parameter
%
% racf is actual residual ACF (nf-by-1). 
% R is residual forming matrix (nf-by-nf)
% w is the weight of each delay (nf-by-1). If unspecified,
%   no weighting is performed. In simulations, works best
%   without weighting. But those are simulations. Consider
%   weighting with w = 1./[1:nf]';
% modeltype = 0 or not specified means an ARN+W
% modeltype = 1 means an ARN with phi = p;
% 
% If R = [], then assumes R=I
%
% err = mean( abs(racf-racfexp).*w ); % L1 norm
% 
% See also: fast_arnw_acf
%
%
%
% (c) Douglas N. Greve, 2004
%


%
% fast_arnw_fiterr.m
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

err = [];
racfexp = [];

% Number of frames
nf = length(racf);

if(exist('modeltype','var') & modeltype == 1)
  % Pure ARN
  alpha = 0;
  phi = p;
else
  % Pure ARN+White
  alpha = p(end);
  phi = p(1:end-1);
end

% Penalize for alpha being out of range
if(alpha < 0 | alpha > 1)
  err = nf;
  return;
end

% Penalize for poles close to the unit circle
poles = roots([1 -phi(:)']);
if(~isempty(find(abs(poles) > .95)))
  err = nf;
  return;
end

% Theoretical autocor function
nacf = fast_arnw_acf(phi,nf,alpha);

if(exist('R','var') & ~isempty(R))
  % Create the noise covariaance matrix
  Sn = toeplitz(nacf);
  % Create the expected residual covariance matrix
  Srexp = R*Sn*R;
  % Extract the expected residual ACF
  racfexp = fast_cvm2acor(Srexp);
else
  racfexp = nacf;
end

% Error is L1 difference between actual and expected
if(exist('w','var') & ~isempty(w))
  err = mean( abs(racf-racfexp).*w );
else
  err = mean( abs(racf-racfexp) );
end

return;

