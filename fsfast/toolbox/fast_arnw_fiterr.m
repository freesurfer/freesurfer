function [err, racfexp] = fast_arnw_fiterr(p,racf,R,w)
% [err, racfexp] = fast_arnw_fiterr(p,racf,R,w)
%
% error function for fitting an ARN+White noise model while also
% taking into account the bias introduced by GLM fitting. This
% function can be used with fminsearch, which uses Nelder-Mead
% unconstrained non-linear search. 
%
% p = [alpha phi1 phi2 ... phiN]; 
%   phi1-phiN are the ARN parameters consistent with aryule.
%   Note: these are negative relative to the "standard" 
%   interpretation. Eg, for an AR1, phi1 = -0.5 corresponds
%   to an acf = (0.5).^[1:nf]
%
% racf is actual residual ACF (nf-by-1). 
% R is residual forming matrix (nf-by-nf)
% w is the weight of each delay (nf-by-1). If unspecified,
%   no weighting is performed. In simulations, works best
%   without weighting. But those are simulations. Consider
%   weighting with w = 1./[1:nf]';
% 
% err = sum( abs(racf-racfexp).*w ); % L1 norm
% 
% See also: fast_arn2_acf
%
% $Id: fast_arnw_fiterr.m,v 1.1 2004/05/24 00:09:29 greve Exp $
%
% (c) Douglas N. Greve, 2004
%

% Number of frames
nf = length(racf);

% Penalize for alpha being out of range
alpha = p(1);
if(alpha < 0 | alpha > 1)
  err = nf;
  return;
end

% Penalize for poles close to the unit circle
phi = p(2:end);
phi = [1 phi(:)'];
poles = roots(phi);
if(~isempty(find(abs(poles) > .95)))
  err = nf;
  return;
end

% Theoretical autocor function
nacf = fast_arnw_acf(p,nf);

% Create the noise covariaance matrix
Sn = toeplitz(nacf);

% Create the expected residual covariance matrix
Srexp = R*Sn*R;

% Extract the expected residual ACF
racfexp = fast_cvm2acor(Srexp);

% Error is L1 difference between actual and expected
if(exist('w','var'))
  err = sum( abs(racf-racfexp).*w );
else
  err = sum( abs(racf-racfexp) );
end

return;

