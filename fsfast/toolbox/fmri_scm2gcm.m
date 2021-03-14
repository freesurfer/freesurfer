function gcm = fmri_scm2gcm(X,Nnnc,TR,tPreStim,delta,tau,alpha)
%
% gcm = fmri_scm2gcm(X,Nnnc,TR,tPreStim,delta,tau,<alpha>)
%
% Produces a Gamma Convolution Matrix from a stimulus
% convolution matrix (X) and parameters of the gamma
% function (delta, tau, alpha).  The gamma functions
% are interpreted as basis vectors. 
%
% If alpha is not specified, then it is set to 2.
%
% See also: fmri_hemodyn.m
%
%


%
% fmri_scm2gcm.m
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


if(nargin ~= 6 & nargin ~= 7)
  msg = 'gcm = fmri_scm2gcm(X,Nnnc,TR,tPreStim,delta,tau,<alpha>)';
  qoe(msg);error(msg);
end

if(nargin == 6) alpha = 2; end


[Ntp Nch Nr] = size(X);
Nh = Nch/Nnnc;
Ng = length(delta);

t = TR*[0:Nh-1] - tPreStim;
h = fmri_hemodyn(t,delta,tau,alpha);
h = h./(repmat(max(h),[Nh 1]));

h_all = zeros(Nch,Nnnc*Ng);
h0 = zeros(Nh,Nnnc*Ng);
h0(1:Nh,1:Ng) = h;
for c = 1:Nnnc,
  r1 = Nh*(c-1)+1;
  r2 = r1 + Nh - 1;
  h_all(r1:r2,:) = fmri_shiftcol(h0,Ng*(c-1));
end

gcm = zeros(Ntp,Nnnc*Ng,Nr);

for r = 1:Nr,
  gcm(:,:,r) = X(:,:,r)*h_all;
end

return;
