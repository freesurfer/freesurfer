function gcm = fmri_scm2gcm(X,Nnnc,TR,tPreStim,delta,tau)
%
% gcm = fmri_scm2gcm(X,delta,tau)
%
% Produces a Gamma Convolution Matrix from a stimulus
% convolution matrix (X) and parameters of the gamma
% function (delta, tau).  The gamma functions
% are interpreted as basis vectors.
%
% $Id: fmri_scm2gcm.m,v 1.1 2003/03/04 20:47:40 greve Exp $


[Ntp Nch Nr] = size(X);
Nh = Nch/Nnnc;
Ng = length(delta);

t = TR*[0:Nh-1] - tPreStim;
h = fmri_hemodyn(t,delta,tau);
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
