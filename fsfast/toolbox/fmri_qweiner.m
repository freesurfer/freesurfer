function s = fmri_qweiner(Ur, Sr, Ux, Sx, y)
%
% s = fmri_qweiner(Ur, Sr, Ux, Sx, y)
%
% Computes s = R * inv(X) * y using (reduced) eigenvector
% representations of R and X.  When R is the signal
% covariance matrix and X is the signal+noise covariance
% matrix, and y is an observable, then s is the weiner 
% estimate of the signal.
%
% $Id: fmri_qweiner.m,v 1.1 2003/03/04 20:47:40 greve Exp $

if(nargin ~= 5)
  msg = 'USAGE: s = fmri_qweiner(Ur, Sr, Ux, Sx, y)';
  qoe(msg); error(msg);
end

[Nv Nr] = size(Ur);
Nt = size(y,2);

if(length(Sr) ~= Nr)
  msg = 'Ur and Sr have inconsistent dimensions';
  qoe(msg); error(msg);
end

if(size(Ux,1) ~= Nv)
  msg = 'Ur and Ux have inconsistent dimensions';
  qoe(msg); error(msg);
end

Nx = size(Ux,2);
if(length(Sx) ~= Nx)
  msg = 'Ux and Sx have inconsistent dimensions';
  qoe(msg); error(msg);
end

if(size(y,1) ~= Nv)
  msg = 'Ur and y have inconsistent dimensions';
  qoe(msg); error(msg);
end

s = Ur * ((diag(Sr) * Ur' * Ux * diag(1./Sx)) * (Ux' * y));

return;