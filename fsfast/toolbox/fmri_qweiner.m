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
%


%
% fmri_qweiner.m
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
