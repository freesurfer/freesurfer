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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
