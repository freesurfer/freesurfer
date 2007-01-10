function [Urc, Src] = fmri_rplusc(Ur,Sr,Uc,Sc,cutofftype,cutoffthresh)
%
% [Urc, Src] = fmri_rplusc(Ur,Sr,Uc,Sc,cutofftype,cutoffthresh)
%
% Computes the (reduced) eigenvectors/values for R+C given the
% (reduced) eigenvectors/values for R and C individually.  R and
% C are symetrical, positive-def matricies.
%
% Algorithm: Let Q = R+C, for which [A X A'] = svd(Q).  Let Ahat
% and Xhat be the reduced eigencomponents Q.  Then, Urc = Ahat,
% and Src = Xhat.  However, for R and C with a large number of rows,
% computing Q and svd(Q) is not feasible.  Ahat and Xhat can still
% be determined in a round-about way by first constructing the matrix
% P = [Ur*diag(Sr) Uc*diag(Sc)].  This is a matrix whose columns are
% the (reduced) eigenvectors of R and C (weighted by the corresponding
% eigenvalues).  Let Up, Sp be the reduced eigencomponents of P, then
% Ahat = P*Up*diag(Sp), and Xhat = Sp;
% 
%


%
% fmri_qrplusc.m
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

if(nargin ~= 6)
  msg = 'USAGE: [Urc, Src] = fmri_rplusc(Ur,Sr,Uc,Sc,cutofftype,cutoffthresh)';
  qoe(msg);error(msg);
end

if(size(Ur,2) ~= length(Sr))
  msg = 'Ur and Sr have incompatable number of items';
  qoe(msg); error(msg);
end

if(size(Uc,2) ~= length(Sc))
  msg = 'Uc and Sc have incompatable number of items';
  qoe(msg); error(msg);
end

if(size(Ur,1) ~= size(Uc,1))
  msg = 'Ur and Uc have incompatable number of rows';
  qoe(msg); error(msg);
end

% Create a new matrix of the concatenated eigenvectors
P = [Ur*diag(Sr) Uc*diag(Sc)];

% Do a quick svd of this matrix %
[Urc Src] = fmri_qsvd(P,cutofftype,cutoffthresh);

% Compute eigencomponents for R+C %
%tmp = P*Up*diag(1./Sp);
%Urc = fmri_norm(tmp,2);
%Src = Sp;

return;
