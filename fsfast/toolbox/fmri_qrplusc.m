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
% $Id: fmri_qrplusc.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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