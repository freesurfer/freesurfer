function [X C] = fast_anova_rm_oneway(Ny,Nr)
% [X C] = fast_anova_rm_oneway(Ny,Nr)
%
% Design and contrast matrices for a One-way Repeated Measures
% ANOVA with Nr replicants 
%
% Ny is the number of data points which must be Ns*Nr, where Ns is the
% number of subjects.  The data must be arranged in y such that the
% first Nr rows are the replicants of subject 1, the next Nr rows are
% the replicants of subject 2, etc. When Nr=1, this reduces to an
% unsigned paired t.
%
% X is the design matrix
% C is the contrast matrix
%
%    [beta, rvar, vdof, r] = fast_glmfit(y,X);
%    [F, Fsig, con] = fast_fratio(beta,X,rvar,C);
%
% $Id: fast_anova_rm_oneway.m,v 1.1 2009/05/26 22:03:02 greve Exp $

X = [];
C = [];
if(nargin ~= 2)
  fprintf('[X C] = fast_anova_rm_oneway(Ny,Nr)\n');
  return;
end

if(rem(Ny,Nr) ~= 0) 
  fprintf('ERROR: number of rows in y (%d) is not integer\n',Ny);
  fprintf('multiple of the number of repeats (%d)\n',Nr);
  return;
end

% Number of subjects
Ns = Ny/Nr;
  
% Create a matrix to model the subject means
Xsubj = fast_blockdiag(ones(Nr,1,Ns));

% Create a matrix to model the difference between replicants
XrepDiff = [];
for nthRep = 1:Nr-1
  M = zeros(Nr,1);
  M(1) = 1;
  M(nthRep+1) = -1;
  XM = repmat(M,[Ns 1]);
  XrepDiff = [XrepDiff XM];
end

% Final design matrix
X = [Xsubj XrepDiff];

% Contrast matrix tests for a difference between replicants
C = [zeros(Nr-1,Ns) eye(Nr-1)];

return;

