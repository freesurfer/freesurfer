function [thetahat, beta, rvar, nover] = tdr_phdist_est(Ref,Mov,rthresh,order)
% [thetahat beta rvar nover] = tdr_phdist_est(Ref,Mov,rthresh,order)
%
% Estimate the phase distortion between two waveforms.
%
% Ref - DFT of reference waveform (column vector). May have mult cols.
% Mov - DFT of movable waveform (column vector). May have mult cols.
% rthresh - relative threshold used for segmentation. The actual
%   threshold is rthresh times the mean of abs of Ref and Mov. 
%   Estimation is only performed on suprathreshold frames. 
%   Default is 2. All columns use the same segmentation.
% order - use a polynomial of given order for estimation. 
%   Default is 2.
%
% thetahat - model fit of the phase distortion.
% beta - coefficients of polynomial regressors
% rvar - residual error variance of fit
% nover - number of columns exceeding threshold
%
% RefHat = Mov .* exp(+i*thetahat);
%
% $Id: tdr_phdist_est.m,v 1.1 2003/11/06 19:47:02 greve Exp $

beta  = [];
rvar  = [];
nover = [];
thetahat = [];

if(nargin < 2 | nargin > 4)
  fprintf('[thetahat beta rvar nover] = tdr_phdist_est(Ref,Mov,rthresh,order)\n');
  return;
end

if(exist('rthresh') ~= 1) rthresh = []; end 
if(isempty(rthresh)) rthresh = 2; end 
if(exist('order') ~= 1) order = []; end 
if(isempty(order)) order = 2; end 
  
if(max(abs(size(Ref)-size(Mov))) ~= 0)
  fprintf('ERROR: size mismatch between Ref and Mov\n');
  return;
end
[nf nv] = size(Ref);

RefMn = mean(abs(Ref),2);
MovMn = mean(abs(Mov),2);
RefMovMn = (RefMn+MovMn)/2;

% Compute the absolute segmentation threshold %
mn = mean(RefMovMn);
athresh = rthresh * mn; 

% Segment - all columns get the same segmentation %
indover = find(RefMovMn > athresh);
nover = length(indover);
if(nover == 0) 
  fprintf('ERROR: no voxels over threshold\n');
  return;
end

% Compute the phase distortion %
theta = unwrap(angle(Ref)-angle(Mov)); 
thetaover = theta(indover,:);

% Create a polynomial design matrix %
X = fast_polytrendmtx(1,nf,1,order);
X = fast_polytrendmtx(1,nf,1,3);
X = X(:,[1 2 4]);

% Use only the suprathreshold frames for estimation %
Xover = X(indover,:);

% Perform the estimation on suprathreshold frames %
[beta rvar vdof] = fast_glmfit(thetaover,Xover);

% Compute the estimate of the phase distortion across all frames
thetahat = X*beta;

return;
