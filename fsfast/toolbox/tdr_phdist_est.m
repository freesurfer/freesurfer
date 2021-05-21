function [thetahat, theta, Fsig, beta] = tdr_phdist_est(Ref,Mov,rthresh,order)
% [thetahat, theta, Fsig, beta] = tdr_phdist_est(Ref,Mov,rthresh,order)
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
% theta  - phase distortion.
% Fsig - -log10(p) of the fit of the slope component
%
% RefHat = Mov .* exp(+i*thetahat);
%
%


%
% tdr_phdist_est.m
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

beta  = [];
rvar  = [];
nover = [];
thetahat = [];

if(nargin < 2 | nargin > 4)
  fprintf('[thetahat theta Fsig beta] = tdr_phdist_est(Ref,Mov,rthresh,order)\n');
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
theta = unwrap(angle(Ref./Mov));

% Create a polynomial design matrix %
X = fast_polytrendmtx(1,nf,1,order);
X = [ones(nf,1) [1:nf]'];

% Use only the suprathreshold frames for estimation %
thetaover = theta(indover,:);
Xover = X(indover,:);

% Perform the estimation on suprathreshold frames %
[beta rvar vdof] = fast_glmfit(thetaover,Xover);

% Compute the estimate of the phase distortion across all frames
thetahat = X*beta;

% Do a test to make sure things are ok
C = zeros(1,order+1);
C(2) = 1;
[F Fsig ces] = fast_fratio(beta,X,rvar,C);
Fsig = -log10(Fsig);

if(min(Fsig) < 100)
  fprintf('WARNING: tdr_phdist_est(): min(Fsig) = %g < 100  ',min(Fsig));
  fprintf('(threshold is %g)\n',rthresh);
end

%plot(indover,thetaover,indover,Xover*beta,'+')
%keyboard

return;
