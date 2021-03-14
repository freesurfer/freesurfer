function [kimg2, beta, phhat] = tdr_deghost(kimg,Rrow,perev,synth,kpcn)
% [kimg2 beta phhat] = tdr_deghost(kimg,Rrow,<perev>,<synth>,<kpcn>)
%
% Recons the rows of kimg with deghosting.
%
% kimg must have had readout reversals applied.
%
% perev indicates that kimg was collected with reversed
%  phase encode. Rrow must already have this taken
%  into account. perev controls whether the odd lines
%  (perev=0) or even lines (perev=1) will be used as 
%  reference. If perev is not passed, perev=0 is assumed.
%
% synth ~= 0 replaces the values in kimg with complex white noise with
% std=1 after the ghosting correction has been computed but before it
% has been applied, and before the final recon of the rows.
%   kimg = (100+randn(size(kimg))) + i*(100+randn(size(kimg)));
% Offset by 100 so that it stays gaussian after recon.
%
%
%


%
% tdr_deghost.m
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

rsubdel = 3; % region around center
rthresh = 1; % fraction of mean to use as segmentation thresh, not used

kimg2 = [];
if(nargin < 2 | nargin > 5)
  fprintf('kimg2 = tdr_deghost(kimg,Rrow,<perev>,<synth>,<kpcn>)\n');
  return;
end

if(exist('perev')==0) perev = []; end
if(isempty(perev)) perev = 0; end

if(exist('synth')==0) synth = []; end
if(isempty(synth)) synth = 0; end

if(exist('kpcn')==0) kpcn = []; end

nrows = size(kimg,1);
ncols = size(Rrow,2);

[nkrows nkcols] = size(kimg);
nkv = nkrows*nkcols;

if(~perev)
  refrows = [1:2:nrows]; % odd
  movrows = [2:2:nrows]; % even
else
  refrows = [2:2:nrows]; % even
  movrows = [1:2:nrows]; % odd
end
nrefrows = length(refrows);
nmovrows = length(movrows);

if(isempty(kpcn))
  % Recon separate images based on every-other kspace
  % line. Both will have wrap-around but neither will
  % have ghosting (may be hard to tell the difference).
  %img2ref = Rcol(:,refrows)*kimg(refrows,:)*Rrow;
  %img2mov = Rcol(:,movrows)*kimg(movrows,:)*Rrow;
  ktmp = kimg*Rrow;
  ktmpref = ktmp;
  ktmpref(movrows,:) = 0;
  img2ref = fftshift(ifft(fftshift(ktmpref,1),[],1),1);
  ktmpmov = ktmp;
  ktmpmov(refrows,:) = 0;
  img2mov = fftshift(ifft(fftshift(ktmpmov,1),[],1),1);
  
  % Choose a range of rows around the center, away from 
  % wrap-around.
  rsub = (round(nrows/2) +1)+[-rsubdel:rsubdel];
  
  % Segment the images, choose columns that exceed thresh
  img2mn    = (abs(img2ref)+abs(img2mov))/2;
  imgrsubmn = mean(img2mn(rsub,:));
  thresh    = rthresh*mean(imgrsubmn);
  csub      = find(imgrsubmn > thresh);
  csub = [round(ncols/2)-10:round(ncols/2)+10];

  % Fit the phase difference at the segmented columns.
  % Model is offset and slope.
  ph = mean(angle(img2ref(rsub,csub)./img2mov(rsub,csub)))';
  X = [ones(length(csub),1) csub'];
  beta = (inv(X'*X)*X')*ph;

  % Generate the estimate of the phase diff at all columns
  X2 = [ones(ncols,1) [1:ncols]'];
  phhat = (X2*beta)';
else
  [gcomp phhat thetahat Fsig beta Rrow] = tdr_ghostcomp(kpcn,Rrow);
end


% Split phase diff between ref and mov
vref = exp(-i*phhat/2);
vmov = exp(+i*phhat/2);

% Replace with white noise if desired. In image space with a simple
% fourier recon, each image voxel will have an mean offest of 10
% and std=1.
if(synth)
  rk0 = round(nkrows/2)+1;
  ck0 = round(nkcols/2)+1;
  kimg = sqrt(nkv)*(randn(size(kimg)) + i*randn(size(kimg)));
  kimg(rk0,ck0,:) = 10*nkv;
end

% Recon rows of all lines
kimg2 = kimg*Rrow;

% Apply phase shift to appropriate lines
kimg2(refrows,:) = repmat(vref,[nrefrows 1]) .* kimg2(refrows,:);
kimg2(movrows,:) = repmat(vmov,[nmovrows 1]) .* kimg2(movrows,:);

%nn=1:nkcols;
%plot(nn(csub),ph,nn(csub),phhat(csub));
%keyboard


%img = Rcol*kimg2;
%figure; imagesc(abs(img)); colormap(gray);
%keyboard

return;

