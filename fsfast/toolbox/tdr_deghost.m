function [kimg2, beta] = tdr_deghost(kimg,Rrow,perev)
% kimg2 = tdr_deghost(kimg,Rrow,<perev>)
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
% $Id: tdr_deghost.m,v 1.4 2003/10/30 19:29:26 greve Exp $
%

rsubdel = 3; % region around center
rthresh = 1; % fraction of mean to use as segmentation thresh

kimg2 = [];
if(nargin ~= 2 & nargin ~= 3)
  fprintf('kimg2 = tdr_deghost(kimg,Rrow,<perev>)\n');
  return;
end

if(exist('perev')==0) perev = []; end
if(isempty(perev)) perev = 0; end

[nrows ncols] = size(kimg);

if(~perev)
  refrows = [1:2:nrows]; % odd
  movrows = [2:2:nrows]; % even
else
  refrows = [2:2:nrows]; % even
  movrows = [1:2:nrows]; % odd
end

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
rsub = (nrows/2 +1)+[-rsubdel:rsubdel];

% Segment the images, choose columns that exceed thresh
img2mn = (abs(img2ref)+abs(img2mov))/2;
imgrsubmn = mean(img2mn(rsub,:));
thresh = rthresh*mean(imgrsubmn);
csub = find(imgrsubmn > thresh);

% Fit the phase difference at the segmented columns.
% Model is offset and slope.
ph = mean(angle(img2ref(rsub,csub)./img2mov(rsub,csub)))';
X = [ones(length(csub),1) csub'];
beta = (inv(X'*X)*X')*ph;

% Generate the estimate of the phase diff at all columns
X2 = [ones(ncols,1) [1:ncols]'];
phsynth = (X2*beta)';

% Split phase diff between ref and mov
vref = exp(-i*phsynth/2);
vmov = exp(+i*phsynth/2);

% Recon rows of all lines
kimg2 = kimg*Rrow;

% Apply phase shift to appropriate lines
kimg2(refrows,:) = repmat(vref,[nrows/2 1]) .* kimg2(refrows,:);
kimg2(movrows,:) = repmat(vmov,[nrows/2 1]) .* kimg2(movrows,:);

%img = Rcol*kimg2;
%figure; imagesc(abs(img)); colormap(gray);
%keyboard

return;

