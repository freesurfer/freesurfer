% tdr_rcolmat.m - computes the matrix that will reconstruct the
% columns based on the FID map and according to the time-domain
% reconstruction method.
%
% $Id: tdr_rcolmat.m,v 1.1 2003/09/22 06:00:09 greve Exp $

rcolmatfile = '/space/greve/2/users/greve/dng072203/R50.nomask.mat';
fidmatfile = '/space/greve/2/users/greve/dng072203/D50.mat';

% Normalize FID with abs value at first FID echo
dnorm = 1;

% Segmentation theshold the first echo of the FID Map
rthreshfid = 0.0;

% SVD regularization percentage
regmethod = 2; % tikhonov = 1, svdpct = 2
svdregpct = 90;
tikregfact = 0.5; % Tikhonov regularization

%-----------------------------------------%
tic; % Start Timer

fprintf('Loading fidmat file ... ');
load(fidmatfile);
fprintf('Done (%g)\n',toc);
[nrows ncols nslices] = size(D0);
nv = prod([nrows ncols nslices]);

% Compute the Ideal col and row DFT reconstruction matrices
Fcol = fast_dftmtx(nrows);
%Rcol = inv(Fcol);
%Frow = fast_dftmtx(ncols);
%Rrow = transpose(inv(Frow));

fidvol1 = abs(D0);
if(rthreshfid > 0)
  % set decay map to 1 for subthresh voxels
  fprintf('Segmenting Decay Map ... ');
  fidvol1mn = mean(reshape1d(fidvol1));
  indsubthresh = find(fidvol1 < rthreshfid*fidvol1mn);
  nsubthresh = length(indsubthresh);
  DTE = reshape(DTE,[nrows nv]);
  DTE(:,indsubthresh) = 1;
  DTE = reshape(DTE,[nrows nrows ncols nslices]);
  fprintf('Done (%g)\n',toc);
  fprintf('Found %d (%g %%) voxels below thresh\n',...
	  nsubthresh,100*nsubthresh/nv);
else
  fidvol1mn = 0;
end

Rtdr = zeros(nrows,nrows,ncols,nslices);
RtdrCond = zeros(ncols,nslices);
RtdrDim  = zeros(ncols,nslices);
for sliceno = 1:nslices
  fprintf('sliceno = %d (%g) ----------------\n',sliceno,toc);

  for imgcol = 1:ncols

    Dcol = DTE(:,:,imgcol,sliceno);

    if(dnorm)
      % Normalize the decay waveform wrt the value at 
      % the first FID echo.
      Dcol0 = transpose(D0(:,imgcol,sliceno));
      Dcol = Dcol ./ repmat(abs(Dcol0),[nrows 1]);
      %Dcol = Dcol ./ repmat(abs(Dcol(1,:)),[nrows 1]);
    end

    % Compute the encoding matrix
    FcolTDR = Dcol .* Fcol;

    % Regularize and compute reconstruction matrix %
    switch(regmethod)
     case 1 % Tikonov
      tmp = FcolTDR'*FcolTDR;
      RcolTDR = inv( tmp + tikregfact*eye(size(tmp))*mean(diag(tmp)) )*FcolTDR';
      FcolTDRdim = 0; % just needs to be set to something
     case 2 % svd
      [FcolTDRreg FcolTDRdim] = fast_svdregpct(FcolTDR,svdregpct);
      RcolTDR = inv(FcolTDRreg);
    end
    
    % Normalize the rows of the recon matrix %
    RcolTDRmag = sqrt(abs(sum(RcolTDR.*conj(RcolTDR),2)));
    RcolTDR = RcolTDR./repmat(RcolTDRmag,[1 nrows]);

    Rtdr(:,:,imgcol,sliceno) = RcolTDR;
    RtdrCond(imgcol,sliceno) = cond(RcolTDR);
    RtdrDim(imgcol,sliceno) = FcolTDRdim;
  
  end

end % slice

fprintf('Saving ...\n');
save(rcolmatfile,'fidmatfile','dnorm','regmethod',...
     'fidvol1mn','fidvol1','TE','perev',...
     'svdregpct','tikregfact','Rtdr','RtdrCond','RtdrDim');
fprintf('Done (%g)\n',toc);
