% tdr_rcolmat.m - computes the matrix that will reconstruct the
% columns based on the FID map and according to the time-domain
% reconstruction method.
%
% $Id: tdr_rcolmat.m,v 1.4 2003/11/07 00:42:42 greve Exp $

if(0) 
  % Input and output files
  fidmatfile  = sprintf('%s/D%2d.2.mat',topdir,TE0);
  rcolmatfile = sprintf('%s/R%2d.2.mat',topdir,TE0);
  
  % Normalize FID with abs value at first FID echo
  dnorm = 1;
  
  % Segmentation theshold the first echo of the FID Map
  rthreshfid = 0.2;
  
  boldweight = 0;
  
  % SVD regularization percentage
  regmethod = 'svdpct'; % tikhonov or svdpct 
  svdregpct = 90;
  tikregfact = 0.5; % Tikhonov regularization
end

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
    
    if(boldweight)
      pedcol = pedmat(:,imgcol);
      Dcol = Dcol .* repmat(pedcol,[1 nrows]);
    end
    
    % Compute the encoding matrix
    FcolTDR = Dcol .* Fcol;

    % Regularize and compute reconstruction matrix %
    switch(regmethod)
     case 'tikhonov' 
      tmp = FcolTDR'*FcolTDR;
      RcolTDR = inv(tmp + tikregfact*eye(size(tmp))*mean(diag(tmp)))*FcolTDR';
      FcolTDRdim = 0; % just needs to be set to something
     case 'svdpct'
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

fprintf('Saving to %s\n',rcolmatfile);
save(rcolmatfile,'fidmatfile','dnorm','regmethod',...
     'fidvol1mn','fidvol1','TE','TEList','perev',...
     'svdregpct','tikregfact','Rtdr','RtdrCond','RtdrDim',...
     'T2s','B0','epiref_dist','epiref_undist','boldweight');
fprintf('Done (%g)\n',toc);

