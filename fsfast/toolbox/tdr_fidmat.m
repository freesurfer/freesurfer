% tdr_fidmat.m - program to compute the FID matrix suitable 
% for use with the time-domain reconstruction method. The
% reconstruction matrix is not computed here.
%
% $Id: tdr_fidmat.m,v 1.2 2003/09/25 02:12:39 greve Exp $

addpath /homes/4/greve/sg2/dale-matlab/utils; % for smoother

if(0)

fiddirlist = '';
fiddirlist = strvcat(fiddirlist,'/space/greve/2/users/greve/fb-104.2/rawk/fid2/mgh');

fidfwhm = 0; % smoothing for fid

% Echo times to generate the FID for 
TEList      = [20 30 50];     % Echo time in ms
TEPERevList = [ 0  0  1];  % Flag indicating reversal in PE direction
outlist = '';
outlist = strvcat(outlist,'/space/greve/2/users/greve/fb-104.2/D20.2.mat');
outlist = strvcat(outlist,'/space/greve/2/users/greve/fb-104.2/D30.2.mat');
outlist = strvcat(outlist,'/space/greve/2/users/greve/fb-104.2/D50.2.mat');

end

% Free induction decay (FID) timing parameters
fidecho1ped = 1810;   % PED of first echo (us)
fidechospacing = 820; % Time between FID echoes (Actually use double)
nfidechoes = 99; % But we'll only use half (odd echoes)

% EPI timing parameters - note EPI data not needed here
epiechospacing = 470; % us
delsamp = 30;         % us
tDwell = 3.2;         % us

% Dimensions: applies to both EPI and FID - should just get this from data
nrows = 64;
ncols = 128;
nslices = 35;
nv = nrows*ncols;
evenrows = [2:2:nrows];
oddrows  = [1:2:nrows];

% Start Timer
tic;

% Only use odd FID echoes
fidechospacing_odd = 2*fidechospacing;
nfidechoes_odd = length([1:2:nfidechoes]);

% Compute the Ideal col and row DFT reconstruction matrices
Fcol = fast_dftmtx(nrows);
Rcol = inv(Fcol);
Frow = fast_dftmtx(ncols);
Rrow = transpose(inv(Frow));

nTEs = length(TEList);
nFIDMaps = size(fiddirlist,1);

D = zeros(nrows,nrows,ncols,nslices,nTEs);
D0 = zeros(nrows,ncols,nslices);

for sliceno = 1:nslices
  fprintf('sliceno = %d (%g) ----------------\n',sliceno,toc);
  fid = zeros(nrows,ncols,nfidechoes_odd);

  %----------------------------------------------------------------%
  % Load in all odd echos of the FID for the given slice %
  nthFIDEcho = 1;
  for FIDEcho = 1:2:nfidechoes
    %fprintf('echo = %d, toc = %g\n',e,toc);
    
    kfidr = 0;
    kfidi = 0;
    for nthFIDMap = 1:nFIDMaps
      fiddir = deblank(fiddirlist(nthFIDMap,:));

      % Load the real part of this FID echo 
      mghname = sprintf('%s/echo%03dr.mgh',fiddir,FIDEcho);
      kfidrtmp = load_mgh(mghname,sliceno); 
      if(isempty(kfidrtmp)) return; end
      kfidr = kfidr + kfidrtmp';% transpose for row major
    
      % Load the imagniary part of this FID echo 
      mghname = sprintf('%s/echo%03di.mgh',fiddir,FIDEcho);
      kfiditmp = load_mgh(mghname,sliceno);
      if(isempty(kfiditmp)) return; end
      kfidi = kfidi + kfiditmp'; % transpose for row major
    end
    
    % Compute the complex
    kfid = (kfidr + i*kfidi)/nFIDMaps;
    
    % Recon the FID image
    fidimg = Rcol * kfid * Rrow;

    if(fidfwhm > 0)
      fidimg = smooth2d(fidimg,fidfwhm,fidfwhm);
    end

    fid(:,:,nthFIDEcho) = fidimg;
    
    nthFIDEcho = nthFIDEcho + 1;
  end % finished loading all FID echoes

  % Save the first echo separately
  D0(:,:,sliceno) = fid(:,:,1);

  %----------------------------------------------------------------%
  % Compute the Decay Map for each echo for this slice %
  nthTE = 1;
  for TE = TEList
    fprintf('  TE = %g (%g) \n',TE,toc);

    % Get the PED of each sample in the EPI
    perev = TEPERevList(nthTE);
    pedmat = pedmatrix2(TE*1000,epiechospacing,delsamp,tDwell,nrows,ncols,perev);

    % Apply the same transforms as will be applied to the EPI kspace data
    % to make it match the first echo of the FID 
    pedmat(evenrows,:) = fliplr(pedmat(evenrows,:));
    if(perev)   pedmat = flipud(pedmat);  end

    % Compute the FID echo indices corresponding to the times at which
    % the EPI k-space samples were taken. Weights are for lin interp.
    [indA, indB, wA, wB] = fidinterp(pedmat,fidecho1ped,...
				     fidechospacing_odd,nfidechoes_odd);
  
    % Compute the decay map for each column %
    for imgcol = 1:ncols
      Dcol  = zeros(nrows,nrows);
      for imgrow = 1:nrows
	fidA = squeeze(fid(imgrow,imgcol,indA(:,imgcol)));
	fidB = squeeze(fid(imgrow,imgcol,indB(:,imgcol)));
	d = wA(:,imgcol) .* fidA + wB(:,imgcol) .* fidB;
	Dcol(:,imgrow) = d;
      end
      %fprintf('sliceno = %d, %g\n',sliceno,mar(Dcol));
      %if(sliceno==10 & imgcol == 65) keyboard; end
      D(:,:,imgcol,sliceno,nthTE) = Dcol;
    end % img col

    nthTE = nthTE + 1;
  end % TE

end %  slice

% Save the results %
for n = 1:nTEs

  fprintf('Saving TE=%g (%g)\n',TEList(n),toc);
  outfile = deblank(outlist(n,:));
  TE = TEList(n);
  perev = TEPERevList(n);
  DTE = D(:,:,:,:,n);
  save(outfile,'fiddirlist','fidecho1ped','fidechospacing','nfidechoes',...
     'fidfwhm','epiechospacing','delsamp','tDwell','TE','perev', ...
     'DTE','D0');
end

toc