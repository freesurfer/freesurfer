% tdr_fidmat.m - program to compute the FID matrix suitable 
% for use with the time-domain reconstruction method. The
% reconstruction matrix is not computed here.
%
%


%
% tdr_fidmat.m
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

%addpath /homes/4/greve/sg2/dale-matlab/utils; % for smoother

if(0)

fiddirlist = '';
fiddirlist = strvcat(fiddirlist,...
     '/space/greve/2/users/greve/fb-104.2/rawk/fid2/mgh');

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
nfidechoes = 81; % But we'll only use half (odd echoes)

% EPI timing parameters - note EPI data not needed here
%epiechospacing = 470; % us
%delsamp = 30;         % us
%tDwell = 3.2;         % us
% For Martin
epiechospacing = 530; % us
delsamp = 60;         % us
tDwell = 1.6;         % us
% Echo time is 40 ms - set in tdr-fidmat script

% Dimensions: applies to both EPI and FID - should just get this from data
nrows = 128;  % martin
ncols = 256;  % martin
nslices = 23; % martin 
nv = nrows*ncols;
evenrows = [2:2:nrows];
oddrows  = [1:2:nrows];

% Start Timer
tic;

% Only use odd FID echoes
fidechospacing_odd = 2*fidechospacing;
nfidechoes_odd = length([1:2:nfidechoes]);

% Times at which the odd fid echoes were aquired %
tfid = (fidecho1ped + fidechospacing_odd*[0:nfidechoes_odd-1]');
tfid = tfid/1000; % convert to ms
nT2sFit = 5; % number of echoes to use to fit the T2s

% Compute the Ideal col and row DFT reconstruction matrices
Fcol = fast_dftmtx(nrows);
Rcol = inv(Fcol);
Frow = fast_dftmtx(ncols);
Rrow = transpose(inv(Frow));

nTEs = length(TEList);
nFIDMaps = size(fiddirlist,1);

D   = zeros(nrows,nrows,ncols,nslices,nTEs);
D0  = zeros(nrows,ncols,nslices);
T2s = zeros(nrows,ncols,nslices);
B0  = zeros(nrows,ncols,nslices);
epiref_undist = zeros(nrows,ncols,nslices,nTEs);

sliceorder = [1:2:nslices 2:2:nslices];
for acqsliceno = 1:nslices
  %sliceno = sliceorder(acqsliceno); 
  sliceno = acqsliceno;
  fprintf('acq = %d, sliceno = %d (%g) ----------------\n',...
	  acqsliceno,sliceno,toc);

  fid  = zeros(nrows,ncols,nfidechoes_odd);
  kfid = zeros(nrows,ncols,nfidechoes_odd);

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
      kfidrtmp = load_mgh(mghname,acqsliceno); 
      if(isempty(kfidrtmp)) return; end
      kfidr = kfidr + kfidrtmp';% transpose for row major
    
      % Load the imagniary part of this FID echo 
      mghname = sprintf('%s/echo%03di.mgh',fiddir,FIDEcho);
      kfiditmp = load_mgh(mghname,acqsliceno);
      if(isempty(kfiditmp)) return; end
      kfidi = kfidi + kfiditmp'; % transpose for row major
    end
    
    % Compute the complex
    kfid_echo = (kfidr + i*kfidi)/nFIDMaps;
    
    % Recon the FID image
    fidimg = Rcol * kfid_echo * Rrow;

    if(fidfwhm > 0)
      fidimg = tdr_smooth2d(fidimg,fidfwhm,fidfwhm);
    end

    fid(:,:,nthFIDEcho)  = fidimg;
    kfid(:,:,nthFIDEcho) = kfid_echo;
    
    nthFIDEcho = nthFIDEcho + 1;
  end % finished loading all FID echoes

  % Save the first echo separately
  D0(:,:,sliceno) = fid(:,:,1);

  % Compute the T2* map
  T2s(:,:,sliceno) = tdr_fidt2star(abs(fid),tfid,nT2sFit);
  
  % Compute the B0 map in radians/sec
  % Divide by 2*pi*123 to get parts-per-million at 3T
  B0(:,:,sliceno) = tdr_fidb0(fid,tfid/1000,nT2sFit);
  
  %----------------------------------------------------------------%
  % Compute the Decay Map for each echo for this slice %
  nthTE = 1;
  for TE = TEList
    fprintf('  TE = %g (%g) \n',TE,toc);

    % Get the PED of each sample in the EPI
    perev = TEPERevList(nthTE);
    pedmat = tdr_pedmatrix(TE*1000,epiechospacing,delsamp,tDwell,...
			   nrows,ncols,perev);

    % Apply the same transforms as will be applied to the EPI kspace data
    % to make it match the first echo of the FID 
    pedmat(evenrows,:) = fliplr(pedmat(evenrows,:));
    if(perev)   pedmat = flipud(pedmat);  end

    % Compute the FID echo indices corresponding to the times at which
    % the EPI k-space samples were taken. Weights are for lin interp.
    [indA, indB, wA, wB] = fidinterp(pedmat,fidecho1ped,...
				     fidechospacing_odd,nfidechoes_odd);
  
    % Simulate distorted and undistorted epis
    kepiref_dist  = zeros(nrows,ncols);
    for r = 1:nrows
      for c = 1:ncols
	nA = indA(r,c);
	nB = indB(r,c);
	kepiref_dist(r,c) = wA(r,c) * kfid(r,c,nA) + wB(r,c) * kfid(r,c,nB);
	epiref_undist(r,c,sliceno,nthTE) = ...
	    abs(wA(r,c) *  fid(r,c,nA) + wB(r,c) *  fid(r,c,nB));
      end
    end
    epiref_dist(:,:,sliceno,nthTE) = abs(Rcol * kepiref_dist * Rrow);
    
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
       'fidfwhm','epiechospacing','delsamp','tDwell','TE','TEList', ...
       'perev','DTE','D0','T2s','B0','epiref_dist','epiref_undist',...
       'sliceorder','pedmat');
  clear DTE;
end

toc
