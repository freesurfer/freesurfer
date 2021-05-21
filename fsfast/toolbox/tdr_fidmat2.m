% tdr_fidmat2.m - program to compute the FID matrix suitable 
% for use with the time-domain reconstruction method. The
% reconstruction matrix is not computed here. Basically the
% same as tdr_fidmat.m with the following exceptions:
%    1. Only one TE is handled at a time
%    2. Only one FID map used as input
%    3. Only the central half of the columns are kept
%    4. Saves data in anatomical order
%
% Note: the images are still not flipped up-down
%
%


%
% tdr_fidmat2.m
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

fidmatversion = 2;

if(0)
  fiddir = '/home/greve/sg1/dng072203/fid1/mgh';
  TE    = 20;     % Echo time in ms
  perev = 0;
  rorev = 0;
  outmat = '/home/greve/sg1/dng072203/D20.1.B.mat';
  fidfwhm = 0;

  % Free induction decay (FID) timing parameters
  fidecho1ped = 1810;   % PED of first echo (us)
  fidechospacing = 820; % Time between FID echoes (Actually use double)
  nfidechoes = 99; % But we'll only use half (odd echoes), should get from file
  
  % EPI timing parameters - note EPI data not needed here
  epiechospacing = 470; % us
  delsamp = 30;         % us
  tDwell = 3.2;         % us
  
  % Dimensions: applies to both EPI and FID - should just get this from data
  nrows = 64;
  ncols = 128;
  nslices = 35;

end

nT2sFit = 5; % number of echoes to use to fit the T2s and B0

nv = nrows*ncols;
evenrows = [2:2:nrows];
oddrows  = [1:2:nrows];

% Start Timer
tic;

% Only keep the interior half of the columns
c1 = ncols/4 + 1;
c2 = c1 + ncols/2 -1;
colkeep = [c1:c2];
ncolskeep = length(colkeep);
nkcols = ncols;

% Only use odd FID echoes
fidechospacing_odd = 2*fidechospacing;
nfidechoes_odd = length([1:2:nfidechoes]);

% Times at which the odd fid echoes were aquired %
tfid = (fidecho1ped + fidechospacing_odd*[0:nfidechoes_odd-1]');
tfid = tfid/1000; % convert to ms

% Compute the Ideal col and row DFT reconstruction matrices
Fcol = fast_dftmtx(nrows);
Rcol = inv(Fcol);
Frow = fast_dftmtx(ncols);
Rrow = transpose(inv(Frow));

% Get the PED of each sample in the EPI
pedmat = tdr_pedmatrix(TE*1000,epiechospacing,delsamp,tDwell,...
		       nrows,ncols,perev);

% Apply the same transforms as will be applied to the EPI kspace data
% to make it match the first echo of the FID image
pedmat(evenrows,:) = fliplr(pedmat(evenrows,:));
if(perev) pedmat = flipud(pedmat);  end
if(rorev) pedmat = fliplr(pedmat);  end

% Compute the FID echo indices corresponding to the times at which
% the EPI k-space samples were taken. Weights are for lin interp.
[indA, indB, wA, wB] = fidinterp(pedmat,fidecho1ped,...
				 fidechospacing_odd,nfidechoes_odd);

D   = zeros(nrows,nrows,ncolskeep,nslices);
D0  = zeros(nrows,ncolskeep,nslices);
T2s = zeros(nrows,ncolskeep,nslices);
B0  = zeros(nrows,ncolskeep,nslices);
epiref_undist = zeros(nrows,ncolskeep,nslices);
kepiref_dist  = zeros(nrows,ncols,nslices);

if(interleaved)
  sliceorder = tdr_sliceorder(nslices,1);
else
  sliceorder = 1:nslices;
end
for acqsliceno = 1:nslices
  sliceno = sliceorder(acqsliceno); 
  %sliceno = acqsliceno;
  fprintf('acq = %d, sliceno = %d (%g) ----------------\n',...
	  acqsliceno,sliceno,toc);

  fid  = zeros(nrows,ncols,nfidechoes_odd);
  kfid = zeros(nrows,ncols,nfidechoes_odd);

  %----------------------------------------------------------------%
  % Load in all odd echos of the FID for the given slice %
  nthFIDEcho = 1;
  for FIDEcho = 1:2:nfidechoes
    %fprintf(' fid echo = %d, toc = %g\n',FIDEcho,toc);
    
    % Load the real part of this FID echo 
    mghname = sprintf('%s/echo%03dr.mgh',fiddir,FIDEcho);
    kfidr = load_mgh(mghname,acqsliceno); 
    if(isempty(kfidr)) return; end
    kfidr = transpose(kfidr); % transpose for row major
    
    % Load the imagniary part of this FID echo 
    mghname = sprintf('%s/echo%03di.mgh',fiddir,FIDEcho);
    kfidi = load_mgh(mghname,acqsliceno);
    if(isempty(kfidi)) return; end
    kfidi = transpose(kfidi); % transpose for row major
    
    % Compute the complex
    kfid_echo = kfidr + i*kfidi;
    
    % Recon the FID image
    %fidimg0 = Rcol * kfid_echo * Rrow;
    fidimg = fftshift(ifft2(fftshift(kfid_echo)));
    
    if(fidfwhm > 0)
      fidimg = tdr_smooth2d(fidimg,fidfwhm,fidfwhm);
    end

    fid(:,:,nthFIDEcho)  = fidimg;
    kfid(:,:,nthFIDEcho) = kfid_echo;
    
    nthFIDEcho = nthFIDEcho + 1;
  end % finished loading all FID echoes
  fprintf('  Finished loading fid (%g)\n',toc);
  
  % Save the first echo separately
  D0(:,:,sliceno) = abs(fid(:,colkeep,1));

  % Compute the T2* map
  nnfit = 1:nT2sFit;
  T2s(:,:,sliceno) = tdr_fidt2star(abs(fid(:,colkeep,nnfit)),tfid(nnfit),nT2sFit);
  
  % Compute the B0 map in radians/sec
  % Divide by 2*pi*123 to get parts-per-million at 3T
  B0(:,:,sliceno) = tdr_fidb0(fid(:,colkeep,1:nT2sFit),tfid(nnfit)/1000,nT2sFit);
  
  %----------------------------------------------------------------%
  % Compute the Decay Map for the given TE for this slice %

  % Simulate distorted and undistorted epis
  fprintf('  Computing simulation (%g)\n',toc);
  kepiref_dist_tmp  = zeros(nrows,ncols);
  epiref_undist_tmp  = zeros(nrows,ncols);
  for r = 1:nrows
    for c = 1:ncols
      nA = indA(r,c);
      nB = indB(r,c);
      kepiref_dist_tmp(r,c) = wA(r,c) * kfid(r,c,nA) + wB(r,c) * kfid(r,c,nB);
      epiref_undist_tmp(r,c) = ...
	  abs(wA(r,c) *  fid(r,c,nA) + wB(r,c) *  fid(r,c,nB));
    end
  end
  kepiref_dist(:,:,sliceno) = kepiref_dist_tmp;
  tmp = abs(Rcol * kepiref_dist_tmp * Rrow);
  epiref_dist(:,:,sliceno) = tmp(:,colkeep);
  epiref_undist(:,:,sliceno) = epiref_undist_tmp(:,colkeep);
    
  % Compute the decay map for each column %
  fprintf('  Computing decay map (%g)\n',toc);
  nthcol = 1;
  for imgcol = colkeep
    Dcol  = zeros(nrows,nrows);
    for imgrow = 1:nrows
      fidA = squeeze(fid(imgrow,imgcol,indA(:,imgcol)));
      fidB = squeeze(fid(imgrow,imgcol,indB(:,imgcol)));
      d = wA(:,imgcol) .* fidA + wB(:,imgcol) .* fidB;
      Dcol(:,imgrow) = d;
    end
    D(:,:,nthcol,sliceno) = Dcol;
    nthcol = nthcol + 1;
  end % img col

end %  slice

pedmat = pedmat(:,colkeep);

fprintf('\n');
fprintf('Saving mat file (%g)\n',toc);
save(outmat,'fiddir','fidecho1ped','fidechospacing','nfidechoes',...
     'fidfwhm','epiechospacing','delsamp','tDwell','TE',...
     'perev','rorev','D','D0','T2s','B0','epiref_dist','epiref_undist',...
     'sliceorder','interleaved','pedmat','colkeep','nT2sFit',...
     'kepiref_dist','fidmatversion','nkcols');

fprintf('tdr_fidmat2 done (%g)\n',toc);

return;
