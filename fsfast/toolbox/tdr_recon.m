% tdr_recon.m - routine to actually perform the k-space
% reconstruction according to the time-domain reconstruction
% algorithm 
%
% $Id: tdr_recon.m,v 1.1 2003/09/22 06:00:09 greve Exp $

rcolmatfile = '/space/greve/2/users/greve/dng072203/R50.nomask.mat';
kepidir = '/space/greve/2/users/greve/dng072203/de1/mgh';
epiecho = 2;

usefid = 1;

% EPI gradient and ADC timing
tDwell     = 3.2; % usec
tRampUp    = 140; % usec
tFlat      = 190; % usec
tRampDown  = 140; % usec
tDelSamp   = 30;  % usec
tdelay     = 1.0;

%----------------------------------------------------%
tic;

fprintf('Loading rcolmat file ... ');
load(rcolmatfile);
fprintf('Done (%g)\n',toc);
[nrows ncols nslices] = size(fidvol1);
nv = prod([nrows ncols nslices]);
sliceorder = [1:2:nslices 2:2:nslices];

% Compute the Ideal col and row DFT reconstruction matrices
Fcol = fast_dftmtx(nrows);
Rcol = inv(Fcol);
Frow = fast_dftmtx(ncols);
Rrow = transpose(inv(Frow));

[kvec0 gvec0] = kspacevector2(ncols,tDwell,tRampUp,tFlat,...
			      tRampDown,tDelSamp,tdelay);
kvec = kvec0;

nn = 1:ncols;
nnrev = ncols:-1:1;
oddrows  = 1:2:nrows;
evenrows = 2:2:nrows;

kepi_rfile = sprintf('%s/echor%d.mgh',kepidir,epiecho-1);
kepi_ifile = sprintf('%s/echoi%d.mgh',kepidir,epiecho-1);

vol = zeros(nrows,ncols/2,nslices);
for nthAcqSlice = 1:nslices
  fprintf('nthAcqSlice = %d, toc = %g\n',nthAcqSlice,toc);
  
  kepi_r = load_mgh(kepi_rfile,nthAcqSlice,1);
  if(isempty(kepi_r)) return; end
  kepi_r = kepi_r';
  
  kepi_i = load_mgh(kepi_ifile,nthAcqSlice,1);
  if(isempty(kepi_i)) return; end
  kepi_i = kepi_i';
  
  kepi = kepi_r + i * kepi_i;
  
  % Apply transforms to make EPI k-space image match the 
  % that of the first echo of the FID. These same transforms
  % must have been applied to the PED matrix (see tdr_fidmat).
  kepi(evenrows,:) = fliplr(kepi(evenrows,:));
  % Flipud already done in unpack (mri_convert_mdh)
  %if(perev) kepi = flipud(kepi);  end 

  % Recon the rows %
  kepi2 = fast_deghost(kepi,kvec,perev);
  %kepi2 = kepi*Rrow; % without deghosting

  if(~usefid)
    img = abs(Rcol * kepi2);
  else
    img = zeros(64,128);
    for imgcol = 1:ncols
      Rcoltdr = Rtdr(:,:,imgcol,nthAcqSlice);
      img(:,imgcol) = abs(Rcoltdr*kepi2(:,imgcol));
    end
  end

  % Extract the center, and flip
  img2 = flipud(img(:,33:end-32));

  % Shift so that it matches the siemens recon %
  img3(2:64,:) = img2(1:63,:);
  frame = 1;
  sliceno = sliceorder(nthAcqSlice);
  vol(:,:,sliceno,frame) = img3;
  
end

fidvol10 = fidvol1;
fidvol1 = fidvol10;
fidvol1 = flipdim( fidvol1(:,33:end-32,:,:), 1);
fidvol1(2:64,:,:) = fidvol1(1:63,:,:);

h = ceil(nslices/2);
a = [1:h]';
b = [h+1:2*h]';
c = reshape1d([a b]');
s = c(1:nslices);
fidvol1 = fidvol1(:,:,s);

swapview('-init','-v1',fidvol1,'-v2',vol(:,:,:,1),'-title','fid-vs-tdr');


