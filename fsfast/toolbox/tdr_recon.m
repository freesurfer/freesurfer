% tdr_recon.m - routine to actually perform the k-space
% reconstruction according to the time-domain reconstruction
% algorithm 
%
% $Id: tdr_recon.m,v 1.5 2003/11/06 19:43:41 greve Exp $
tic;


if(0)
  %topdir = '/space/greve/2/users/greve/dng072203';
  topdir = '/space/greve/2/users/greve/fb-105.2/';
  TE0 = 50;
  run = 1;
  usefid = 1;
  nframes = 85;
  %nframes = 1;
  
  if(TE0 == 30)
    kepidir = sprintf('%s/rawk/sing-echo-r%d/mgh',topdir,run);
    epiecho = 1;
  else
    kepidir = sprintf('%s/sm%d/mgh',topdir,run);
    if(TE0 == 20) epiecho = 1;
    else          epiecho = 2;
    end
  end
  
  %rcolmatfile = sprintf('%s/R%2d.1.mat',topdir,TE0);
  rcolmatfile = sprintf('%s/tdr/R%2d.1.mat',topdir,TE0);

  sessdir = sprintf('%s/tdr-te%2d-r%d',topdir,TE0,run);
  %bhdrfile = sprintf('%s/siemens-te30/bold/001/f.bhdr',topdir);
  bhdrfile = [];
  
  funcdir = sprintf('%s/bold/001',sessdir);
  funcstem = sprintf('%s/f',funcdir);
  mkdirp(funcdir);

  fidoutdir = sprintf('%s/fid/001',sessdir);
  mkdirp(fidoutdir);
  fidstem = sprintf('%s/f',fidoutdir);
end

if(~isempty(bhdrfile))
  mristruct = fast_ldbhdr(bhdrfile);
  if(isempty(mristruct))
    fprintf('ERROR: could not load bhdrfile %s\n',bhdrfile);
    return;
  end
else
  mristruct = [];
end

% EPI gradient and ADC timing
tDwell     = 3.2; % usec
tRampUp    = 140; % usec
tFlat      = 190; % usec
tRampDown  = 140; % usec
tDelSamp   = 30;  % usec
tdelay     = 1.0;

%----------------------------------------------------%

fprintf('Loading rcolmat file ... ');
load(rcolmatfile);
fprintf('Done (%g)\n',toc);
[nrows ncols nslices] = size(fidvol1);
nv = prod([nrows ncols nslices]);
sliceorder = [1:2:nslices 2:2:nslices];

[kvec0 gvec0] = kspacevector2(ncols,tDwell,tRampUp,tFlat,...
			      tRampDown,tDelSamp,tdelay);
kvec = kvec0;
if(perev) kvec = fliplr(kvec); end

% Compute the Ideal col and row DFT reconstruction matrices
Frow = fast_dftmtx(kvec);
Frow = fast_svdregpct(Frow,90);
Rrow = transpose(inv(Frow));
Fcol = fast_dftmtx(nrows);
Rcol = inv(Fcol);


nn = 1:ncols;
nnrev = ncols:-1:1;
oddrows  = 1:2:nrows;
evenrows = 2:2:nrows;

kepi_rfile = sprintf('%s/echo%03dr.mgh',kepidir,epiecho);
kepi_ifile = sprintf('%s/echo%03di.mgh',kepidir,epiecho);

% Get Phase Encode Shift for all frames %
fprintf('Estimating shift in phase encode at each frame (%g)\n',toc);
nthecho = find(TEList == TE);
vref = epiref_dist(:,:,:,nthecho);
pevoxshift = zeros(nframes,1);
for frame = 1:nframes
  kr = load_mgh(kepi_rfile,[],frame);
  kr = permute(kr,[2 1 3]);
  ki = load_mgh(kepi_ifile,[],frame);
  ki = permute(ki,[2 1 3]);
  kvol = kr + i*ki;
  [krvol, dgbeta] = tdr_recon_rows(kvol,Rrow,perev);
  ind = 40:87;
  pevoxshift(frame) = tdr_peshift(vref(:,ind,:),krvol(:,ind,:));
end
fprintf('   done (%g).\n',toc);
fname = sprintf('%s-peshift.dat',funcstem);
fp = fopen(fname,'w');
fprintf(fp,'%g\n',pevoxshift);
fclose(fp);

vol = zeros(nrows,ncols/2,nslices,nframes);
for nthAcqSlice = 1:nslices
  fprintf('nthAcqSlice = %d, toc = %g\n',nthAcqSlice,toc);
  Rtdrslice = Rtdr(:,:,:,nthAcqSlice);
  
  for frame = 1:nframes
  
    kepi_r = load_mgh(kepi_rfile,nthAcqSlice,frame);
    if(isempty(kepi_r)) return; end
    kepi_r = kepi_r';
    
    kepi_i = load_mgh(kepi_ifile,nthAcqSlice,frame);
    if(isempty(kepi_i)) return; end
    kepi_i = kepi_i';
    
    kepi = kepi_r + i * kepi_i;
    kepi = tdr_kshift(kepi,pevoxshift(frame));
    
    % Apply transforms to make EPI k-space image match the 
    % that of the first echo of the FID. These same transforms
    % must have been applied to the PED matrix (see tdr_fidmat).
    kepi(evenrows,:) = fliplr(kepi(evenrows,:));
    % Flipud already done in unpack (mri_convert_mdh)
    %if(perev) kepi = flipud(kepi);  end 
    
    % Recon the rows %
    kepi2 = tdr_deghost(kepi,Rrow,perev);
    %kepi2 = kepi*Rrow; % without deghosting
    
    if(~usefid)
      img = abs(Rcol * kepi2);
    else
      img = zeros(64,128);
      for imgcol = 1:ncols
	Rtdrcol = Rtdrslice(:,:,imgcol);
	img(:,imgcol) = abs(Rtdrcol*kepi2(:,imgcol));
      end
    end
    
    % Extract the center, and flip
    img2 = flipud(img(:,33:end-32));
    
    % Shift so that it matches the siemens recon %
    img3(2:64,:) = img2(1:63,:);
    sliceno = sliceorder(nthAcqSlice);
    vol(:,:,sliceno,frame) = img3;
  
  end % frame
  
end

% Rescale for bshort, save scaling
volmin = min(reshape1d(vol));
volmax = max(reshape1d(vol));
fscale = (2^14-1)/(volmax-volmin);
vol = fscale*(vol-volmin);
scalefile = sprintf('%s.scale.dat',funcstem);
fp = fopen(scalefile,'w');
fprintf(fp,'%g %g %g\n',volmin,volmax,fscale);
fclose(fp);

% Save the functional volume
fast_svbslice(vol,funcstem,-1,'bshort',mristruct);

if(~isempty(fidstem))
  % Prep the fidimg
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
  
  % Save the fidimg
  fast_svbslice(fidvol1,fidstem,-1,[],mristruct);
end

fprintf('tdr_recon done %g\n',toc);

return;


tit = sprintf('tdr-vs-fid: TE=%g, usefid=%d',TE0,usefid);
swapview('-init','-v1',vol(:,:,:,1),'-v2',fidvol1,'-title',tit);


