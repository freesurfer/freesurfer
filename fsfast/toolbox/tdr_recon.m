% tdr_recon.m - routine to actually perform the k-space
% reconstruction according to the time-domain reconstruction
% algorithm 
%
% $Id: tdr_recon.m,v 1.8 2004/01/22 00:52:02 greve Exp $
tic;


if(0)
  topdir = '/space/greve/1/users/greve/dng072203';
  kepidir = sprintf('%s/rawk/230855/mgh',topdir);
  rcolmatfile  = sprintf('%s/tdr2/R20.mat',topdir);  
  epiecho = 1;

  nframes = 80;
  usefid = 1;

  sessdir = sprintf('%s/dng-tdr20b',topdir);
  bhdrfile = 'dng-tdr/bold/001/f.bhdr';
  
  funcdir = sprintf('%s/bold/001',sessdir);
  funcstem = sprintf('%s/f',funcdir);
  mkdirp(funcdir);

  fidoutdir = sprintf('%s/fid/001',sessdir);
  mkdirp(fidoutdir);
  fidstem = sprintf('%s/f',fidoutdir);

  % EPI gradient and ADC timing
  tDwell     = 3.2; % usec
  tRampUp    = 140; % usec
  tFlat      = 190; % usec
  tRampDown  = 140; % usec
  tDelSamp   = 30;  % usec
  tdelay     = 1.0;
  
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


%----------------------------------------------------%

fprintf('Loading rcolmat file ... ');
load(rcolmatfile);
if(exist('rorev') ~= 1) rorev = 0; end
fprintf('Done (%g)\n',toc);
[nrows ncols nslices] = size(fidvol1);
nv = prod([nrows ncols nslices]);
sliceorder = [1:2:nslices 2:2:nslices];
nkcols = 2*ncols;

kvec = kspacevector2(nkcols,tDwell,tRampUp,tFlat,...
		     tRampDown,tDelSamp,tdelay);
if(perev & ~rorev) kvec = fliplr(kvec); end

% Determine Deghosting Reference:
%  0 = use first  (odd) lines
%  1 = use second (even) lines
DeghostRef = xor(perev,rorev);
fprintf('DeghostRef = %d (perev = %d, rorev = %d)\n',...
	DeghostRef,perev,rorev);

% Compute the Ideal col and row DFT reconstruction matrices
Frow = fast_dftmtx(kvec);
Frow = Frow(:,colkeep);
Rrow = transpose(inv(Frow'*Frow)*Frow');
%Frow = fast_svdregpct(Frow,90);
%Rrow = transpose(inv(Frow));
Fcol = fast_dftmtx(nrows);
Rcol = inv(Fcol);

nn = 1:ncols;
nnrev = ncols:-1:1;
oddrows  = 1:2:nrows;
evenrows = 2:2:nrows;

kepi_rfile = sprintf('%s/echo%03dr.mgh',kepidir,epiecho);
kepi_ifile = sprintf('%s/echo%03di.mgh',kepidir,epiecho);

if(fixpedrift)
  % Get Phase Encode Shift for all frames %
  fprintf('Estimating shift in phase encode at each frame (%g)\n',toc);
  vref = epiref_dist;
  pevoxshift = zeros(nframes,1);
  for frame = 1:nframes
    if(nframes > 100) 
      fprintf('  frame = %d (%g)  pevoxshift\n',frame,toc);
    end
    kr = load_mgh(kepi_rfile,[],frame);
    kr = permute(kr,[2 1 3]);
    ki = load_mgh(kepi_ifile,[],frame);
    ki = permute(ki,[2 1 3]);
    kvol = kr + i*ki;
    [krvol, dgbeta] = tdr_recon_rows(kvol,Rrow,DeghostRef);
    c1 = round(ncols/4);
    c2 = c1 + round(ncols/2) - 1;
    ind = c1:c2;
    pevoxshift(frame) = tdr_peshift(vref(:,ind,:),krvol(:,ind,:));
  end
  clear kr ki kvol;
  fprintf('   done (%g).\n',toc);
  fname = sprintf('%s-peshift.dat',funcstem);
  fp = fopen(fname,'w');
  fprintf(fp,'%g\n',pevoxshift);
  fclose(fp);
end

dgbeta = zeros(2,1,nslices,nframes);
vol = zeros(nrows,ncols,nslices,nframes);
%vol = zeros(nrows,ncols,nslices,nframes);
for nthAcqSlice = 1:nslices
  sliceno = sliceorder(nthAcqSlice);
  fprintf('nthAcqSlice = %d, slice = %d, toc = %g\n',...
	  nthAcqSlice,sliceno,toc);
  %Rtdrslice = Rtdr(:,:,:,nthAcqSlice);
  Rtdrslice = Rtdr(:,:,:,sliceno);
  
  %fslice = zeros(nrows,ncols,nframes);
  for frame = 1:nframes
  
    kepi_r = load_mgh(kepi_rfile,nthAcqSlice,frame);
    if(isempty(kepi_r)) return; end
    kepi_r = kepi_r';
    
    kepi_i = load_mgh(kepi_ifile,nthAcqSlice,frame);
    if(isempty(kepi_i)) return; end
    kepi_i = kepi_i';
    
    kepi = kepi_r + i * kepi_i;
    if(fixpedrift)
      kepi = tdr_kshift(kepi,pevoxshift(frame));
    end
    
    % Apply transforms to make EPI k-space image match the 
    % that of the first echo of the FID. These same transforms
    % must have been applied to the PED matrix (see tdr_fidmat).
    kepi(evenrows,:) = fliplr(kepi(evenrows,:));
    % Flipud already done in unpack (mri_convert_mdh)
    %if(perev) kepi = flipud(kepi);  end 
    
    % Recon the rows %
    [kepi2 dgbetatmp] = tdr_deghost(kepi,Rrow,DeghostRef);
    dgbeta(:,1,sliceno,frame) = dgbetatmp;
    %kepi2 = kepi*Rrow; % without deghosting
    
    if(~usefid)
      img = abs(Rcol * kepi2);
    else
      %img = zeros(64,128);
      img = zeros(nrows,ncols);
      for imgcol = 1:ncols
	Rtdrcol = Rtdrslice(:,:,imgcol);
	img(:,imgcol) = abs(Rtdrcol*kepi2(:,imgcol));
      end
    end
    
    % Extract the center, and flip
    %img2 = flipud(img(:,33:end-32));
    img2 = flipud(img);
    
    % Shift so that it matches the siemens recon %
    img3 = img2;
    img3(2:end,:) = img2(1:end-1,:);
    vol(:,:,sliceno,frame) = img3;
    %fslice(:,:,frame) = img3 * 1e5;
  
  end % frame
  
  %fprintf('Saving slice to %s (%g)\n',funcstem,toc);
  %fast_svbslice(fslice,funcstem,sliceno-1,'bfloat',mristruct);
end

if(~isempty(fidstem))
  if(0)
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
  else
    fidvol1 = flipdim(fidvol1,1);
    fidvol1(2:end,:,:) = fidvol1(1:end-1,:,:);
  end
  
  % Save the fidimg
  fast_svbslice(fidvol1,fidstem,-1,[],mristruct);
end

% Save the deghosting parameters %
fname = sprintf('%s.dgbeta',funcstem);
fast_svbslice(dgbeta,fname,-1,'bfloat',mristruct);
if(1)
fname = sprintf('%s.dgbeta.int',funcstem);
fp = fopen(fname,'w');
fmt = repmat('%g ',[1 nslices]);
fmt = sprintf('%s\\n',fmt);
fprintf(fp,fmt,squeeze(dgbeta(1,1,:,:)));
fclose(fp);
fname = sprintf('%s.dgbeta.slope',funcstem);
fp = fopen(fname,'w');
fprintf(fp,fmt,squeeze(dgbeta(2,1,:,:)));
fclose(fp);
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

fprintf('tdr_recon done %g\n',toc);
return;


tit = sprintf('tdr-vs-fid: TE=%g, usefid=%d',TE,usefid);
swapview('-init','-v1',vol(:,:,:,1),'-v2',fidvol1,'-title',tit);


