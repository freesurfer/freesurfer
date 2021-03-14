function err = myMRIwrite(mri,fstring,datatype,tempdir)
%
% Version modified by Eugenio to avoid collisions in temp directory (and
% also to allow using a temp directory other than /tmp/')
%
% err = MRIwrite(mri,fstring,datatype)
%
% Writes an MRI volume based on the fstring. fstring can be:
%  1. MGH file. Eg, f.mgh or f.mgz
%  2. BHDR file Eg, f.bhdr. Result will be written to a bfloat
%     volume, eg f_000.bfloat.
%  3. NIFIT file, Eg, f.nii, f.nii.gz (uncompressed and compressed)
%
% mri should be a structure like that read by MRIread.m The goemetry
% (ie, direction cosines, voxel resolution, and P0 are all recomputed
% from mri.vox2ras0. So, if in the course of analysis, you changed
% mri.x_r, this change will not be reflected in the output volume.
%
% The only thing you need to fill-in in the mri struct is the mri.vol.
% All other fields will be filled in with defaulted values. 
% Fields are: vol, tr, vox2ras0, te, ti, flip_angle.
%
% 
% When writing in bhdr format, the default will be bfloat. If you want
% bshort, then set mri.outbext = 'bshort'. When a bhdr file is read in
% with MRIread(), mri.srcbext is set to either bshort or bfloat, so to
% keep the same precision set mri.outbext = mri.srcbext.  This only
% applies to bhdr format.
% 
% datatype can be uchar, short, int, float, double, ushort,
% uint. Only applies to nifti.

%
% MRIwrite.m
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


err = 1;

if(nargin < 2 || nargin > 4)
  fprintf('err = MRIwrite(mri,fstring,<datatype>,<tempdir>)\n');
  return;
end
if(exist('tempdir','var')~=1) || isempty(tempdir), tempdir = '/tmp/'; end
if(exist('datatype','var')~=1) || isempty(datatype), datatype = 'float'; end

if(~isfield(mri,'vol'))
  fprintf('ERROR: MRIwrite: structure does not have a vol field\n');
  return;
end

if tempdir(end)~='/', tempdir=[tempdir '/']; end


vsz = size(mri.vol);
nvsz = length(vsz);
if(nvsz ~= 4) vsz = [vsz ones(1,4-nvsz)]; end
if(~isfield(mri,'volsize'))
  mri.volsize = [vsz(1) vsz(2) vsz(3)];
end
if(~isfield(mri,'nframes'))  mri.nframes = vsz(4); end
if(~isfield(mri,'volres'))  mri.volres = [1 1 1];end
if(~isfield(mri,'tr')) mri.tr = 0; end
if(~isfield(mri,'te')) mri.te = 0; end
if(~isfield(mri,'ti')) mri.ti = 0; end
if(~isfield(mri,'flip_angle')) mri.flip_angle = 0;end
if(~isfield(mri,'vox2ras0'))  mri.vox2ras0 = eye(4);end

  
[fspec fstem fmt] = MRIfspec(fstring,0); % 0 = turn off checkdisk
if(isempty(fspec))
  fprintf('ERROR: could not determine format of %s\n',fstring);
  return;
end

switch(fmt)
 case {'mgh','mgz'} %----------- MGH/MGZ ------------%
  M = mri.vox2ras0;
  mr_parms = [mri.tr mri.flip_angle mri.te mri.ti];
  err = save_mgh(permute(mri.vol,[2 1 3 4]), fspec, M, mr_parms,tempdir);  
  return;
 case {'bhdr'} %----------- BHDR ------------%
  bmri.te = mri.te;
  bmri.tr = mri.tr;
  bmri.ti = mri.ti;
  bmri.flip_angle = mri.flip_angle;
  bmri.voldim = [size(mri.vol,1) size(mri.vol,2) size(mri.vol,3)];
  bmri.nframes = size(mri.vol,4);
  bmri.T = mri.vox2ras0;
  % Recompute voxel size based on vox2ras, to assure that
  % things only depend upon vox2ras0.
  xsize = sqrt(sum(mri.vox2ras0(:,1).^2)); 
  ysize = sqrt(sum(mri.vox2ras0(:,2).^2));
  zsize = sqrt(sum(mri.vox2ras0(:,3).^2));
  bmri.volres = [xsize ysize zsize];
  outbext = 'bfloat';
  if(isfield(mri,'outbext'))
    if(strcmp(mri.outbext,'bshort')) outbext = 'bshort'; end
  end
  err = fast_svbslice(mri.vol,fstem,[],outbext,bmri);
 
 case {'nii','nii.gz'} %----------- NIFTI! ------------%
  hdr.data_type       = '';
  hdr.db_name         = '';
  hdr.extents         = 0;
  hdr.session_error   = 0;
  hdr.regular         = '';
  hdr.dim_info        = '';
  
  % Note that the order is 2 1 3 4
  hdr.dim = [mri.volsize(2) mri.volsize(1) mri.volsize(3)  mri.nframes];
  hdr.intent_p1       = 0;
  hdr.intent_p2       = 0;
  hdr.intent_p3       = 0;
  hdr.intent_code     = 0;
  
  switch(datatype)
   case 'uchar',  hdr.datatype = 2;   hdr.bitpix = 8*1;
   case 'short',  hdr.datatype = 4;   hdr.bitpix = 8*2;
   case 'int',    hdr.datatype = 8;   hdr.bitpix = 8*4;
   case 'float',  hdr.datatype = 16;  hdr.bitpix = 8*4;
   case 'double', hdr.datatype = 64;  hdr.bitpix = 8*8;
   case 'ushort', hdr.datatype = 512; hdr.bitpix = 8*2;
   case 'uint',   hdr.datatype = 768; hdr.bitpix = 8*4;
   otherwise,
    fprintf('ERROR: unrecognized data type %s\n',datatype);
    return;
  end
  hdr.slice_start     = 0;

  % volres is not permuted in MRIread()
  %hdr.pixdim          = [0 mri.volres([2 1 3]) mri.tr]; % physical units
  hdr.pixdim          = [0 mri.volres mri.tr]; % physical units
  hdr.vox_offset      = 348; % will be set again
  hdr.scl_slope       = 0;
  hdr.scl_inter       = 0;
  
  hdr.slice_end       = 0;
  
  hdr.slice_code      = 0;
  hdr.xyzt_units = bitor(2,16); % 2=mm, 16=msec

  hdr.cal_max         = max(mri.vol(:));
  hdr.cal_min         = min(mri.vol(:));
  hdr.slice_duration  = 0;
  hdr.toffset         = 0;
  hdr.glmax           = 0;
  hdr.glmin           = 0;
  hdr.descrip         = sprintf('%-80s','FreeSurfer matlab');
  hdr.aux_file        = '';
  hdr.qform_code      = 1; % 1=NIFTI_XFORM_SCANNER_ANAT
  hdr.sform_code      = 1; % 1=NIFTI_XFORM_SCANNER_ANAT
  
  % Qform (must be 6dof)
  [b,c,d,x,y,z,qfac] = vox2rasToQform(mri.vox2ras0);
  hdr.pixdim(1)       = qfac;
  hdr.quatern_b       = b;
  hdr.quatern_c       = c;
  hdr.quatern_d       = d;
  hdr.quatern_x       = x;
  hdr.quatern_y       = y;
  hdr.quatern_z       = z;
  % Sform (can by any affine)
  hdr.srow_x          = mri.vox2ras0(1,:);
  hdr.srow_y          = mri.vox2ras0(2,:);
  hdr.srow_z          = mri.vox2ras0(3,:);
  hdr.intent_name     = 'huh?';
  hdr.magic           = 'n+1';

  % Note that the order is 2 1 3 4
  hdr.vol = permute(mri.vol,[2 1 3 4]);
  err = save_nifti(hdr,fspec);
 otherwise
  fprintf('ERROR: format %s not supported\n',fmt);
  return;
end

if(err) fprintf('ERROR: saving %s \n',fstring); end

return;













function r = save_mgh(vol, fname, M, mr_parms,tempdir); 
%
% save_mgh(vol,fname, M, <mr_parms>);
%
% M is the 4x4 vox2ras transform such that
% y(i1,i2,i3), xyz = M*[i1 i2 i3 1] where the
% indicies are 0-based
%
% mr_parms = [tr flipangle te ti]
%
% See also: load_mgh, vox2ras_0to1
%
%


%
% save_mgh.m
%
% Original Author: Bruce Fischl
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

r = 1;

if(nargin < 2 || nargin > 5)
  msg = 'USAGE: save_mgh(vol,fname,M,<mr_params>,<tempdir>)';
  fprintf('%s\n',msg);
  return;
end

if exist(tempdir,'dir')==0
    error('Error in myMRIread: temporary directory does not exist')
end

if(exist('mr_parms','var')~=1), mr_parms = []; end
if(isempty(mr_parms)),   mr_parms = [0 0 0 0]; end
if(length(mr_parms) < 4)
  fprintf('ERROR: mr_parms length = %d, must be 4 or 5\n', ...
	  length(mr_parms));
  return;
end

% These dont appear to be used %
MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;
MRI_TENSOR = 6 ;

fid = fopen(fname, 'wb', 'b') ;
if(fid == -1)
  fprintf('ERROR: could not open %s for writing\n',fname);
  return;
end


[ndim1,ndim2,ndim3,frames] = size(vol) ;
fwrite(fid, 1, 'int') ;		% magic #
fwrite(fid, ndim1, 'int') ; 
fwrite(fid, ndim2, 'int') ; 
fwrite(fid, ndim3, 'int') ; 
fwrite(fid, frames, 'int') ;	% # of frames
if(ndims(vol) == 5)
  is_tensor = 1 ;
  fwrite(fid, MRI_TENSOR, 'int') ; % type = MRI_TENSOR
else
  is_tensor = 0 ;
  fwrite(fid, MRI_FLOAT, 'int') ;  % type = MRI_FLOAT
end

%%?????????????%%%
fwrite(fid, 1, 'int') ;          % dof (not used)
dof = fread(fid, 1, 'int') ; 

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

MdcD = M(1:3,1:3);
delta = sqrt(sum(MdcD.^2));

Mdc = MdcD./repmat(delta,[3 1]);
Pcrs_c = [ndim1/2 ndim2/2 ndim3/2 1]'; %'
Pxyz_c = M*Pcrs_c;
Pxyz_c = Pxyz_c(1:3);

fwrite(fid, 1,      'short') ;       % ras_good_flag = 1
fwrite(fid, delta,  'float32') ; 
fwrite(fid, Mdc,    'float32') ; 
fwrite(fid, Pxyz_c, 'float32') ; 

unused_space_size = UNUSED_SPACE_SIZE-2 ;
unused_space_size = unused_space_size - USED_SPACE_SIZE ;
fwrite(fid, zeros(unused_space_size,1), 'char') ;

fwrite(fid,vol,'float32');

fwrite(fid, mr_parms, 'float32') ; 
fclose(fid) ;

r = 0;

if (strcmpi(fname((length(fname)-3):length(fname)), '.MGZ') | ...
        strcmpi(fname((length(fname)-3):length(fname)), '.GZ'))
    
    cl=clock();
    p1=num2str(round(1e6*cl(end)));
    [~,aux]=fileparts(fname);
    aux(aux=='.')='_';
    p2=aux;
    p3=getenv('PBS_JOBID');
    cl=clock();
    rng(round(1e6*cl(end))+sum(double(fname)));
    p4=num2str(round(1000000*rand(1)));
    fname2 = [tempdir p1 '_' p2 '_' p3 '_' p4 '.mgh'];
    
    
    while exist(fname2,'file')>0
        pause(0.1*rand(1));
        cl=clock();
        p1=num2str(round(1e6*cl(end)));
        cl=clock();
        rng(round(1e6*cl(end))+sum(double(fname)));
        p4=num2str(round(1000000*rand(1)));
        fname2 = [tempdir p1 '_' p2 '_' p3 '_' p4 '.mgh'];
    end
    
    unix(sprintf('mv %s %s ; gzip %s ; mv %s.gz %s', fname, fname2, fname2, fname2, fname)) ;
    fname = fname2 ;
end
return;







function err = save_nifti(hdr,niftifile)
% err = save_nifti(hdr,niftifile)
%
% Pixel data should be in hdr.vol
%
% Handles data structures with more than 32k cols by setting
% hdr.dim(2) = -1 and hdr.glmin = ncols. This is FreeSurfer specific,
% for handling surfaces. The exception to this is when the total
% number of spatial voxels equals 163842, then the volume is 
% reshaped to 27307x1x6xnframes. This is for handling the 7th
% order icosahedron used by FS group analysis.
%

%
% save_nifti.m
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

err = 1;
if(nargin ~= 2)
  fprintf('err = save_nifti(hdr,niftifile)\n');
  return;
end

ext = niftifile((strlen(niftifile)-2):strlen(niftifile));
if(strcmpi(ext,'.gz'))
  gzip_needed = 1;
  niftifile = niftifile(1:strlen(niftifile)-3);
  %fprintf('First, saving to %s before compressing\n',niftifile);
else
  gzip_needed = 0;
end

% Check for ico7
sz = size(hdr.vol);
if(sz(1) == 163842)
  %fprintf('save_nifti: ico7 reshaping\n');
  dim = [27307 1 6 size(hdr.vol,4)];
  hdr.vol = reshape(hdr.vol, dim);
end

fp = fopen(niftifile,'w');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',niftifile);
  return;
end

hdr.data_type = [hdr.data_type(:)' repmat(' ',[1 10])];
hdr.data_type = hdr.data_type(1:10);

hdr.db_name = [hdr.db_name(:)' repmat(' ',[1 18])];
hdr.db_name = hdr.db_name(1:18);

hdr.dim = ones(1,8);
hdr.dim(1) = 4;
hdr.dim(2) = size(hdr.vol,1);
hdr.dim(3) = size(hdr.vol,2);
hdr.dim(4) = size(hdr.vol,3);
hdr.dim(5) = size(hdr.vol,4);

% This is to accomodate structures with more than 32k cols
% FreeSurfer specific. See also mriio.c.
if(hdr.dim(2) > 2^15)
  hdr.glmin = hdr.dim(2);
  hdr.dim(2) = -1;
end

hdr.pixdim = [hdr.pixdim(:)' repmat(0,[1 8])];
hdr.pixdim = hdr.pixdim(1:8);

hdr.descrip = [hdr.descrip(:)' repmat(' ',[1 80])];
hdr.descrip = hdr.descrip(1:80);

hdr.aux_file = [hdr.aux_file(:)' repmat(' ',[1 24])];
hdr.aux_file = hdr.aux_file(1:24);

hdr.intent_name = [hdr.intent_name(:)' repmat(' ',[1 16])];
hdr.intent_name = hdr.intent_name(1:16);

hdr.magic = [hdr.magic(:)' repmat(0,[1 4])];
hdr.magic = hdr.magic(1:4);

if(isempty(hdr.regular))  hdr.regular  = ' '; end
if(isempty(hdr.dim_info)) hdr.dim_info = ' '; end
if(isempty(hdr.slice_code)) hdr.slice_code = ' '; end
if(isempty(hdr.xyzt_units)) hdr.xyzt_units = ' '; end % should be err

hdr.vox_offset = 352; % not 348
fwrite(fp,348,'int');
fwrite(fp,hdr.data_type,    'char');
fwrite(fp,hdr.db_name,      'char');
fwrite(fp,hdr.extents,      'int');
fwrite(fp,hdr.session_error,'short');
fwrite(fp,hdr.regular,      'char');
fwrite(fp,hdr.dim_info,     'char');
fwrite(fp,hdr.dim,          'short');
fwrite(fp,hdr.intent_p1,    'float');
fwrite(fp,hdr.intent_p2,    'float');
fwrite(fp,hdr.intent_p3,    'float');
fwrite(fp,hdr.intent_code,  'short');
fwrite(fp,hdr.datatype,     'short');
fwrite(fp,hdr.bitpix,       'short');
fwrite(fp,hdr.slice_start,  'short');
fwrite(fp,hdr.pixdim,       'float');
fwrite(fp,hdr.vox_offset,   'float');
fwrite(fp,hdr.scl_slope,    'float');
fwrite(fp,hdr.scl_inter,    'float');
fwrite(fp,hdr.slice_end,    'short');
fwrite(fp,hdr.slice_code,   'char');
fwrite(fp,hdr.xyzt_units,   'char');
fwrite(fp,hdr.cal_max,      'float');
fwrite(fp,hdr.cal_min,      'float');
fwrite(fp,hdr.slice_duration,'float');
fwrite(fp,hdr.toffset,       'float');
fwrite(fp,hdr.glmax,         'int');
fwrite(fp,hdr.glmin,         'int');
fwrite(fp,hdr.descrip,       'char');
fwrite(fp,hdr.aux_file,      'char');
fwrite(fp,hdr.qform_code,    'short');
fwrite(fp,hdr.sform_code,    'short');
fwrite(fp,hdr.quatern_b,     'float');
fwrite(fp,hdr.quatern_c,     'float');
fwrite(fp,hdr.quatern_d,     'float');
fwrite(fp,hdr.quatern_x,     'float');
fwrite(fp,hdr.quatern_y,     'float');
fwrite(fp,hdr.quatern_z,     'float');
fwrite(fp,hdr.srow_x,        'float');
fwrite(fp,hdr.srow_y,        'float');
fwrite(fp,hdr.srow_z,        'float');
fwrite(fp,hdr.intent_name,   'char');
fwrite(fp,hdr.magic,         'char');

% Pad to get to 352 bytes (header size is 348)
fwrite(fp,0,'char');
fwrite(fp,0,'char');
fwrite(fp,0,'char');
fwrite(fp,0,'char');

npix = prod(size(hdr.vol));
switch(hdr.datatype)
 case   2, nitemswritten = fwrite(fp,hdr.vol,'uchar'); % dont use char
 case   4, nitemswritten = fwrite(fp,hdr.vol,'short');
 case   8, nitemswritten = fwrite(fp,hdr.vol,'int');
 case  16, nitemswritten = fwrite(fp,hdr.vol,'float');
 case  64, nitemswritten = fwrite(fp,hdr.vol,'double');
 case 512, nitemswritten = fwrite(fp,hdr.vol,'ushort');
 case 768, nitemswritten = fwrite(fp,hdr.vol,'uint');
 otherwise,
  fprintf('WARNING: data type %d not supported, but writing as float',...
	  hdr.datatype);
  nitemswritten = fwrite(fp,hdr.vol,'float');
  return;
end
fclose(fp);

if(npix ~= nitemswritten)
  fprintf('ERROR: tried to write %d, but only wrote %d',npix,nitemswritten);
  return;
end

if(gzip_needed)
  cmd = sprintf('gzip -f %s', niftifile);
  %fprintf('Compressing with\n');
  %fprintf('   %s\n',cmd);
  unix(cmd);
end


err = 0;
return;







