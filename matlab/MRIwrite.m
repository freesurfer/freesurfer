function err = MRIwrite(mri,fstring,datatype)
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
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2008/11/18 20:49:52 $
%    $Revision: 1.11 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%


err = 1;

if(nargin < 2 | nargin > 3)
  fprintf('err = MRIwrite(mri,fstring,<datatype>)\n');
  return;
end

if(nargin == 2) datatype = 'float'; end

[fspec fstem fmt] = MRIfspec(fstring,0); % 0 = turn off checkdisk
if(isempty(fspec))
  fprintf('ERROR: could not determine format of %s\n',fstring);
  return;
end

switch(fmt)
 case {'mgh','mgz'} %----------- MGH/MGZ ------------%
  M = mri.vox2ras0;
  mr_parms = [mri.tr mri.flip_angle mri.te mri.ti];
  err = save_mgh(permute(mri.vol,[2 1 3 4]), fspec, M, mr_parms);  
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
  bmri.volres = [mri.xsize mri.ysize mri.zsize];
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











