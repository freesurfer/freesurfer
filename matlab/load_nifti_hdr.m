function hdr = load_nifti_hdr(niftifile)
% hdr = load_nifti_hdr(niftifile)
%
% Changes units to mm and msec.
% Creates hdr.vox2ras based on sform.
% Does not use qform (yet).
% Does not and will not handle compressed. Compression is handled
% in load_nifti.m, which calls load_nifti_hdr.m after any
% decompression. 
%
% Endianness is returned as hdr.endian, which is either 'l' or 'b'. 
% When opening again, use fp = fopen(niftifile,'r',hdr.endian);
%
% $Id: load_nifti_hdr.m,v 1.2 2006/03/30 07:52:48 greve Exp $

hdr = [];

if(nargin ~= 1)
  fprintf('hdr = load_nifti_hdr(niftifile)\n');
  return;
end

% Try opening as big endian first
fp = fopen(niftifile,'r','b');
if(fp == -1) 
  niftifile0 = niftifile;
  niftifile = sprintf('%s.hdr',niftifile);
  fp = fopen(niftifile,'rb');
  if(fp == -1) 
    fprintf('ERROR: could not read %s or %s\n',niftifile0,niftifile);
    return;
  end
end

hdr.sizeof_hdr  = fread(fp,1,'int');
if(hdr.sizeof_hdr ~= 348)
  fclose(fp);
  % Now try opening as little endian
  fp = fopen(niftifile,'r','l');
  hdr.sizeof_hdr  = fread(fp,1,'int');
  if(hdr.sizeof_hdr ~= 348)
    fclose(fp);
    fprintf('ERROR: %s: hdr size = %d, should be 348\n',...
	    niftifile,hdr.sizeof_hdr);
    hdr = [];
    return;
  end
  hdr.endian = 'l';
else
  hdr.endian = 'b';
end

hdr.data_type       = fscanf(fp,'%c',10);
hdr.db_name         = fread(fp,18,'char');
hdr.extents         = fread(fp, 1,'int');
hdr.session_error   = fread(fp, 1,'short');
hdr.regular         = fread(fp, 1,'char');
hdr.dim_info        = fread(fp, 1,'char');

hdr.dim             = fread(fp, 8,'short');
hdr.intent_p1       = fread(fp, 1,'float');
hdr.intent_p2       = fread(fp, 1,'float');
hdr.intent_p3       = fread(fp, 1,'float');
hdr.intent_code     = fread(fp, 1,'short');

hdr.datatype        = fread(fp, 1,'short');
hdr.bitpix          = fread(fp, 1,'short');
hdr.slice_start     = fread(fp, 1,'short');

hdr.pixdim          = fread(fp, 8,'float'); % physical units
hdr.vox_offset      = fread(fp, 1,'float');
hdr.scl_slope       = fread(fp, 1,'float');
hdr.scl_inter       = fread(fp, 1,'float');

hdr.slice_end       = fread(fp, 1,'short');

hdr.slice_code      = fread(fp, 1,'char');
hdr.xyzt_units      = fread(fp, 1,'char');
hdr.cal_max         = fread(fp, 1,'float');
hdr.cal_min         = fread(fp, 1,'float');
hdr.slice_duration  = fread(fp, 1,'float');
hdr.toffset         = fread(fp, 1,'float');
hdr.glmax           = fread(fp, 1,'int');
hdr.glmin           = fread(fp, 1,'int');
hdr.descrip         = fscanf(fp,'%c',80);
hdr.aux_file        = fscanf(fp,'%c',24);
hdr.qform_code      = fread(fp, 1,'short');
hdr.sform_code      = fread(fp, 1,'short');

hdr.quatern_b       = fread(fp, 1,'float');
hdr.quatern_c       = fread(fp, 1,'float');
hdr.quatern_d       = fread(fp, 1,'float');
hdr.quatern_x       = fread(fp, 1,'float');
hdr.quatern_y       = fread(fp, 1,'float');
hdr.quatern_z       = fread(fp, 1,'float');
hdr.srow_x          = fread(fp, 4,'float');
hdr.srow_y          = fread(fp, 4,'float');
hdr.srow_z          = fread(fp, 4,'float');
hdr.intent_name     = fscanf(fp,'%c',16);
hdr.magic           = fscanf(fp,'%c',4);

fclose(fp);

% look at xyz units and convert to mm if needed
xyzunits = bitand(hdr.xyzt_units,7); % 0x7
switch(xyzunits)
  case 1, xyzscale = 1000.000; % meters
  case 2, xyzscale =    1.000; % mm
  case 3, xyzscale =     .001; % microns
end
hdr.pixdim(2:4) = hdr.pixdim(2:4) * xyzscale;
hdr.srow_x = hdr.srow_x * xyzscale;
hdr.srow_y = hdr.srow_y * xyzscale;
hdr.srow_z = hdr.srow_z * xyzscale;

% look at time units and convert to msec if needed
tunits = bitand(hdr.xyzt_units,3*16+8); % 0x38 
switch(tunits)
  case  8, tscale = 1000.000; % seconds
  case 16, tscale =    1.000; % msec
  case 32, tscale =     .001; % microsec
end
hdr.pixdim(5) = hdr.pixdim(5) * tscale;

% Change value in xyzt_units to reflect scale change
hdr.xyzt_units = bitor(2,16); % 2=mm, 16=msec


% should look at slc slope and interp and rescale

if(hdr.sform_code ~= 0)
  hdr.vox2ras = [hdr.srow_x'; 
		 hdr.srow_y'; 
		 hdr.srow_z';
		0 0 0 1];
% Else check qform

end

return;





