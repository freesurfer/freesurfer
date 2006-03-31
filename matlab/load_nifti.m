function hdr = load_nifti(niftifile,hdronly)
% hdr = load_nifti(niftifile,hdronly)
%
% Loads nifti header and volume. The volume is stored
% in hdr.vol. Columns and rows are not swapped.
%
% Handles compressed nifti (nii.gz) by issuing a unix command to
% uncompress the file to a temporary file, which is then deleted.
%
% Dimensions are in mm and msec
% hdr.pixdim(1) = physical size of first dim (eg, 3.125 mm or 2000 ms)
% hdr.pixdim(2) = ...
% 
% hdr.vox2ras is the vox2ras matrix based on sform
%
% See also: load_nifti_hdr.m
%
% $Id: load_nifti.m,v 1.3 2006/03/31 06:24:12 greve Exp $

hdr = [];

if(nargin < 1 | nargin > 2)
  fprintf('hdr = load_nifti(niftifile,<hdronly>)\n');
  return;
end

if(~exist('hdronly','var')) hdronly = []; end
if(isempty(hdronly)) hdronly = 0; end

% unzip if it is compressed 
ext = niftifile((strlen(niftifile)-2):strlen(niftifile));
if(strcmpi(ext,'.gz'))
  gzipped =  round(rand(1)*10000000);
  ind = findstr(niftifile, '.');
  new_niftifile = sprintf('/tmp/tmp%d.nii', gzipped);
  fprintf('Uncompressing %s to %s\n',niftifile,new_niftifile);
  unix(sprintf('zcat %s > %s', niftifile, new_niftifile)) ;
  niftifile = new_niftifile ;
else
  gzipped = -1 ;
end

hdr = load_nifti_hdr(niftifile);
if(isempty(hdr)) 
  if(gzipped >=0) unix(sprintf('rm %s', niftifile)); end
  return; 
end

% If only header is desired, return now
if(hdronly) 
  if(gzipped >=0) unix(sprintf('rm %s', niftifile)); end
  return; 
end

% Open to read the pixel data
fp = fopen(niftifile,'r',hdr.endian);

% Get past the header
fseek(fp,round(hdr.vox_offset),'bof');

switch(hdr.datatype)
 case   2, [hdr.vol nitemsread] = fread(fp,inf,'char');
 case   4, [hdr.vol nitemsread] = fread(fp,inf,'short');
 case   8, [hdr.vol nitemsread] = fread(fp,inf,'int');
 case  16, [hdr.vol nitemsread] = fread(fp,inf,'float');
 case  64, [hdr.vol nitemsread] = fread(fp,inf,'double');
 case 512, [hdr.vol nitemsread] = fread(fp,inf,'ushort');
 case 768, [hdr.vol nitemsread] = fread(fp,inf,'uint');
 otherwise,
  fprintf('ERROR: data type %d not supported',hdr.datatype);
  hdr = [];
  return;
end

fclose(fp);
if(gzipped >=0) unix(sprintf('rm %s', niftifile)); end

% Get total number of voxels
dim = hdr.dim(2:end);
ind0 = find(dim==0);
dim(ind0) = 1;
nvoxels = prod(dim);

% Check that that many voxels were read in
if(nitemsread ~= nvoxels) 
  fprintf('ERROR: %s, read in %d voxels, expected %d\n',...
	  niftifile,nitemsread,nvoxels);
  hdr = [];
  return;
end

hdr.vol = reshape(hdr.vol, dim');
if(hdr.scl_slope ~= 0)
  fprintf('nifti rescale: slope = %g, intercept = %g\n',...
	  hdr.scl_slope,hdr.inter);
  fprintf('Good luck, this has never been tested ... \n');
  hdr.vol = hdr.vol * hdr.scl_slope  + hdr.scl_inter;
end

return;





