function hdr = load_nifti(niftifile,hdronly)
% hdr = load_nifti(niftifile,hdronly)
%
% Loads nifti header and volume. The volume is stored
% in hdr.vol. It is not permuted to swap columns and rows.
%
% hdr.pixdim(1) = physical size of first dim (eg, 3.125 mm or 2000 ms)
% hdr.pixdim(2) = ...
% 
% hdr.vox2ras is the vox2ras matrix based on sform
%
% See also: load_nifti_hdr.m
%
% $Id: load_nifti.m,v 1.1 2006/03/30 07:00:27 greve Exp $

hdr = [];

if(nargin < 1 | nargin >2)
  fprintf('hdr = load_nifti(niftifile,<hdronly>)\n');
  return;
end

if(~exist('hdronly','var')) hdronly = []; end
if(isempty(hdronly)) hdronly = 0; end

hdr = load_nifti_hdr(niftifile);
if(isempty(hdr)) return; end

% If only header is desired, return now
if(hdronly) return; end

% Open to read the pixel data
fp = fopen(niftifile,'r',hdr.endian);

% Get past the header
fseek(fp,hdr.sizeof_hdr,'bof');

switch(hdr.datatype)
 case   2, [hdr.vol nitemsread] = fread(fp,inf,'char');
 case   4, [hdr.vol nitemsread] = fread(fp,inf,'short');
 case   8, [hdr.vol nitemsread] = fread(fp,inf,'int');
 case  16, [hdr.vol nitemsread] = fread(fp,inf,'float');
 case  64, [hdr.vol nitemsread] = fread(fp,inf,'double');
 case 512, [hdr.vol nitemsread] = fread(fp,inf,'ushort');
 case 768, [hdr.vol nitemsread] = fread(fp,inf,'uint');
 otherwise,
  fprintf('ERROR: data type %d not supported',hdr.dime.datatype);
  hdr = [];
  return;
end

fclose(fp);

% Get total number of voxels
dim = hdr.dim(2:end);
ind0 = find(dim==0);
dim(ind0) = 1;
nvoxels = prod(dim);

% Check that that many voxels were read in
if(nitemsread ~= nvoxels) 
  fprintf('ERROR: %s, read in %d voxels, expected %d\n',...
	  imgfile,nitemsread,nvoxels);
  hdr = [];
  return;
end

hdr.vol = reshape(hdr.vol, dim');

return;





