function [csr] = tal_xyz2csr(x, y, z, res, fov)
%
% [csr] = tal_xyz2csr(x, y, z, res, fov)
%
% Returns the corronal slice number and the row and column numbers 
% in the corronal slice given the Talairach coords x,y,z (mm).
%
% res is the voxel resolution in mm (defalut 8mm)
% fov is the field of view (default 256mm)

if(nargin < 3 | nargin > 6)
  msg = 'USAGE: [xyz] = [csr] = tal_xyz2csr(x, y, z, <res, <fov>>)';
  qoe(msg);error(msg);
end

if(nargin == 3) 
  res =   8; 
  fov = 256;
elseif(nargin == 4) 
  fov = 256;
end

col_offset   = fov/(2*res) + 0.5;
slice_offset = fov/(2*res) + 0.5 - 1.0;
row_offset   = fov/(2*res) + 0.5;

col   =  -x/res + col_offset ;
slice =  y/res + slice_offset ;
row   = -z/res + row_offset ;

csr(1) = col;
csr(2) = slice;
csr(3) = row;

return
