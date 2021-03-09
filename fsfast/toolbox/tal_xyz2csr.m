function [csr] = tal_xyz2csr(x, y, z, res, fov)
%
% [csr] = tal_xyz2csr(x, y, z, res, fov)
%
% Returns the corronal slice number and the row and column numbers 
% in the corronal slice given the Talairach coords x,y,z (mm).
%
% res is the voxel resolution in mm (defalut 8mm)
% fov is the field of view (default 256mm)


%
% tal_xyz2csr.m
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
