function [xyz] = tal_csr2xyz(col, slice, row, res, fov)
%
% [xyz] = tal_csr2xyz(col, slice, row, res, fov)
%
% Returns Talairach coords x,y,z (mm) given the slice 
% number and the row and column in the corronal slice.  
%
% res is the voxel resolution in mm (defalut 8mm)
% fov is the field of view (default 256mm)


%
% tal_csr2xyz.m
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
  msg = 'USAGE: [xyz] = tal_csr2xyz(col, slice, row, <res, <fov>>)';
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

x = res*(   col - col_offset);
y = res*( slice - slice_offset);
z = res*(  -row + row_offset);

xyz(1) = x;
xyz(2) = y;
xyz(3) = z;


return
