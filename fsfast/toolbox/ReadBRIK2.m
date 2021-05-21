function V = ReadBRIK2(fname)
%
% V = ReadBRIK2(fname)
%
% ReadBRIK function to read 3D+time AFNI BRIKs
% of any type and dimensionality
%
% Example:
% brikfname = 'ASt3avvr+orig.BRIK'
% BRIKDATA = ReadBRIK2(brikfname)
% BRIKDATA is now a 4D matrix (xdim,ydim,zdim,timepoint)
%
% Written December 4, 1998 by Timothy M. Ellmore
% Laboratory of Brain and Cognition, NIMH
%
% Modified by Douglas N. Greve
% greve@nmr.mgh.harvard.edu
% MGH-NMR Center, Boston Mass
% Now does not require any info other than the brikname,
% ie, automatically detects size, type, endianness, etc


%
% ReadBRIK2.m
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

V = [];

if(nargin ~= 1)
  fprintf('USAGE: V = ReadBRIK2(fname)\n');
  return;
end

% Strip off .BRIK
l = length(fname);
if(l < 5)
  fprintf('ERROR: fname = %s, must be something.BRIK\n',fname);
  return;
end
stem = fname(1:l-5);

headerfile = sprintf('%s.HEAD',stem);
header = ReadAFNIHead2(headerfile);
if(isempty(header)) return; end

byteorder = header.byteorder;
type = header.type;
xdim = header.xdim;
ydim = header.ydim;
zdim = header.zdim;
tdim = header.nvols;

% set permission string 
permission = 'r';
machine_type = sprintf('%s', header.byteorder);
%permission = sprintf('r%s',header.byteorder);

disp('ReadBRIK2: reading raw data . . . .')
%fid = fopen(fname, permission);
fid = fopen(fname, permission,machine_type); 
if(fid == -1)
  fprintf('ERROR: cannot open %s\n',fname);
  return;
end

[V count] = fread(fid, [xdim * ydim * zdim * tdim], header.type);
fclose(fid);
disp('ReadBRIK2: done reading raw data')

if(count ~= xdim * ydim * zdim * tdim)
  fprintf('ERROR: read %d, expected %d\n',count,xdim * ydim * zdim * tdim);
  V = [];
  return;
end

% Reshape to cols, rows, slices, frames
V = reshape(V,xdim,ydim,zdim,tdim);

return;
