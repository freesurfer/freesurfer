function [vol, M, mr_parms] = load_mgh2(fname)
fprintf('This program, load_mgh2, is obsolete. Use MRIread instead\n');
return;

% [vol, M, mr_parms] = load_mgh2(fname)
%
% M is the 4x4 vox2ras transform such that
% y(i1,i2,i3), xyz1 = M*[i1 i2 i3 1] where the
% indicies are 0-based. If the input has multiple frames,
% only the first frame is read.
%
% mr_parms = [tr flipangle te ti]
%
% See also: save_mgh, vox2ras_0to1
%


%
% load_mgh2.m
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


vol = [];
M = [];
mr_parms = [];

if(nargin < 1 | nargin > 1)
  msg = 'USAGE: [vol M] = load_mgh(fname)';
  return;
end

fid    = fopen(fname, 'rb', 'b') ;
if(fid == -1)
  fprintf('ERROR: could not open %s for reading\n',fname);
  return;
end
v       = fread(fid, 1, 'int') ; 
ndim1   = fread(fid, 1, 'int') ; 
ndim2   = fread(fid, 1, 'int') ; 
ndim3   = fread(fid, 1, 'int') ; 
nframes = fread(fid, 1, 'int') ; 
type    = fread(fid, 1, 'int') ; 
dof     = fread(fid, 1, 'int') ; 

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

unused_space_size = UNUSED_SPACE_SIZE-2 ;
ras_good_flag = fread(fid, 1, 'short') ; 
if (ras_good_flag)
  delta  = fread(fid, 3, 'float32') ; 
  Mdc    = fread(fid, 9, 'float32') ; 
  Mdc    = reshape(Mdc,[3 3]);
  Pxyz_c = fread(fid, 3, 'float32') ; 

  D = diag(delta);

  Pcrs_c = [ndim1/2 ndim2/2 ndim3/2]'; %'

  Pxyz_0 = Pxyz_c - Mdc*D*Pcrs_c;

  M = [Mdc*D Pxyz_0;  ...
       0 0 0 1];
  unused_space_size = unused_space_size - USED_SPACE_SIZE ;
end

MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;

fseek(fid, unused_space_size, 'cof') ;

nv = ndim1 * ndim2 * ndim3;  

switch type
  case MRI_FLOAT,
    vol = fread(fid, nv, 'float32') ; 
  case MRI_UCHAR,
    vol = fread(fid, nv, 'uchar') ; 
  case MRI_SHORT,
    vol = fread(fid, nv, 'short') ; 
  case MRI_INT,
    vol = fread(fid, nv, 'int') ; 
end
if(~feof(fid))
  [mr_parms count] = fread(fid,4,'float32');
  if(count ~= 4) 
    fprintf('WARNING: error reading MR params\n');
  end
end
fclose(fid) ;

nread = prod(size(vol));
if(nread ~= nv)
  fprintf('ERROR: tried to read %d, actually read %d\n',nv,nread);
  vol = [];
  return;
end
  
vol = reshape(vol,[ndim1 ndim2 ndim3]);


return;
