function v = fast_ldbvoxel(stem,c,r,s,base)
% v = fast_ldbvoxel(stem,c,r,s,base)
%
% Loads a voxel value from a bshort/bfloat volume. If there
% is more than one frame, all frames are read in.
% crs are assumed to be zero-based, unless base is set to 1.
%
% $Id: fast_ldbvoxel.m,v 1.3 2003/09/28 21:46:58 greve Exp $

v=[];

if(nargin ~= 4 & nargin ~= 5)
  fprintf('v = fast_ldbvoxel(stem,c,r,s,<base>)\n');
  return;
end

if(nargin == 4) base = 0; end
if(base == 1)
  c = c - 1;
  r = r - 1;
  s = s - 1;
end

[nrows ncols nframes fs nslices endian bext] = fmri_bfiledim(stem);
if(isempty(nrows))
  fprintf('ERROR: fast_ldbvoxel: %s\n',stem);
  return;
end

if(c < 0 | c >= ncols)
  fprintf('ERROR: fast_ldbvoxel: column %d is out of range\n',c);
  return;
end
if(r < 0 | r >= nrows)
  fprintf('ERROR: fast_ldbvoxel: row %d is out of range (%d,%d)\n',r,1,nrows);
  return;
end
if(s < 0 | s >= nslices)
  fprintf('ERROR: fast_ldbvoxel: slice %d is out of range\n',s);
  return;
end

if( strcmp(bext,'bshort'))
  precision = 'int16';
  nbytes = 2;
else               
  precision = 'float32';
  nbytes = 4;
end

bfile = sprintf('%s_%03d.%s',stem,s,bext);

%%%% Open the bfile %%%%%
if(endian == 0) fid=fopen(bfile,'r','b'); % Big-Endian
else            fid=fopen(bfile,'r','l'); % Little-Endian
end
if(fid == -1)
  fprintf('ERROR: could not open %s\n',bfile);
  return;
end

baseoffset = r*ncols + c;
nv = nrows*ncols;

v = zeros(nframes,1);
for frame = 0:nframes-1
  offset = nbytes*(baseoffset + nv*frame);

  status = fseek(fid,offset,'bof');
  if(status == -1)
    fprintf('ERROR: error seeking to %d in %s\n',offset,bfile);
    fclose(fid); v = [];
    return;
  end

  [v(frame+1) count] = fread(fid,1,precision);
  if(count ~= 1)
    fprintf('ERROR: error reading %d frame from %s\n',frame+1,bfile);
    fclose(fid); v = [];
    return;
  end
end
fclose(fid); 

return;

