function [vol, ext, endian] = fmri_ldbvolume(stem,ext)
% [vol ext endian] = fmri_ldbvolume(stem,<ext>)
%
% Loads a volume in bfile format and returns a 4D structure
% of dimension Nslices X Nrows X Ncols X Ndepth. Automatically
% counts the number of slices and determines the extension 
% (unless extension is specified).
%


%
% fmri_ldbvolume.m
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
ext = [];
endian = [];

if(nargin ~= 1 & nargin ~= 2)  
  fprintf(2,'USAGE: fmri_ldbvolume(stem,<ext>)\n');
  qoe; error; return;
end

% First, check whether stem is really an entire filename %
l = length(stem);
if(l > 6) 
  tmp = stem(l-6:l);
  i = findstr(tmp,'bshort');
  j = findstr(tmp,'bfloat');
  if(~isempty(i) | ~isempty(j))
    fprintf('INFO: loading %s as a simple bfile\n',stem);
    vol = fmri_ldbfile(stem);
    return;
  end
end

firstslice = 0;
nslices = 0;
slice = firstslice;
fname = sprintf('%s_%03d.hdr',stem,slice);
fid = fopen(fname,'r');
while(fid ~= -1)
  hdr = fscanf(fid,'%d',[1,4]);
  Nrows  = hdr(1);
  Ncols  = hdr(2);
  Ndepth = hdr(3);
  endian = hdr(4);
  fclose(fid);
  nslices = nslices + 1;
  slice   = slice + 1;
  fname = sprintf('%s_%03d.hdr',stem,slice);
  fid = fopen(fname,'r');
end

if(nslices == 0)
  fprintf('ERROR: cannot find volume matching %s\n',stem);
  return;
end

% fprintf('fmri_ldbvolume: found %d slices\n',nslices);

if(nargin == 1)
  ext = 'bshort';
  fname = sprintf('%s_%03d.%s',stem,firstslice,ext);
  fid = fopen(fname,'r');
  if(fid == -1)
    ext = 'bfloat';
    fname = sprintf('%s_%03d.%s',stem,firstslice,ext);
    fid = fopen(fname,'r');
    if(~fid)
      msg = sprintf('fmri_ldbvolume: could not determine extension for stem',...
                    stem);
      qoe(msg); error(msg);
    end
  end
  fclose(fid);
end
  
vol = zeros(nslices,Nrows,Ncols,Ndepth);
for slice = firstslice : firstslice + nslices - 1
  % fprintf('Loading Slice %3d\n',slice);
  fname = sprintf('%s_%03d.%s',stem,slice,ext);
  n = slice-firstslice+1;
  z = fmri_ldbfile(fname);
  if(size(z,1) == 1 & size(z,3) == 1)
    vol(n,:,:)  = z;
  else
    vol(n,:,:,:) = z;
  end

end

return;
