function [nrows, ncols, ntp, fs, ns, endian, bext] = fmri_bfiledim(stem)
% [nrows ncols ntp fs ns endian bext] = fmri_bfiledim(stem)
%
% $Id: fmri_bfiledim.m,v 1.1 2003/03/04 20:47:39 greve Exp $

if(nargin ~= 1) 
  msg = '[nrows ncols ntp fs ns endian bext] = fmri_bfiledim(stem)';
  qoe(msg); error(msg);
end

stem = deblank(stem);

fs = 0;
fname = sprintf('%s_%03d.hdr',stem,fs);
fid = fopen(fname,'r');
while(fid == -1 & fs < 100)
  fs = fs + 1;
  fname = sprintf('%s_%03d.hdr',stem,fs);
  fid = fopen(fname,'r');
end

if(fs == 100)
  msg = sprintf('Could not find any slices with stem %s',stem);
  qoe(msg); error(msg);
end

hdr = fscanf(fid,'%d',[1,4]);
fclose(fid);
nrows  = hdr(1);
ncols  = hdr(2);
ntp    = hdr(3);
endian = hdr(4);

ns = 1;
slice = fs + ns;
fname = sprintf('%s_%03d.hdr',stem,slice);
fid = fopen(fname,'r');
while(fid ~= -1)
  fclose(fid);
  ns = ns + 1;
  slice = fs + ns;
  fname = sprintf('%s_%03d.hdr',stem,slice);
  fid = fopen(fname,'r');
end


% Get extension %
bext = 'bshort';
fname = sprintf('%s_%03d.%s',stem,fs,bext);
fid = fopen(fname,'r');
if(fid == -1)
  bext = 'bfloat';
  fname = sprintf('%s_%03d.%s',stem,fs,bext);
  fid = fopen(fname,'r');
  if(fid == -1)
     msg = sprintf('Could not find bshort or bfloat for %s',stem);
     qoe(msg);error(msg);
  end
end

fclose(fid);

return;
