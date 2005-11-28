function [nslices, nrows, ncols, nt, endian, bext, hdrdat] = fmri_bvoldim(stem)
% [nslices nrows ncols nt endian bext hdrdat] = fmri_bvoldim(stem)
% 
% $Id: fmri_bvoldim.m,v 1.2 2005/11/28 22:40:29 greve Exp $

nslices = 0;
nrows   = 0;
ncols   = 0;
nt      = 0;
endian  = -1;
bext    = '';
hdrdat  = [];

if(nargin ~= 1)
  msg = 'USAGE: [nslices nrows ncols nt] = fmri_bvoldim(stem)';
  qoe(msg); error(msg);
end

stem = deblank(stem);

firstslice = -1;
nslices = 0;
for slice = 0:1000
  fname = sprintf('%s_%03d.hdr',stem,nslices);
  fid = fopen(fname,'r');
  if(fid == -1) break;  end
  hdr = fscanf(fid,'%d',[1,4]);
  if(isempty(hdr))
    fprintf('ERROR: reading %s\n',fname);
    nslices = 0;
    return;
  end
  nrows = hdr(1);
  ncols = hdr(2);
  nt = hdr(3);
  endian = hdr(4);
  fclose(fid);

  % Get the extension %
  if(firstslice == -1) 
    firstslice = nslices ;
    bext = 'bshort';
    fname = sprintf('%s_%03d.%s',stem,nslices,bext);
    fid = fopen(fname,'r');
    if(fid == -1) 
      bext = 'bfloat';
      fname = sprintf('%s_%03d.%s',stem,nslices,bext);
      fid = fopen(fname,'r');
      if(fid == -1) 
        nslices = 0;
        nrows   = 0;
        ncols   = 0;
        nt      = 0;
        bext = '';
        endian = -1;
        return;
      end
    end
    fclose(fid);
  end

  nslices = nslices + 1;
end

fname = sprintf('%s.dat',stem);
fid = fopen(fname,'r');
if(fid ~= -1)
  fclose(fid);
  hdrdat = fmri_lddat3(fname);
end

return;
