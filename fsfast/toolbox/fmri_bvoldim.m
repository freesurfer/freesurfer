function [nslices, nrows, ncols, nt, endian, bext, hdrdat] = fmri_bvoldim(stem)
% [nslices nrows ncols nt endian bext hdrdat] = fmri_bvoldim(stem)
% 
%


%
% fmri_bvoldim.m
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

% First try to get all the info from the bhdr. This circumvents the
% problem of having extra slices laying around from some larger 
% volume. The extra slices do not get overwritten when a smaller
% volume is created.
bhdr = sprintf('%s.bhdr',stem);
if(fast_fileexists(bhdr))
  mri = fast_ldbhdr(bhdr);
  if(isempty(mri)) return; end
  ncols   = mri.voldim(1);
  nrows   = mri.voldim(2);
  nslices = mri.voldim(3);
  bext = 'bshort';
  fname = sprintf('%s_000.%s',stem,bext);  
  if(~fast_fileexists(fname))
    bext = 'bfloat';
    fname = sprintf('%s_000.%s',stem,bext);  
    if(~fast_fileexists(fname))
      fprintf('ERROR: cannot find slice 000 for %s\n',stem);
      nslices = 0;
      nrows   = 0;
      ncols   = 0;
      nt      = 0;
      return;
    end
  end
  fname = sprintf('%s_000.hdr',stem);
  fid = fopen(fname,'r');
  if(fid == -1) 
    nslices = 0;
    fprintf('ERROR: cannot find %s\n',fname);
  end
  hdr = fscanf(fid,'%d',[1,4]);
  if(isempty(hdr))
    fprintf('ERROR: reading %s\n',fname);
    nslices = 0;
    return;
  end
  nt = hdr(3);
  endian = hdr(4);
  fclose(fid);
  return;
end

% If it gets here, there is not bhdr file
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
