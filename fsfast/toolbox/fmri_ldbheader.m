function [nrows, ncols, ntp, endian, fs] = fmri_ldbheader(stem)
% [nrows ncols ntp endian fs ns] = fmri_ldbheader(stem)


%
% fmri_ldbheader.m
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

if(nargin ~= 1) 
  msg = '[nrows ncols ntp endian fs] = fmri_ldbheader(stem)';
  qoe(msg); error(msg);
end


nslices = 0;
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
nrows  = hdr(1);
ncols  = hdr(2);
ntp    = hdr(3);
endian = hdr(4);


return;
