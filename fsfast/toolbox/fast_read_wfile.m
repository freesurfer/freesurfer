function w = fast_read_wfile(fname)
%
% w = fast_read_wfile(fname)
% reads a vector from a binary 'w' file
%	fname - name of file to read from
%	w     - vector of values 
%


%
% fast_read_wfile.m
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

w = [];
v = [];

if(nargin ~= 1)
  fprintf('USAGE: [w,v] = fast_read_wfile(fname, w) \n');
  return;
end

% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
  str = sprintf('could not open w file %s.', fname) ;
  error(str) ;
end

fread(fid, 1, 'int16') ;  % Skip ilat
vnum = fast_fread3(fid) ; % Number of non-zero values
v = zeros(vnum,1) ;
w0 = zeros(vnum,1) ;
for i=1:vnum
  v(i) = fast_fread3(fid) ;
  w0(i) = fread(fid, 1, 'float') ;
end

fclose(fid) ;

w = zeros(max(v),1);
w(v+1) = w0;

return;




