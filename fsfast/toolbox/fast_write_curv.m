function [curv] = write_curv(fname, curv, fnum)
%
% [curv] = write_curv(fname, curv, fnum)
%
% writes a curvature vector into a binary file
% fname - name of file to write to
% curv  - vector of curvatures
% fnum  - # of faces in surface.
%


%
% fast_write_curv.m
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

if(nargin ~= 3)
  fprintf('USAGE: curv = write_curv(fname, curv, fnum)\n');
  return;
end

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;
vnum = length(curv) ;
NEW_VERSION_MAGIC_NUMBER = 16777215;
fast_fwrite3(fid, NEW_VERSION_MAGIC_NUMBER ) ;
fwrite(fid, vnum,'int32') ;
fwrite(fid, fnum,'int32') ;
fwrite(fid, 1, 'int32');
fwrite(fid, curv, 'float') ;
fclose(fid) ;

