%
% read_patch.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/01/22 20:59:09 $
%    $Revision: 1.3.2.2 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


function patch = read_patch(fname)
% function patch = read_patch(fname)


fid = fopen(fname,'r');
if (fid == -1)
   error('could not open file %s', fname) ;
end

ver = fread(fid, 1, 'int', 0, 'b');
if (ver ~= -1)
   error('incorrect version # %d (not -1) found in file',ver) ;
end

patch.npts = fread(fid, 1, 'int', 0, 'b') ;

for i=1:patch.npts
    ind = fread(fid, 1, 'int', 0, 'b') ;
    if (ind < 0)
       ind = -ind - 1 ;
    else
       ind = ind - 1 ;
    end
    patch.ind(i) = ind ;
    patch.x(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.y(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.z(i) = fread(fid, 1, 'float', 0, 'b') ;
    patch.vno(i) = ind ;
end

fclose(fid);
