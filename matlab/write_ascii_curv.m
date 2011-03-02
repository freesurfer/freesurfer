function [fid] = write_curv(curv, fname)
%
% writes a curvature vector into an ascii file
%


%
% write_ascii_curv.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.3 $
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


fid = fopen(fname, 'w') ;
nvertices = size(curv,1) ;
fprintf(fid, '%d\n', nvertices) ;
for i=1:nvertices
		fprintf(fid, '0 0.0 0.0 0.0 %f\n', curv(i,1)) ;
end
fclose(fid) ;
