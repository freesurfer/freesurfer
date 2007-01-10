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
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%


fid = fopen(fname, 'w') ;
nvertices = size(curv,1) ;
fprintf(fid, '%d\n', nvertices) ;
for i=1:nvertices
		fprintf(fid, '0 0.0 0.0 0.0 %f\n', curv(i,1)) ;
end
fclose(fid) ;
