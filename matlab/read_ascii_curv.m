function [curv] = read_curv(fname)
%
% reads an ascii curvature into a vector
%


%
% read_ascii_curv.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
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


fid = fopen(fname, 'r') ;
nvertices = fscanf(fid, '%d', 1);
all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
curv = all(5, :)' ;
