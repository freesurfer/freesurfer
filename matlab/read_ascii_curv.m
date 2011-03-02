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
%    $Date: 2011/03/02 00:04:12 $
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


fid = fopen(fname, 'r') ;
nvertices = fscanf(fid, '%d', 1);
all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
curv = all(5, :)' ;
