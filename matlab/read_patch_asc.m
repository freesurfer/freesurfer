function [surf] = read_patch_asc(fname)
% function [ surf] = read_asc(fname)
%
% Reads asc-file with patch-information


%
% read_asc.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: fischl $
%    $Date: 2009/06/26 13:54:04 $
%    $Revision: 1.1 $
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

fp = fopen(fname, 'r');

% Dump first line
fgets(fp);

% Nr of vertices and faces
[a b] = fscanf(fp, '%d', 2);
nv = a(1) ;
nf = a(2) ;

fclose(fp);

% Read vertices and its indices
[ind vno x y z] = textread(fname, '%d vno=%d %f %f %f', nv, 'headerlines', 2, 'whitespace', '\t\b\n');

% Read faces and its indices
[find v1 v2 v3] = textread(fname, '%d %d %d %d', nf, 'headerlines', 2*nv+2, 'whitespace', '\t\b\n');

surf.faces = [v1 v2 v3];
surf.vertices = vno ;
surf.x = x;
surf.y = y;
surf.z = z;
