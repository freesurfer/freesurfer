function [S, v, f] = read_asc(fname)
% function [S, v, f] = read_asc(fname)
%
% Reads asc-file with patch-information


%
% read_asc.m
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

fp = fopen(fname, 'r');

% Dump first line
fgets(fp);

% Nr of vertices and faces
S = zeros(1, 3);
S(1) = 1;
[S(2:3)] = fscanf(fp, '%d', 2);

% Read vertices and its indices
v = fscanf(fp, '%f', S(2)*4);
v = reshape(v, [4 S(2)])';
v(:, 1) = v(:, 1)+1;

% Read faces and its indices
f = fscanf(fp, '%d', S(3)*5);
f = reshape(f, [5 S(3)])';
f = f+1;

fclose(fp);
