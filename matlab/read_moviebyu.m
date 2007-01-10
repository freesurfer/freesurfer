function [S, v, f] = read_moviebyu(fname)
% function [v, f] = read_moviebyu(fname)
% Reads an ascii-file (movie.byu format). Returns a vertices and
% faces matrix, ready to input in patches.


%
% read_moviebyu.m
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

fp = fopen(fname, 'r');

% Read in first two lines
S = fscanf(fp, '%d', 4);
fscanf(fp, '%d', 2);

v = zeros(6, ceil(S(2)/2));
% Read vertices
v(:, 1:floor(S(2)/2)) = fscanf(fp, '%f ', [6, floor(S(2)/2)]);

if S(2)/2 ~= floor(S(2)/2)
   v(1:3, ceil(S(2)/2)) = fscanf(fp, '%f', 3);
   v = reshape(v, [3 S(2)+1])';
   v = v(1:S(2), :);
else
	v = reshape(v, [3 S(2)])';
end

% Read faces
f = fscanf(fp, '%d', [4 S(3)]);
f(4, :) = -f(4, :);
f = f+1;
f = f';

S = S(1:3);

fclose(fp);
