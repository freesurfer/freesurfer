function [S, v, f] = read_moviebyu(fname)
% function [v, f] = read_moviebyu(fname)
% Reads an ascii-file (movie.byu format). Returns a vertices and
% faces matrix, ready to input in patches.


%
% read_moviebyu.m
%
% Original Author: Bruce Fischl
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
