function [S, v, f] = read_asc(fname)
% function [S, v, f] = read_asc(fname)
%
% Reads asc-file with patch-information

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
