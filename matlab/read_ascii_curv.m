function [curv] = read_curv(fname)

%
% reads an ascii curvature into a vector
%


fid = fopen(fname, 'r') ;
nvertices = fscanf(fid, '%d', 1);
all = fscanf(fid, '%d %f %f %f %f\n', [5, nvertices]) ;
curv = all(5, :)' ;
