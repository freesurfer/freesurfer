function [fid] = write_curv(curv, fname)

%
% writes a curvature vector into an ascii file
%


fid = fopen(fname, 'w') ;
nvertices = size(curv,1) ;
fprintf(fid, '%d\n', nvertices) ;
for i=1:nvertices
		fprintf(fid, '0 0.0 0.0 0.0 %f\n', curv(i,1)) ;
end
fclose(fid) ;