function [annots] = read_annotation(fname)

%
% [annots] = read_annotation(fname)
% reads a binary annotation file into a vector
%



% open it as a big-endian file
fid = fopen(fname, 'rb', 'b') ;
if (fid < 0)
	 str = sprintf('could not open annotation file %s.', fname) ;
	 error(str) ;
end
vnum = fread(fid, 1, 'int32') ;
tmp = fread(fid, vnum*2, 'int') ; 
annots = tmp(2:2:vnum*2) ;
	   	
fclose(fid) ;
