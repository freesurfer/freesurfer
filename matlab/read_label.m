function [l] = read_label(sname, lname)

% l = read_label(sname, lname)

%
% reads the label file 'lname' from the subject 'sname' 
% in the subject's label directory into the vector l
%

sdir = getenv('SUBJECTS_DIR') ;
fname = sprintf('%s/%s/label/%s.label', sdir, sname, lname) ;

% open it as an ascii file
fid = fopen(fname, 'r') ;
fgets(fid) ;
line = fgets(fid) ;
nv = sscanf(line, '%d') ;
l = fscanf(fid, '%d %f %f %f %f\n') ;
l = reshape(l, 5, nv) ;
l = l' ;
fclose(fid) ;

