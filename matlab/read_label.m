function [l] = read_label(sname, lname)
% l = read_label(<sname>, lname)
%
% reads the label file 'lname' from the subject 'sname' 
% in the subject's label directory into the vector l
% l will be nvertices-by-5, where each column means:
% (1) vertex number, (2-4) xyz at each vertex, (5) stat
%
% IMPORTANT: the vertex number is 0-based.
% 
% $Id: read_label.m,v 1.4 2005/01/23 17:19:23 greve Exp $

l = [];

if(nargin ~= 2)
  fprintf('l = read_label(<sname>, lname)\n');
  return;
end

if(~isempty(sname))
  sdir = getenv('SUBJECTS_DIR') ;
  fname = sprintf('%s/%s/label/%s.label', sdir, sname, lname) ;
else
  fname = lname;
end

% open it as an ascii file
fid = fopen(fname, 'r') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

fgets(fid) ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

line = fgets(fid) ;
nv = sscanf(line, '%d') ;
l = fscanf(fid, '%d %f %f %f %f\n') ;
l = reshape(l, 5, nv) ;
l = l' ;

fclose(fid) ;

