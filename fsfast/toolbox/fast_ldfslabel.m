function l = fast_ldfslabel(lname, sname)
% l = fast_ldfslabel(lname, sname)
%
% Reads the label file lname. If sname is included,
% then the label path will be SUBJECTS_DIR/sname/label/lname.
% 
% Returns the matrix l which will be npoints X 5, where
% the five columns are: (1) vertex number, (2-4) xyz,
% (5) value associated with label point.
%
% Based on read_label.m by Bruce F.

l = [];

if(nargin ~=1 & nargin ~=2)
  fprintf('USAGE: l = fast_ldfslabel(lname, sname)\n');
  return;
end

if(nargin ==2)
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
line = fgets(fid) ;
nv = sscanf(line, '%d') ;
l = fscanf(fid, '%d %f %f %f %f\n') ;
fclose(fid) ;

l = reshape(l, 5, nv) ;
l = l' ; %'


return;
