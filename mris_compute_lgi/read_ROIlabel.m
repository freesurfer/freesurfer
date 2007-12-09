function [ pialROI ] = read_ROIlabel (fname)


l = [];

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

pialROI=l(:,1);
pialROI=pialROI+1;