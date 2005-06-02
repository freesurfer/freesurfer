function [code, name, rgbv] = read_fscolorlut(fname)
% [code name rgb] = read_fscolorlut(fname)
%
% Reads a freesurfer color lookup table. By default
% reads $FREESURFER_HOME/FreeSurferColorLUT.txt.
%
% $Id: read_fscolorlut.m,v 1.1 2005/06/02 22:42:34 greve Exp $

code = [];

if(nargin == 0)
  FREESURFER_HOME = getenv('FREESURFER_HOME');
  fname = sprintf('%s/FreeSurferColorLUT.txt',FREESURFER_HOME);
end

fp = fopen(fname,'r');
if(fp == -1) 
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

name = '';
nthitem = 1;
while(1)
  
  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fp);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end
    
  c = sscanf(tline,'%d',1);
  n = sscanf(tline,'%*d %s',1);
  r = sscanf(tline,'%*d %*s %d',1);
  g = sscanf(tline,'%*d %*s %*d %d',1);
  b = sscanf(tline,'%*d %*s %*d %*d %d',1);
  v = sscanf(tline,'%*d %*s %*d %*d %*d %d',1);
  code(nthitem,1) = c;
  name = strvcat(name,n);
  rgbv(nthitem,:) = [r g b v];

  nthitem = nthitem + 1;
end

fclose(fp);

return;







