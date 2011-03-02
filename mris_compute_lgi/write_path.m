function ok = write_path (lxyz, lindex, pathfile)
% ok = write_label(lxyz, lindex, pathfile)

%
% modified from write_label.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:30 $
%    $Revision: 1.2 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

ok = 0;

if(isempty(lindex) & isempty(lxyz))
  fprintf('ERROR: both lindex and lxyz are empty.\n');
  return;
end

if(~isempty(lindex) & ~isempty(lxyz))
  npoints1 = length(lindex); 
  npoints2 = size(lxyz,1); 
  if(npoints1 ~= npoints2)
    fprintf('ERROR: lindex and lxyz have different lengths.\n');
    return;
  end
  npoints = npoints1;
elseif(~isempty(lindex))
  npoints = length(lindex); 
  lxyz = zeros(npoints,3); 
elseif(~isempty(lxyz))
  npoints = length(lxyz); 
  lindex = zeros(npoints,1); 
end

if(size(lxyz,2) ~= 3)
  fprintf('ERROR: lxyz does not have 3 columns\n');
  return;
end

% open as an ascii file
fid = fopen(pathfile, 'w') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',labelfile);
  return;
end

fprintf(fid,'# Path file\n');
fprintf(fid,'VERSION 2\n');
fprintf(fid,'BEGINPATH\n') ;
fprintf(fid,'NUMPOINTS %d\n',npoints);

% Make sure they are npoints by 1 %
lindex = reshape(lindex,[npoints 1]);
lxyz   = reshape(lxyz,[npoints 3]);

l = [lxyz lindex];
fprintf(fid,'%f %f %f %d \n',l');
fprintf(fid,'ENDPATH');

fclose(fid) ;

ok = 1;

return;
