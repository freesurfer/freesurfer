function ok = write_label(lindex, lxyz, lvals, labelfile, subjid,space)
% ok = write_label(lindex, lxzy, lvals, labelfile, <subjid>,<space name>)


%
% write_label.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%
% where space name can be "voxel", "tkreg", "scanner"

ok = 0;

if(nargin ~= 4 & nargin ~= 5 & nargin ~= 6)
  fprintf('ok = write_label(lindex, lxzy, lvals, labelfile, <subjid> <space>)\n');
  return;
end

if(exist('subjid') ~= 1) subjid = 'matlab'; end
if(exist('space') ~= 1) space = 'tkReg'; end

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

if(~isempty(lvals)) 
  if(npoints ~= length(lvals))
    fprintf('ERROR: length of lvals inconsistent\n');
    return;
  end
else
  lvals  = zeros(npoints,1); 
end

% open as an ascii file
fid = fopen(labelfile, 'w') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',labelfile);
  return;
end

fprintf(fid,'#!ascii label, from subject %s, vox2ras=%s \n',subjid,space);
fprintf(fid,'%d\n',npoints);

% Make sure they are npoints by 1 %
lindex = reshape(lindex,[npoints 1]);
lxyz   = reshape(lxyz,[npoints 3]);
lvals  = reshape(lvals,[npoints 1]);

l = [lindex lxyz lvals];
fprintf(fid,'%d %f %f %f %f\n',l') ;

fclose(fid) ;

ok = 1;

return;
