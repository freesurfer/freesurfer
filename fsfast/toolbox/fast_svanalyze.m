function err = fast_svanalyze(vol,volid,voxsize,M)
% err = fast_svanalyze(vol,volid,voxsize,M)
%
% Uses spm functions, so you have to have spm installed and
% in your path.
%
% volid is the name of the volume WITHOUT extension
% voxsize is the size of each dimension of vol. Set to []
%   for all ones
% M is the registration matrix. Set to [] for spm default.
%
% vol CANNOT have multiple time points because spm cannot
%   write them.


%
% fast_svanalyze.m
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

err = 1;

if(nargin ~= 4)
  fprintf('USAGE: err = fast_svanalyze(vol,volid,voxsize,M)\n');
  return;
end

if(isempty(voxsize)) voxsize = ones(length(size(vol)),1); end
voxsize = reshape1d(voxsize);

if(length(size(vol)) ~= length(voxsize))
  fprintf('ERROR: the voxsize must have the same number of elements\n');
  fprintf('       as the volume has dimensions\n');
  return;
end

anafile = sprintf('%s.img',volid);
DIM = size(vol);
VOX = voxsize;
SCALE = 1;
TYPE = 4; % float=16, short=4
OFFSET = 0;
ORIGIN = [0 0 0];
DESCRIP = '';

if(isempty(M)) 
  matfile = sprintf('%s.mat',volid);
  delete(matfile);
end

nwritten = spm_hwrite(anafile,DIM,VOX,SCALE,TYPE,OFFSET,ORIGIN,DESCRIP);
if(nwritten ~= 348)
  fprintf('ERROR: writing header for %s\n',anafile);
  return;
end

volspec = spm_vol(anafile);
if(isempty(volspec))
  fprintf('ERROR: could not read %s voxel data\n',anafile);
  return;
end

fprintf('INFO: Writing %s\n',anafile);
spm_write_vol(volspec,vol);

if(~isempty(M)) 
  matfile = sprintf('%s.mat',volid);
  save(matfile,'M');
end

err = 0;

return;
