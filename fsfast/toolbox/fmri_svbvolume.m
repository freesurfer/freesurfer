function fmri_ldbvolume(vol,stem,voldim,ext)
% fmri_ldbvolume(vol,stem,voldim,ext)
%
% Saves a volume in bfile format where vol is a 4D structure
% of dimension Nslices X Nrows X Ncols X Ndepth.
%


%
% fmri_svbvolume.m
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

if(nargin ~= 2 & nargin ~= 3 & nargin ~= 4)    
  fprintf(2,'USAGE: fmri_svbvolume(vol,stem,<voldim>,<ext>)\n');
  qoe; error; return;
end

if(nargin > 2)  
  if(~isempty(voldim))  
    vol = reshape(vol, voldim); 
  else
    voldim = size(vol);
  end
else
  voldim = size(vol);
end

if(nargin ~= 4)  ext = 'bfloat'; end

nslices = voldim(1);
for slice = 0:nslices-1
  % fprintf('Saving Slice %3d\n',slice);
  fname = sprintf('%s_%03d.%s',stem,slice,ext);
  z = shiftdim(vol(slice+1,:,:,:),1);
  fmri_svbfile(z, fname);
end

return;
