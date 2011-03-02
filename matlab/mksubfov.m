% mksubfov.m
%
% The purpose of this matlab script is to create the
% mni305.cor.subfov1.mgz and mni305.cor.subfov2.mgz volumes.  These
% are the same data as found in mni305.cor.mgz, but the field-of-view
% is much smaller and only covers the brain. The purpose of this is to
% reduce the amount of space needed to store the data. This is
% especially important when doing group functional analysis where
% there might be many subjects combined. In one subfov (subfov1), the
% voxel size is 1mm isotropic, 151 x 151 x 186. In the other (yep,
% subfov2), its 2mm isotropic,  76 x 76 x 93. These volumes are in
% register with mni305.cor.mgz in that the share the same RAS
% space, ie, you can run:
%   tkregister2 --targ mni305.cor.mgz --mov mni305.cor.subfov1.mgz \
%        --regheader --reg /tmp/reg
% And the volumes will be in register.
%
% After these files are created by this program, run the following:
%   mri_convert mni305.cor.subfov1.mgz mni305.cor.subfov1.mgz --odt uchar
%   mri_convert mni305.cor.subfov2.mgz mni305.cor.subfov2.mgz --odt uchar
% This just reduces the size by a factor of 4.
%
% You should be in $DEV/distribution/averages when running this script


%
% mksubfov.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.5 $
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

% Load in the pixel values. These are 1mm, isotropic
cor = MRIread('mni305.cor.mgz');
if(isempty(cor)) 
  fprintf('ERROR: cannot find mni305.cor.mgz\n');
  fprintf('Make sure you are in $DEV/distribution/averages\n');
  return; 
end

% Binary mask of the brain.
m = MRIread('mni305.mask.cor.mgz');
if(isempty(m)) 
  fprintf('ERROR: cannot find mni305.mask.cor.mgz\n');
  fprintf('Make sure you are in $DEV/distribution/averages\n');
  return; 
end

% Create a bounding box from the mask, expanded by 2 voxels
ind = find(m.vol);
[r c s] = ind2sub(m.volsize,ind); % These are all 1-based
rmin = min(r)-4; % do 4 here 'cause it's close to the edge
rmax = max(r)+2;
fovr = rmax - rmin;
cmin = min(c)-2;
cmax = max(c)+2;
fovc = cmax - cmin;
smin = min(s)-2;
smax = max(s)+2;
fovs = smax - smin;
fovrc = max(fovr,fovc);

% 1mm isotropic SubFOV ---------------------------------------
crsP0 = [cmin rmin smin 1]';  % crs at first voxel of FOV
P0 = m.vox2ras1*crsP0; % use vox2ras1 since crs are 1-based
D = diag(m.volres); % Diagonal matrix of voxel sizes
vox2ras = [m.Mdc*D P0(1:3); 0 0 0 1]; % but this is 0-based
cor2 = cor;
cor2.vol = cor.vol(rmin:rmin+fovrc,cmin:cmin+fovrc,smin:smax);
cor2.vox2ras0 = vox2ras;
MRIwrite(cor2,'mni305.cor.subfov1.mgz');

% 2mm isotropic SubFOV ---------------------------------------
% Subsample, don't interpolate
crsP0 = [cmin rmin smin 1]';
P0 = m.vox2ras1*crsP0; % use vox2ras1 since crs are 1-based
D = 2*diag(m.volres); % use 2-times volres
vox2ras = [m.Mdc*D P0(1:3); 0 0 0 1];
cor2 = cor;
cor2.vol = cor.vol(rmin:2:rmin+fovrc,cmin:2:cmin+fovrc,smin:2:smax);%skip2
cor2.vox2ras0 = vox2ras;
MRIwrite(cor2,'mni305.cor.subfov2.mgz');

