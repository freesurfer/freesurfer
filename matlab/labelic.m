function err = labelic(subject,apslice)
% err = labelic(subject,apslice)
%
% Simple program that attempts to label the internal capsule based
% upon the automatic segmentation. It is very crude, to say the
% least, but it might serve as a good starting point. 
%
% Uses SUBJECTS_DIR/subject/mri/aseg.mgz as input.
%
% Writes SUBJECTS_DIR/subject/mri/asegic.mgz as output using
% IC Labels are as found in FreeSurferColors.txt
%  155  Left-IntCapsule-Ant    
%  156  Right-IntCapsule-Ant   
%  157  Left-IntCapsule-Pos    
%  158  Right-IntCapsule-Pos   
%
% Anterior/Posterior is defined by apslice, ie, voxels anterior to
% coronal apslice are labeled anterior, etc. If unspecified, apslice
% defaults to 128.
%
% Also writes SUBJECTS_DIR/subject/mri/icmask.mgz as output. This
% is a binary mask of what the full IC label and may be helpful 
% when editing the asegic.mgz.
% 
% 
% tkmedit subject orig.mgz -aux icmask.mgh -segmentation ./asegic.mgz
%
%


%
% labelic.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.3 $
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



ndil = 6;

err = 1;
if(nargin < 1 | nargin > 2)
  fprintf('err = labelic(subject,apslice)\n');
  return;
end

if(~exist('apslice','var')) apslice = []; end
if(isempty(apslice)) apslice = 128; end

fprintf('Starting labelic\n');
fprintf('apslice = %d\n',apslice);

SUBJECTS_DIR = deblank(getenv('SUBJECTS_DIR'));

fspec = sprintf('%s/%s/mri/aseg.mgz',SUBJECTS_DIR,subject);
aseg = MRIread(fspec);
if(isempty(aseg)) return; end

% First, Create a binary mask of these structures
% Putamen
% Caudate
% Pallidum
% Lateral-Ventricle
% Thalamus-Proper
% VentralDC 28 60
m0 = (aseg.vol == 10 | aseg.vol == 11 | ...
      aseg.vol == 12 | aseg.vol == 13 | ...
      aseg.vol ==  4 | aseg.vol == 43 | ...
      aseg.vol == 49 | aseg.vol == 50 | ...
      aseg.vol == 51 | aseg.vol == 52 | ...
      aseg.vol == 28 | aseg.vol == 60);

% Dilate the mask. Needs to dilate enough to bridge the gap across
% the IC.
fprintf('Dilating %d \n',ndil);
tic;
md0 = fast_dilate(m0,ndil);
fprintf('  done %g\n',toc);

% Now erode by the same amount. Wont erode where mask bridged the
% IC gap.
fprintf('Eroding %d \n',ndil);
tic;
md = fast_dilate(md0,ndil,1);
fprintf('  done %g\n',toc);

wmleft  = (aseg.vol == 2);    % Left WM
icleft  = md & ~m0 & wmleft;  % Left IC
indicleft = find(icleft);
[r c s] = ind2sub(size(icleft),indicleft);
iant = find(s >  apslice);
ipos = find(s <= apslice);
indicleftant = indicleft(iant);
indicleftpos = indicleft(ipos);
aseg.vol(indicleftant) = 155;
aseg.vol(indicleftpos) = 157;

wmright = (aseg.vol == 41);    % Right WM
icright  = md & ~m0 & wmright; % Right IC
indicright = find(icright);
[r c s] = ind2sub(size(icright),indicright);
iant = find(s >  apslice);
ipos = find(s <= apslice);
indicrightant = indicright(iant);
indicrightpos = indicright(ipos);
aseg.vol(indicrightant) = 156;
aseg.vol(indicrightpos) = 158;

fspec = sprintf('%s/%s/mri/asegic.mgz',SUBJECTS_DIR,subject);
MRIwrite(aseg,fspec);

% Create the IC Mask
icmask = aseg;
icmask.vol = icleft | icright;
fspec = sprintf('%s/%s/mri/icmask.mgz',SUBJECTS_DIR,subject);
MRIwrite(icmask,fspec);

err = 0;

return;
