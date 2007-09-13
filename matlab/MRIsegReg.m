function err = MRIsegReg(subject)
% err = MRIsegReg(subject)
% $Id: MRIsegReg.m,v 1.1 2007/09/13 23:26:12 greve Exp $

err = 1;
if(nargin ~= 1)
  fprintf('err = MRIsegReg(subject)\n');
  return;
end

SUBJECTS_DIR = getenv('SUBJECTS_DIR');
sd = sprintf('%s/%s',SUBJECTS_DIR,subject);
if(exist(sd) ~= 7)
  fprintf('ERROR: cannot find subject %s in %s\n',subject, ...
	  SUBJECTS_DIR);
  return;
end

fprintf('Reading aparc+aseg\n');
fspec = sprintf('%s/mri/aparc+aseg.mgz',sd);
apas = MRIread(fspec);
if(isempty(apas)) return; end
aparcaseg = apas.vol;

fprintf('Reading brainmask\n');
fspec = sprintf('%s/mri/brainmask.mgz',sd);
bm = MRIread(fspec);
if(isempty(bm)) return; end
brainmask = bm.vol;

% Get a mask that excludes the edge of the brain
% as well as B0 areas:
%  1. X016 - parahippocampal
%  2. X032 - frontal pole
indmask0 = find(brainmask > 0.5 & ...
		aparcaseg ~= 1016 & aparcaseg ~= 2016 & ...
		aparcaseg ~= 1032 & aparcaseg ~= 2032);
nmask0 = length(indmask0);
mask = zeros(size(aparcaseg));
mask(indmask0) = 1;
nerode = 5;
fprintf('Eroding mask by %d\n',nerode);tic
mask = fast_dilate(mask,5,1); % 3D erode by 5 vox = 5 mm
fprintf(' ... done %g\n',toc);
indmask = find(mask);
nmask = length(indmask);
printf('nmask0 = %d, nmask = %d',nmask0,nmask);

m = apas;
m.vol = mask;
fprintf('Writing regmask\n');
fspec = sprintf('%s/mri/regmask.mgz',sd);
MRIwrite(m,fspec);

% Get WM by eroding WM by 3 voxels (3mm)
indwm = find((aparcaseg == 2 | aparcaseg == 41) & mask);
nwm = length(indwm);
wm = zeros(size(aparcaseg));
wm(indwm) = 1;
fprintf('Eroding WM\n');tic
wm = fast_dilate(wm,2,1);
fprintf(' ... done %g\n',toc);
indwm2 = find(wm);
nwm2 = length(indwm2);
printf('nwm = %d, nwm2 = %d',nwm,nwm2);

% Now get Cortex
indctx = find( ((aparcaseg >= 1000 & aparcaseg <= 1034 ) | ...
		(aparcaseg >= 2000 & aparcaseg <= 2034 ) ) & ...
		mask);
nctx = length(indctx);
ctx = zeros(size(aparcaseg));
ctx(indctx) = 1;

regseg = apas;
regseg.vol = 41*wm + 3*ctx; % 41=green

fspec = sprintf('%s/mri/regseg.mgz',sd);
MRIwrite(regseg,fspec);

err = 0;

return;
