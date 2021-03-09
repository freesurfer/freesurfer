function nstats = fast_fgnoise(fmri,isepi,tpexc)
% nstats = fast_fgnoise(fmri,<isepi>,<tpexc>)
%
% isepi = 1 for epi (default) or 0 (eg, for spiral)
%
% nstats.z - zstat mristruct
% nstats.pz - pvals based on zstat (mristruct)
% nstats.fgmask - foreground/ghost mask mristruct 
% nstats.bgmask - background mask mristruct
% nstats.bgpzthresh - zthresh used to dtermine background
% nstats.nbgmask - nvox in bg
% nstats.bgmn - measured mean of bg
% nstats.bgvar - measured var of bg
% nstats.bgvarexp - expected bg var based on mean
% nstats.fgvarexp - expected fg var based on bgvar
% nstats.fgvar - measured fg var
% nstats.fgmask0 - init fgmask mristruct (> 2*global mean)
%
%


%
% fast_fgnoise.m
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

pzthresh = .01;
isfloored = 1;
isquant = 1;
ndil = 1;

nstats = [];
if(nargin < 1 | nargin > 3)
  fprintf('nstats = fast_fgnoise(fmri,<isepi>,<tpexc>)\n');
  return;
end

if(~exist('isepi','var')) isepi = []; end
if(isempty(isepi)) isepi = 1; end

if(~exist('tpexc','var')) tpexc = []; end
if(~isempty(tpexc))
  if(min(tpexc) < 1 | max(tpexc) > fmri.nframes)
    fprintf('ERROR: tpexc out of bounds\n');
    return;
  end
  tmp = ones(fmri.nframes,1);
  tmp(tpexc) = 0;
  tpinc = find(tmp);
  fmri.vol = fmri.vol(:,:,:,tpinc);
  fmri.nframes = size(fmri.vol,4);
end
  
f1 = fmri.vol(:,:,:,1:2:end);
fmn1  = mean(f1,4);
fstd1 = std(f1,[],4);
fgmn1 = mean(fmn1(:));
indz = find(fmn1==0);
fstd1(indz) = 1;
fmn1(indz) = 1000;

f2 = fmri.vol(:,:,:,2:2:end);
fmn2  = mean(f2,4);
fstd2 = std(f2,[],4);
fstd2(indz) = 1;
fmn2(indz) = 1000;

a = .1379; % emperical
z = (fmn1./fstd1 - 1/sqrt(4/pi-1))/a;
pz = 2*mri_zcdf(abs(z))-1;

fgmask0 = fmn1 > 2*fgmn1;
indfgmask0 = find(fgmask0);

if(isepi)
  % Ghost of the fg mask
  fgghost = fast_ghostmask(fgmask0); 
  % Create a mask of everything not wanted in background
  fgmask = fgghost | fgmask0; % Foreground and ghost
else
  fgmask = fgmask0;
end

if(ndil > 0)
  fprintf('Dilating FG by %d\n',ndil);
  fgmask = fast_dilate(fgmask,ndil);
end

% First pass at Background mask
bgmask0 = ~fgmask;

% assure that volume edges are exlcuded
bgmask0(1,:,:)   = 0;
bgmask0(end,:,:) = 0;
bgmask0(:,1,:)   = 0;
bgmask0(:,end,:) = 0;
bgmask0(:,:,1,:) = 0;
bgmask0(:,:,end,:) = 0;

indbgmask0 = find(bgmask0);

bgmask = pz < pzthresh & bgmask0;
indbgmask = find(bgmask);
nbgmask = length(indbgmask);
if(nbgmask == 0)
  fprintf('ERROR: could not segment background\n');
  return;
end

% Computes stats based on f2
bgmn  = mean(fmn2(indbgmask));
if(isfloored) bgmn  = bgmn + 0.5; end
bgvar = mean(fstd2(indbgmask).^2);
if(isquant) bgvar = bgvar - 1/12; end
bgvarexp = (bgmn*sqrt(4/pi-1)).^2;
fgvarexp =  bgvar/(2-pi/2);

% Foreground variance based on f2 in fgmask0
fgvar = mean(fstd2(indfgmask0).^2);

fprintf('nbgmask = %d, bgmn=%g, bgvar=%g, bgvarexp=%g, fgvarexp=%g,  fgvar=%g\n',...
	nbgmask,bgmn,bgvar,bgvarexp,fgvarexp,fgvar);


nstats.z = fmri;
nstats.z.vol = z;

nstats.pz = nstats.z;
nstats.pz.vol = pz;

nstats.fgmask0 = nstats.z;
nstats.fgmask0.vol = fgmask0;

nstats.fgmask = nstats.z;
nstats.fgmask.vol = fgmask;

nstats.bgmask = nstats.z;
nstats.bgmask.vol = bgmask;
nstats.indbgmask = indbgmask;

nstats.isepi = isepi;
nstats.isfloored = isfloored;
nstats.isquant = isquant;
nstats.ndil = ndil;

nstats.bgpzthresh = pzthresh;
nstats.nbgmask = nbgmask;
nstats.bgmn = bgmn;
nstats.bgvar = bgvar;
nstats.bgvarexp = bgvarexp;
nstats.fgvar = fgvar;
nstats.fgvarexp = fgvarexp;

nstats.fmn1 = nstats.z;
nstats.fmn1 = fmn1;
nstats.fstd1 = nstats.z;
nstats.fstd1 = fstd1;

nstats.fmn2 = nstats.z;
nstats.fmn2 = fmn2;
nstats.fstd2 = nstats.z;
nstats.fstd2 = fstd2;

return;

