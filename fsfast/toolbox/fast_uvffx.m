function [p,t,gffx,gvarffx] = fast_uvffx(gam,gamvar,dof)
% [p,t,gffx,gvarffx] = fast_uvffx(gam,gamvar,dof)
% 
% Univariate fixed-effects model. The inputs are gam, each with
% variance gamvar computed with the given dof. The means from each
% sample are averaged to get the ffx mean, and the vars from each
% sample are averaged and then divided by the number of samples.
%
% gam cat be a matrix nframes-by-nvox or mristruct.
% gamvar must be the same size as gam.
% dof must be either a scalar or a vector of nframes.
%
% p and t are tests of mean(gam) = 0. If gam is an mristruct, then
% p and t will be mristructs.
%
%


%
% fast_uvffx.m
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

p = [];
t = [];
gffx = [];
gvarffx = [];

if(nargin ~= 3)
  fprintf('[p,t,gffx,gvarffx] = fast_uvffx(gam,gamvar,dof)\n');
  return;
end

if(isstruct(gam))
  gammat = fast_vol2mat(gam.vol);
  sz = size(gam.vol);
else
  gammat = gam;
end

if(isstruct(gamvar))
  gamvarmat = fast_vol2mat(gamvar.vol);
else
  gamvarmat = gamvar;
end

nframes = size(gammat,1);
nvox    = size(gammat,2);

if(size(gammat,1) ~= size(gamvarmat,1))
  fprintf('ERROR: gam and gamvar differ in number of frames\n');
  return;
end

if(size(gammat,2) ~= size(gamvarmat,2))
  fprintf('ERROR: gam and gamvar differ in number of voxels\n');
  return;
end

if(length(dof) == 1) dof = repmat(dof,[nframes 1]); end
if(length(dof) ~= nframes)
  fprintf('ERROR: dof list does not equal 1 or nframes\n');
  return;
end

% FFx mean is just mean of inputs
gffxmat = mean(gammat);

% Variance = sum(var*dof)/(nframes*sum(dof));
ffxdof = sum(dof);
rdofmat = repmat(dof,[1 nvox])/(nframes*ffxdof);
gvarffxmat = sum(gamvarmat .* rdofmat);

% t = mean/std(mean)
indnz = find(gvarffxmat ~= 0);
tmat = zeros(1,nvox);
tmat(indnz) = gffxmat(indnz) ./ sqrt(gvarffxmat(indnz));

% pvalues
pmat = tTest(ffxdof, tmat, 300);

if(isstruct(gam))
  gffx = gam;
  gffx.vol = fast_mat2vol(gffxmat,sz);
  gvarffx = gam;
  gvarffx.vol = fast_mat2vol(gvarffxmat,sz);
  t = gam;
  t.vol = fast_mat2vol(tmat,sz);
  p = gam;
  p.vol = fast_mat2vol(pmat,sz);
else
  gffx = gffxmat;
  gvarffx = gvarffxmat;
  t = tmat;
  p = pmat;
end  


return;





















