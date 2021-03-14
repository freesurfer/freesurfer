function [hsum, hsum2, dof] = fmri_accrandeff(h, hsum, hsum2, dof,...
                               rescale,roi)
%
% [hsum hsum2 dof] = fmri_accrandeff(h, hsum, hsum2, dof,rescale,roi)
%
% Accumulate values for random effects model.
%
%


%
% fmri_accrandeff.m
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


if(nargin < 4 | nargin > 6)
  msg = 'USAGE: [hsum hsum2 dof] = fmri_accrandeff(h, hsum, hsum2, dof,rescale,roi)'
  qoe(msg);error(msg);
end

if(nargin < 6)
  roi = [];
end

if(~isempty(roi)) 
  h = h(:,roi);
end

if(nargin >= 5)
  if(length(rescale)==1) 
    if(rescale ~= 1) h = h*rescale; end
  else
    hmin = min(reshape1d(h));
    hmax = max(reshape1d(h));
    hscale = (hmax-hmin)/(rescale(2) - rescale(1));
    h = (h-hmin)./hscale + rescale(1);
  end
end

%%%% Per-Voxel Processing %%%%%
if(isempty(roi)) 

  if(isempty(hsum)) 
    hsum  = zeros(size(h));
    hsum2 = zeros(size(h));
    dof   = 0;
  end

  hsum = hsum + h;
  hsum2 = hsum2 + h.*h;
  dof   = dof + 1;

%%%% Combine across voxels %%%
else             

  if(isempty(hsum)) 
    n = size(h,1);
    hsum  = zeros(n,1);
    hsum2 = zeros(n,1);
    dof   = 0;
  end

  hsum  = hsum + sum(h,2);
  hsum2 = hsum2 + sum(h.*h,2);
  dof   = dof + length(roi);

end


return;






return;
