function [havg, hstd, hdrdat] = fmri_untangle(havgstd, hdrdat, SubtCond0)
% 
% [havg, hstd] = fmri_untangle(havgstd, hdrdat, <SubtCond0>)
%
% havg: rows cols cond delay
%


%
% fmri_untangle.m
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


if(nargin ~= 2 & nargin ~= 3)
  msg = 'Usage: [havg, hstd] = fmri_untangle(havgstd, hdrdat, <SubtCond0>)';
  qoe(msg);error(msg);
end

if(size(havgstd,2) == 1 & length(size(havgstd)) == 2)
  hdrdat.nRows = 1;
  hdrdat.nCols = 1;
else
  [hdrdat.Nrows hdrdat.Ncols Nch2] = size(havgstd);
end

% row, cols, tangled-as ---> tangled-ascond, rows, cols 
hsa2 = permute(havgstd, [3 1 2]);
% tangled-ascond, rows, cols --> stat statid condid rows cols
hsa3 = reshape(hsa2, [hdrdat.Nh 2 hdrdat.Nc hdrdat.Nrows hdrdat.Ncols]); 

havg = squeeze(hsa3(:,1,:,:,:)); % extract average
hstd = squeeze(hsa3(:,2,:,:,:)); % extract stdev
if(nargin == 3) % SubtCond0
  havg = havg(:,[2:hdrdat.Nc],:,:) - repmat(havg(:,1,:,:),[1 hdrdat.Nnnc 1 1]);
  % hstd = sqrt(hstd.^2 + repmat(hstd(:,1,:,:).^2,[1 hdrdat.Nc 1 1]));
end

% average condid rows cols --> rows cols condid average
havg = permute(havg, [3 4 2 1]);
% stdev condid rows cols --> rows cols condid stdev
hstd = permute(hstd, [3 4 2 1]);

return
