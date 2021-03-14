function [cvm, nSupThresh] = fast_fmri2tcvm(f,RmMean,SegThresh)
% [cvm, nSupThresh] = fast_fmri2tcvm(f,RmMean,SegThresh)


%
% fast_fmri2tcvm.m
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


if(nargin < 1 | nargin > 3)
  msg = 'USAGE: [cvm, nSupThresh] = fast_fmri2tcvm(f,RmMean,SegThresh)';
  qoe(msg); error(msg);
end

if(nargin < 2)
  RmMean = 0;
end

if(nargin < 3)
  SegThresh = -1;
end

[ntrs nv] = size(f);

GlobalMean = mean(reshape1d(f));

if(SegThresh > 0 | RmMean)
  fmean = mean(f);
end

if(RmMean)
  f = f - repmat(fmean,[ntrs 1]);
end


if(SegThresh > 0)
  indSupThresh = find(fmean > SegThresh*GlobalMean);
else
  indSupThresh = 1:nv;
end

nSupThresh = length(indSupThresh);
if(nSupThresh == 0)
     fprintf('WARNING: no voxels above threshold\n');
  cvm = zeros(ntrs);
  return;
end

cvm = f(:,indSupThresh) * f(:,indSupThresh)'; %'

cvm = cvm/nv;


return;
