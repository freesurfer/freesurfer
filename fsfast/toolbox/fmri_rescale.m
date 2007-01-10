function [ypost, oldmean] = fmri_rescale(ypre,newmean,TPExclude)
% [ypost oldmean] = fmri_rescale(ypre,newmean,TPExclude)
%
% Rescales slices on a run-by-run basis.  The scaling factor for each
% run is f = newmean/mean(run) where mean(run) is the mean computed
% over all voxels at all time-points (ie, the global mean for the run).
% This forces the global mean to be the same for all runs.
% The oldmean for each slice is returned.
% If TPExclude is passed as an argument, data points that correspond
% to TPExclude==1 are ignored.  They are replaced with newmean in
% ypost.  The dimensions of TPExclude should be nTP by nRuns.
%
% Comments or questions: analysis-bugs@nmr.mgh.harvard.edu
%
%


%
% fmri_rescale.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: [ypost, oldmean] = fmri_rescale(ypre,newmean,<TPExclude>)';
  qoe(msg);  error(msg);
end

%% Get dimensions of functional slices
[nRows nCols nTP nRuns] = size(ypre);

if(nargin == 3)
  if(size(TPExclude,1) ~= nTP)
    msg = 'ypre and TPExclude differ in the number of time points';
    qoe(msg);  error(msg);
  elseif(size(TPExclude,2) ~= nRuns)
    msg = 'ypre and TPExclude differ in the number of runs';
    qoe(msg);  error(msg);
  end
end

%% Set all of the outputs to be the new mean %%%
ypost = newmean*ones(size(ypre));

%% Go through each run %%%
for r = 1:nRuns,

  if(nargin == 3) % Use exclusions %
    iTPInclude = find(TPExclude(:,r)==0);
    oldmean(r) = mean(reshape1d(ypre(:,:,iTPInclude,r)));
    ypost(:,:,iTPInclude,r) = (newmean/oldmean(r))*ypre(:,:,iTPInclude,r);

  else % Dont use exclusions %
    oldmean(r) = mean(reshape1d(ypre(:,:,:,r)));
    ypost(:,:,:,r) = (newmean/oldmean(r))*ypre(:,:,:,r);
  end

end

return;
