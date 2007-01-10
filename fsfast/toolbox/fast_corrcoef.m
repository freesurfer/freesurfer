function [cc, yvar] = fast_corrcoef(y,lag,DOF)
% [cc, yvar] = fast_corrcoef(y,<lag>,<DOF>)
%
% y is nframes by ncols
% lag defaults to 1 if not present or null
% DOF defaults to nframes if not present or null
%
% cc is 1 by ncols set of normalized correlation coefficients
% yvar is 1 by ncols
%


%
% fast_corrcoef.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
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

cc = [];
yvar = [];

if(nargin < 1 | nargin > 3)
  fprintf('[cc, yvar] = fast_corrcoef(y,<lag>,<DOF>)\n');
  return; 
end

[nf nv] = size(y);

if(nargin < 2)   lag = 1; end
if(isempty(lag)) lag = 1; end
if(nargin ~= 3)  DOF = nf; end
if(isempty(DOF)) DOF = nf; end

if(DOF-lag <= 0)
  fprintf('ERROR: DOF-lag <= 0\n');
  return;
end

nn1 = 1:nf-lag;
nn2 = lag+1:nf;

yvar = sum(y.^2,1)/DOF;
ycvar = sum(y(nn1,:).*y(nn2,:),1)/(DOF-lag);
cc = ycvar./yvar;

return;
