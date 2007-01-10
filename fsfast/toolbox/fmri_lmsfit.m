function h = fmri_lmsfit(Ch,ss)
%
% h = fmri_lmsfit(Ch,ss)
%
% Ch: nTotEst x nTotEst
% ss: nRows x nCols x nTotEst
% h:  nRows x nCols x nTotEst
%
%


%
% fmri_lmsfit.m
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

if(nargin ~= 2)
  msg = 'USAGE: fmri_lmsfit(Ch,ss)';
  qoe(msg);
  error(msg);
end

nTotEst = size(Ch,1);

szs = size(ss);
ns = length(szs);

% Determine nRows, nCols, and nTotEst %
if(ns == 2)
  nRows = size(ss,1);
  nCols = 1;
  nTotEst = size(ss,2);
else
  nRows = size(ss,1);
  nCols = size(ss,2);
  nTotEst = size(ss,3);
end

% Compute number of voxels %
nV = nRows*nCols; 

% Check dimensions of Ch and ss %%
if(size(Ch,1) ~= nTotEst | size(Ch,2) ~= nTotEst)
  msg = 'Ch and ss dimensions are inconsistent';
  qoe(msg);
  error(msg);
end

% Reshape ss for matrix mult %
ss = reshape(ss, [nV nTotEst])'; %'

% Heres where the action is %
h = Ch * ss;

% Reshape the output to the proper dim %
h = reshape(h', [nRows nCols nTotEst]); %'


return;

