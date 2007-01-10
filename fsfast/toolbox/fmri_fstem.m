function fstem = fmri_fstem(varargin)
%
% fstem = fmri_fstem(BFileName)
%
% Returns the stem of a file.  The stem is
% the file name/path up to but not including
% the last underscore character '_'. BFileName
% can be vertically concatenated list or a 
% comma-separated list of file names. fstem
% is a vertically concatenated list of the
% stems.
%
%
%


%
% fmri_fstem.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

if(nargin == 0) 
  msg = 'USAGE: fmri_fstem(BFileName)';
  qoe(msg); error(msg);
end

if( length(varargin) == 1)
  BFileList = varargin{1};
  nRuns = size(BFileList,1);
else
  nRuns = length(varargin);
  BFileList = '';
  for r = 1:nRuns,
    BFileList = strvcat(BFileList,varargin{r});
  end
end

fstem = [];

for r = 1:nRuns,
  BFileName = deblank(BFileList(r,:));
  k = max(findstr(BFileName,'_'));

  Base = BFileName(1:k-1);
  fstem = strvcat(fstem,Base);
end
