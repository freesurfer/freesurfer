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
