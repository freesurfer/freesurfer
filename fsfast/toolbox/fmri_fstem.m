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
% $Id: fmri_fstem.m,v 1.1 2003/03/04 20:47:39 greve Exp $
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
