function runlist = fast_runlistfile(runlistfile)
% runlist = fast_runlistfile(runlistfile)


%
% fast_runlistfile.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
%    $Revision: 1.4 $
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

runlist = [];

if(nargin ~= 1)
  msg = 'runlist = fast_runlistfile(runlistfile)';
  qoe(msg);error(msg);
end

fid = fopen(runlistfile);
if(fid == -1)
  msg = sprintf('Could not open %s',runlistfile);
  qoe(msg);error(msg);
end

runid = deblank(fscanf(fid,'%s',1));
while( ~isempty(runid) )
  runlist = strvcat(runlist,runid);
  runid = fscanf(fid,'%s',1);
end

if(isempty(runlist))
  msg = sprintf('ERROR: no runs found in %s\n',runlistfile);
  qoe(msg);error(msg);
end

return;
