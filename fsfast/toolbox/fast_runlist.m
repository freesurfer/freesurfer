function runlist = fast_runlist(dirname,runlistfile)
% runlist = fast_runlist(dirname,<runlistfile>)
% runlistfile, if present, should be relative to dirname


%
% fast_runlist.m
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'runlist = fast_runlist(dirname,<runlistfile>)';
  qoe(msg);error(msg);
end

d = dir(dirname);
if(isempty(d))
  msg = sprintf('No runs found in %s\n',dirname);
  fprintf('%s',msg);
  return;
  %qoe(msg); error(msg);
end

if(~exist('runlistfile')) runlistfile = ''; end
if(~isempty(runlistfile))
  rlf = sprintf('%s/%s',dirname,runlistfile);
  runlist = fast_runlistfile(rlf);
  % Should check that they exist too
  return;
end

for n = 1:length(d);
  dname = d(n).name;
  if(length(dname)==3)
    if(~isempty(str2num(dname))) 
      runlist = strvcat(runlist,dname);
    end
  end
end

if(isempty(runlist))
  msg = sprintf('ERROR: no runs found in %s\n',dirname);
  qoe(msg);error(msg);
end

return;
