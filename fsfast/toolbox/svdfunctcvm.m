% svdfunctcvm.m

% allrunlist, runlist, tcvmstem


%
% svdfunctcvm.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
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

tcvmfile = sprintf('%s.bfloat',tcvmstem);
tcvmall = fmri_ldbfile(tcvmfile);
if(isempty(tcvmall))
  fprintf('ERROR: reading %s\n',tcvmfile);
  error; qoe;
end

nruns = length(runlist);
nallruns = length(allrunlist);

ntalltot = size(tcvmall,1);

if(mod(ntalltot,nallruns) ~= 0)
  fprintf('ERROR: dimension mismatch: ntalltot=%d, nallruns=%d\n',...
          nalltot,nallruns);
  error; qoe;
end  

nt = ntalltot/nallruns;
nttot = nt * nruns;

for run = runlist
  if(isempty(find(allrunlist==run)))
    fprintf('ERROR: run %d not found in allrunlist\n',run);
    error; qoe;
  end  
end

% Remove the runs not in the run list from tcvmall %
rmlist = [];
nthrun = 1;
for run = allrunlist
  if(isempty(find(runlist==run)))
    rmlist = [rmlist nthrun];
  end
  nthrun = nthrun + 1;
end

tcvm = fast_tcvm_rmrun(rmlist,tcvmall,nallruns);

%-------- Construct the design matrix -----------%
