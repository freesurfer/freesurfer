% svdfunctcvm.m

% allrunlist, runlist, tcvmstem

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
