function tcvm = fast_tcvm_rmrun(rmlist,tcvmall,nruns)
% tcvm = fast_tcvm_rmrun(rmlist,tcvmall,nruns)
%
% This will remove the runs listed in rmlist from the temporal
% covariance matrix tcvmall. It is assumed that tcvmall was
% created from nruns, each nt = size(tcvmall)/nruns, so that
% tcvmall constsists of nruns X nruns blocks, with each block 
% being nt X nt. If M is a member of rmlist, then this function
% will remove the Mth block row and column from tcvmall.


%
% fast_tcvm_rmrun.m
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

tcvm = [];

if(nargin ~= 3)
  fprintf('USAGE: tcvm = fast_tcvm_rmrun(rmlist,tcvmall,nruns)\n');
  return;
end

if(max(rmlist)>nruns)
  fprintf('ERROR: max(rmlist) = %d > nruns = %d\n',max(rmlist),nruns);
  return;
end

if(min(rmlist)<1)
  fprintf('ERROR: min(rmlist) = %d < 1\n',min(rmlist));
  return;
end

ntall = size(tcvmall,1);

if(mod(ntall,nruns) ~= 0)
  fprintf('ERROR: dimension mismatch: ntall=%d, nruns=%d\n',...
          ntall,nruns);
  return;
end  

if(isempty(rmlist))
  tcvm = tcvmall;
  return;
end

% number of time-points per run %
nt = ntall/nruns;

nrm = length(rmlist);
nkeep = nruns - nrm;

%ntkeep = nt*nkeep;
%tcvm = zeros(ntkeep);

indkeep = [];
for run = [1:nruns]
  if(isempty(find(rmlist==run)))
    nmin = (run-1)*nt + 1;
    nmax = nmin + nt - 1;
    indkeep = [indkeep [nmin:nmax]];
  end
end


tcvm = tcvmall(indkeep,indkeep);

return;
