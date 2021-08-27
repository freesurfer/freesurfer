function [paropt,trmin,tsearched,nsearched, optstats] = fmri_optseq(oss)


%
% fmri_optseq.m
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

paropt  = [];
trmin   = 10^10*ones(oss.Nsessions,1);
nopt    = zeros(oss.Nsessions,1);

oss.Nc   = length(oss.Npercond) + 1;
oss.Nh   = ceil(oss.TimeWindow/oss.TER);
oss.Nnnc = oss.Nc - 1 ;
oss.Trun = oss.TR*oss.Ntp;
oss.DOF  = oss.Nruns * oss.Ntp - oss.Nh * oss.Nnnc;

if(oss.DOF < 0)
  msg = sprintf('Not enough degrees of freedom: DOF = %d',oss.DOF);
  qoe(msg);error(msg);
end

if(oss.pforder > 0)
  Xpf = [];
  for n = 1:oss.Nruns
    Xpf = [Xpf fast_polytrendmtx(n,oss.Ntp,oss.Nruns,oss.pforder)];
  end
else
  Xpf = [];
end

GammaFit = 0;
if(~isempty(oss.GammaParams)) GammaFit = 1; end

ok = 1;
n  = 1;
t0 = clock;
trsum = 0;
tr2sum = 0;
effsum = 0;
eff2sum = 0;
trmax = 0;
effmin = 10^10;
effmax = 0;
dorand = 0;
if(GammaFit)
  Nbeta_par = oss.Nnnc;
  gfDelta = oss.GammaParams(1);
  gfTau   = oss.GammaParams(2);
else
  Nbeta_par = oss.Nnnc*oss.Nh;
end

while(ok)

  %%% Synthesize a paradigm file %%%
  par = fmri_synthpar3(oss.Npercond,oss.Tpercond,oss.Nruns,...
                       oss.Trun,oss.Tres,oss.TPreScan);

  Xfir = fmri_par2scm(par,oss.Nc,oss.Ntp,oss.TER,oss.Nh, ...
		      oss.TPreStim);
  if(GammaFit)
    Xpar = fmri_scm2gcm(Xfir,oss.Nc-1,oss.TER,oss.TPreStim,gfDelta,gfTau);
  else
    Xpar = Xfir;
  end

  X = [Xpar Xpf];
  Ch = fmri_hcovar(X);
  %Ch = inv(X'*X); %'
  if(~isempty(Xpf)) Ch = Ch(1:Nbeta_par,1:Nbeta_par); end
  tr = trace(Ch);
  if(oss.FindWorst) tr = 1/tr; end % for testing
  trsum  = trsum  + tr;
  tr2sum = tr2sum + tr*tr;
  eff = 1/tr;
  effsum  = effsum  + eff;
  eff2sum = eff2sum + eff*eff;

  if(n==1)
     paropt = zeros(size(par,1),2,oss.Nruns,oss.Nsessions);
  end

  if(eff > effmax) effmax = eff; end
  if(eff < effmin) effmin = eff; end

  if(tr > trmax) trmax = tr; end

  if(eff > 1/trmin(oss.Nsessions) & eff < oss.MaxEffLimit)
    % MaxEffLimit is there for testing purposes; normally
    % it is infinite
    npar = size(par,1);

    if(npar > size(paropt,1))
      % This is a hack for the case where the pars have
      % a different number of entries.
      tmp = zeros(npar,2,oss.Nruns,oss.Nsessions);
      tmp(1:size(paropt,1),:,:,:) = paropt;
      paropt = tmp;
    end

    paropt(:,:,:,oss.Nsessions) = zeros(size(paropt(:,:,:,oss.Nsessions)));
    paropt(1:npar,:,:,oss.Nsessions) = par;
    trmin(oss.Nsessions) = tr;
    nopt(oss.Nsessions)  = n;

    [tmp indx] = sort(trmin);
    trmin = tmp;
    paropt = paropt(:,:,:,indx);
    nopt   = nopt(indx);

  end

  n  = n + 1;
  dt = etime(clock,t0);

  if(oss.Nsearch > 0)
    if(n > oss.Nsearch) break; end
  else
    if(dt > oss.Tsearch) break; end
  end

end

tsearched = dt;
nsearched = n-1;

travg = trsum/n;
trstd = sqrt(tr2sum/n - travg^2);

effavg = effsum/n;
effstd = sqrt(eff2sum/n - effavg^2);

optstats = [travg trstd trmax effavg effstd effmin effmax];

% This is a hack to assure that the presentation times 
% are monotonically increasing.  Non-mono can be caused
% by padding optpar when the individual runs/sessions
% have different number of entries.
for s = 1:oss.Nsessions
  for r = 1:oss.Nruns
    partmp = paropt(:,1,r,s);
    npar = size(partmp,1);
    n0 = min(find(diff(partmp)<0))+1;
    if(~isempty(n0))
      t0 = partmp(n0-1,1);
      t  = t0 + oss.TR*[1:npar-n0+1]'; %'
      paropt(n0:npar,1,r,s) = t;
    end
  end
end

return;
