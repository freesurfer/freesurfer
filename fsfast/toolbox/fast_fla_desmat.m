function X = fast_fla_desmat(flacfg,runlist,jkrun)
% X = fast_fla_desmat(flacfg,<runlist>,<jkrun>)

X = [];

if(nargin < 1 | nargin > 3)
  fprintf('X = fast_fla_desmat(flacfg,<runlist,jkrun>)\n');
  return;
end

if(exist('runlist')~=1) runlist = []; end

if(isempty(flacfg.sesscfg.runlist))
  fprintf('ERROR: no runs in sesscfg run list\n');
  return;
end
nruns = length(flacfg.sesscfg.runlist);

if(exist('jkrun')~=1) jkrun = []; end
if(~isempty(jkrun))
  ind = find(flacfg.sesscfg.runlist == jkrun);
  if(isempty(ind))
    fprintf('ERROR: jkrun = %d, not in run list ');
    fprintf('%d ',flacfg.sesscfg.runlist);
    fprintf('\n');
    return;
  end
end

if(isempty(flacfg.fxlist))
  fprintf('ERROR: no fx defined in fla.\n');
  return;
end
nfx = length(flacfg.fxlist);

Xfe = [];
Xre = [];
for nthrun = 1:nruns
  if(~isempty(jkrun) & nthrun == jkrun) continue; end
  flacfg.nthrun = nthrun;
  Xferun = [];
  Xrerun = [];
  for nthfx = 1:nfx
    flacfg.nthfx = nthfx;
    fx = fast_fxcfg('getfxcfg',flacfg);
    if(isempty(fx)) return; end
    Xtmp = fast_fxcfg('matrix',flacfg);
    if(isempty(Xtmp)) return; end
    
    if(strcmp(fx.fxtype,'fixed'))  Xferun = [Xferun Xtmp];  end
    if(strcmp(fx.fxtype,'random')) Xrerun = [Xrerun Xtmp]; end
  end

  Xfe = [Xfe; Xferun]; % Should check size here
  Za = zeros(size(Xre,1),size(Xrerun,2));
  Zb = zeros(size(Xrerun,1),size(Xre,2));

  Xre = [Xre Za; Zb Xrerun];
  
end

X = [Xfe Xre];


return;