function rt = fast_fxcfg_fir(DoWhat,thing)
% rt = fast_fxcfg_fir(DoWhat,thing)
%
% DoWhat can be:
%  iserm  - returns 1 if erm, 0 otherwise
%  nparams 
%  nregressors
%  matrix
%  amatrix (for erm only)
%  parseline
%  createline
%  (autopsd not available for fir)
%
% thing - line to be parsed or flacfg (see fast_flacfg_struct).
%
% FIR Parameters:
%  1. EventId
%  2. BoxCarWidth
%  3. PSDMin
%  4. dPSD
%  5. PSDMax
%

rt = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('rt = fast_fxcfg_fir(DoWhat,<thing>)\n');
  return;
end

if(exist('thing') ~= 1) 
  thing = []; 
  flacfg = [];
  line = [];
else
  if(isfield(thing,'flaname')) 
    flacfg = thing;
    line = [];
  else 
    line = thing;
    flacfg = [];
  end
end

DoWhat = lower(DoWhat);
switch(DoWhat)
 
 case 'iserm'
  rt = 1;
 
 case 'nparams'
  rt = 5;
 
 case 'nregressors'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_nregressors(flacfg);
 
 case 'autopsd'
  fprintf('ERROR: fir: cannot autopsd\n');
  rt = []; 
 
 case 'amatrix'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  nreg = fast_fxcfg_fir('nregressors',flacfg);
  if(isempty(nreg)) return; end
  rt = eye(nreg);
 
 case 'matrix'

  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  X = get_matrix(flacfg);
  if(isempty(X)) return; end
  rt = X;
 
 case 'parseline'
  if(isempty(line))
    fprintf('ERROR: fir: line needed with parseline\n');
    return;
  end
  rt = parseline(line);
  
 otherwise
  fprintf('ERROR: fast_erm_fir: getwhat = %s, unrecognized\n',getwhat);

end

return;

%------------------------------------------------------------%
function pr_fla_needed(DoWhat)
fprintf('ERROR: flacfg needed with %s\n',DoWhat);
return;

%------------------------------------------------------------%
function fxcfg = parseline(line)
% Read and check input line
% InputLine: Effect Fixed Label FIR EvId BCW PSDMin dPSD PSDMax
fxcfg = [];

nparams = fast_fxcfg_fir('nparams');
[tmp nitems] = sscanf(line,'%s',inf);
if(nitems ~= nparams + 4)
  fprintf('ERROR: fir: line has wrong number of items (%d)\n',nitems);
  fprintf('%s\n',line);
  return;
end

fxcfg = fast_fxcfg_struct;

[fxcfg.fxtype n]  = sscanf(line,'%*s %s ',1);
[fxcfg.label  n]  = sscanf(line,'%*s %*s %s',1);
[fxcfg.model  n]  = sscanf(line,'%*s %*s %*s %s',1);

fxcfg.fxtype = lower(fxcfg.fxtype);
fxcfg.model = lower(fxcfg.model);

if(~strcmp(fxcfg.fxtype,'fixed'))
  fprintf('ERROR: fir: must be fixed effect type\n');
  fprintf('%s\n',line);
  fxcfg = [];
  return;
end
if(~strcmp(fxcfg.model,'fir'))
  fprintf('ERROR: fir: not fir model\n');
  fprintf('%s\n',line);
  fxcfg = [];
  return;
end

fmt0 = '%*s %*s %*s %*s';
for p = 1:nparams
  fmt = sprintf('%s %%f',fmt0);
  fxcfg.params(p) = sscanf(line,fmt,1);
  fmt0 = sprintf('%s %%*f',fmt0);
end

evid = fxcfg.params(1);

if(evid < 1)
  fprintf('ERROR: fir: evid=%d < 1\n',evid);
  fxcfg = [];
  return;
end

psdwin = get_psdwin(fxcfg.params);
if(isempty(psdwin)) 
  fxcfg = [];
  return;
end


return;

%------------------------------------------------------------%
function nthfx = get_nthfx(flacfg)

nthfx = [];

if(isempty(flacfg.nthfx))
  fprintf('ERROR: nthfx is empty\n');
  return;
end

nfx = length(flacfg.fxlist);
if(flacfg.nthfx < 1 | flacfg.nthfx > nfx)
  fprintf('ERROR: nthfx = %d, out of range\n',flacfg.nthfx);
  return;
end

nthfx = flacfg.nthfx;

return;
%------------------------------------------------------------%
function psdwin = get_psdwin(params)
psdwin = [];

bcw     = params(2);
psdmin  = params(3);
dpsd    = params(4);
psdmax  = params(5);

if(rem(bcw,dpsd) ~= 0)
  fprintf('ERROR: fir: bcw=%g not int mult of dpsd=%g\n',bcw,dpsd);
  return; 
end
if(rem(bcw,psdmin) ~= 0)
  fprintf('ERROR: fir: psdmin=%g not int mult of dpsd=%g\n',psdmin,dpsd);
  return; 
end
if(rem(psdmax,dpsd) ~= 0)
  fprintf('ERROR: fir: psdmax=%g not int mult of dpsd=%g\n',psdmax,dpsd);
  return; 
end

psdmax = psdmax + bcw - dpsd;
psdwin = [psdmin dpsd psdmax];

return;

%------------------------------------------------------------%
function fxcfg = get_fxcfg(flacfg)
fxcfg = [];

nthfx = get_nthfx(flacfg);
if(isempty(nthfx)) return; end
fxcfg = flacfg.fxlist(nthfx).fx;

return;

%------------------------------------------------------------%
function nr = get_nregressors(flacfg)

nr = [];

fxcfg = get_fxcfg(flacfg);
if(isempty(fxcfg)) return; end

psdwin = get_psdwin(fxcfg.params);
if(isempty(psdwin)) return; end

psdmin  = psdwin(1);
dpsd    = psdwin(2);
psdmax  = psdwin(3);

nr = (psdmax-psdmin)/dpsd;

return;

%------------------------------------------------------------%
function nthrun = get_nthrun(flacfg)
nthrun = [];

if(isempty(flacfg.nthrun))
  fprintf('ERROR: fir: nthrun is empty\n');
  return;
end
if(isempty(flacfg.sesscfg))
  fprintf('ERROR: fir: sesscfg is empty\n');
  return;
end
nruns = length(flacfg.sesscfg.ntp);
if(flacfg.nthrun < 1 | flacfg.nthrun > nruns)
  fprintf('ERROR: fir: nthrun=%d, out of range\n');
  return;
end
  
nthrun = flacfg.nthrun;

return;

%------------------------------------------------------------%
function ntp = get_ntp(flacfg)
ntp = [];

nthrun = get_nthrun(flacfg);
if(isempty(nthrun)) return; end
ntp = flacfg.sesscfg.ntp(nthrun);

return;

%------------------------------------------------------------%
function  X = get_matrix(flacfg)
X = [];

nr = get_nregressors(flacfg);
if(isempty(nr)) return; end
ntp = get_ntp(flacfg);
if(isempty(ntp)) return; end
fxcfg = get_fxcfg(flacfg);
if(isempty(fxcfg)) return; end

[tPres, EvIdPres, wPres] = get_evsch(flacfg);
if(isempty(tPres)) return; end

evid   = fxcfg.params(1);
bcw    = fxcfg.params(2);
psdmin = fxcfg.params(3);
dpsd   = fxcfg.params(4);
psdmax = fxcfg.params(5);
psd = [psdmin dpsd psdmax];

indEvId = find(EvIdPres == evid);
tPresEvId = tPres(indEvId);
if(~isempty(wPres)) wPresEvId = wPres(indEvId); end

X = fast_sched2Xfir(tPresEvId,ntp,flacfg.TR,psd,bcw,flacfg.tDelay,wPres); 

return;

%------------------------------------------------------------%
function [tpres, evidpres, wpres] = get_evsch(flacfg)
tpres = [];
evidpres = [];
wpres = [];

nthrun = get_nthrun(flacfg);
if(isempty(nthrun)) return; end
% Should check for erm too

nruns_evsch = length(flacfg.sesscfg.evsch);
if(nruns_evsch < nthrun)
  fprintf('ERROR: fir: flacfg.sesscfg.evsch not enough runs\n');
  return;
end

if(isempty(flacfg.sesscfg.evsch(nthrun).tpres))
  fprintf('ERROR: fir: flacfg.sesscfg.evsch(%d).tpres is empty\n',nthrun);
  return;
end

if(isempty(flacfg.sesscfg.evsch(nthrun).evidpres))
  fprintf('ERROR: fir: flacfg.sesscfg.evsch(%d).evidpres is empty\n',nthrun);
  return;
end

tpres = flacfg.sesscfg.evsch(nthrun).tpres;
evidpres = flacfg.sesscfg.evsch(nthrun).evidpres;
wpres = flacfg.sesscfg.evsch(nthrun).wpres;

if(length(tpres) ~= length(evidpres))
  fprintf('ERROR: fir: tpres and evid have a different number of reps\n');
  tpres = []; evidpres = []; wpres = [];
  return;
end

if(~isempty(wpres))
  if(length(tpres) ~= length(wpres))
    fprintf('ERROR: fir: tpres and wpres have a different number of reps\n');
    tpres = []; evidpres = []; wpres = [];
    return;
  end
end

return;





