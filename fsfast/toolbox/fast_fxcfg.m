function rt = fast_fxcfg(DoWhat,thing)
% rt = fast_fxcfg(DoWhat,thing)

% Things to do:
%   getfstempath
%   getevschfile
%   loadsesscfg
%   

rt = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('rt = fast_fxcfg(DoWhat,<thing>)\n');
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
 case 'getfxcfg'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_fxcfg(flacfg);
 case 'getevsch'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_evsch(flacfg);
 case 'getntp'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_ntp(flacfg);
 case 'getrunweight'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_runweight(flacfg);
 case 'checkermid'
  % returns 1 if ok, 0 or [] otherwise
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = check_ermid(flacfg);
 
 case {'iserm','nparams','nregressors','autopsd',...
       'amatrix','matrix','parseline','createline'}
  if(strcmp(DoWhat,'parseline'))
    if(isempty(line))
      fprintf('ERROR: line needed with parseline\n');
      return;
    end
  else
    if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  end
  fxcfg = fast_fxcfg('getfxcfg',flacfg);
  if(isempty(fxcfg)) return; end

  switch(fxcfg.model)
   case 'fir'
    rt = fast_fxcfg_fir(DoWhat,thing);
   otherwise
    fprintf('ERROR: model %s unrecognized\n',fxcfg.model);
    return;
  end
 
 otherwise
  fprintf('ERROR: DoWhat=%s, unregconized\n',DoWhat);
end

return;
%------------------------------------------------------------%
%------------>>>>>>>>>>>>>>><<<<<<<<<<<<<<<------------------%
%------------------------------------------------------------%

%------------------------------------------------------------%
function pr_fla_needed(DoWhat)
fprintf('ERROR: flacfg needed with %s\n',DoWhat);
return;

%---------------------------------------------------------------%
function fxcfg = get_fxcfg(flacfg)
% Get the current (nth) fxcfg
fxcfg = [];
if(isempty(flacfg.nthfx))
  fprintf('ERROR: nthfx is empty\n');
  return;
end
nfx = length(flacfg.fxlist);
if(flacfg.nthfx < 1 | flacfg.nthfx > nfx)
  fprintf('ERROR: nthfx = %d, out of range\n',flacfg.nthfx);
  return;
end
fxcfg = flacfg.fxlist(flacfg.nthfx).fx;

return;

%------------------------------------------------------------%
function ntp = get_ntp(flacfg)
ntp = [];

if(isempty(flacfg.nthrun))
  fprintf('ERROR: nthrun is empty\n');
  return;
end
if(isempty(flacfg.sesscfg))
  fprintf('ERROR: sesscfg is empty\n');
  return;
end
nruns = length(flacfg.sesscfg.ntp);
if(flacfg.nthrun < 1 | flacfg.nthrun > nruns)
  fprintf('ERROR: nthrun=%d, out of range\n');
  return;
end

ntp = flacfg.sesscfg.ntp(flacfg.nthrun);

return;

%------------------------------------------------------------%
function evsch = get_evsch(flacfg)
evsch = [];

% Check here to make sure model is erm?

if(isempty(flacfg.nthrun))
  fprintf('ERROR: nthrun is empty\n');
  return;
end
if(isempty(flacfg.sesscfg))
  fprintf('ERROR: sesscfg is empty\n');
  return;
end
if(isempty(flacfg.sesscfg.evschlist))
  fprintf('ERROR: flacfg.sesscfg.evschlist is empty\n');
  return;
end
nruns = length(flacfg.sesscfg.evschlist);
if(flacfg.nthrun < 1 | flacfg.nthrun > nruns)
  fprintf('ERROR: nthrun=%d, out of range\n');
  return;
end

evsch = flacfg.sesscfg.evschlist(flacfg.nthrun).evsch;
ncols = size(evsch,2);
if(ncols ~= 2 & ncols ~= 3)
  fprintf('ERROR: evsch has %d cols, must be 2 or 3\n',ncols);
  evsch = [];
  return;
end

return;

%---------------------------------------------------------------%
function w = get_runweight(flacfg)
w = [];

if(isempty(flacfg.nthrun))
  fprintf('ERROR: nthrun is empty\n');
  return;
end
if(isempty(flacfg.sesscfg))
  fprintf('ERROR: sesscfg is empty\n');
  return;
end
if(isempty(flacfg.sesscfg.runweight))
  % Or try to load if flacfg.runweightfile is defined?
  w = 1;
  return;
end
nruns = length(flacfg.sesscfg.runweight);
if(flacfg.nthrun < 1 | flacfg.nthrun > nruns)
  fprintf('ERROR: nthrun=%d, out of range\n');
  return;
end

w = flacfg.sesscfg.runweight(flacfg.nthrun);

return;

%---------------------------------------------------------------%
function ok = check_ermid(flacfg)
% checks for replications in ERM EvIds

ok = [];

nfx = length(flacfg.fxlist);
evidlist = [];
for nthfx = 1:nfx
  flacfg.nthfx = nthfx;
  if(fast_fxcfg('iserm',flacfg))
    fxcfg = fast_fxcfg('getfxcfg',flacfg);
    if(isempty(fxcfg)) return; end
  end
  evidlist = [evidlist fxcfg.params(1)];
end

evidlistunique = unique(evidlist);
if(length(evidlist) ~= length(evidlistunique))
  fprintf('ERROR: replications found int event id list\n');
  fprintf('%d ',evidlist);
  fprintf('\n');
  ok = 0;
  return;
end

ok = 1;
return;