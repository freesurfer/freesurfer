function rt = fast_fxcfg(DoWhat,thing)
% rt = fast_fxcfg(DoWhat,thing)

% Things to do:
%   loadsesscfg - load run weights
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
 case 'loadsesscfg'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = load_sesscfg(flacfg);
 
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
    model = sscanf(line,'%*s %*s %*s %s ',1);
    if(isempty(model))
      fprintf('ERROR: line is not formatted correctly (no model)\n');
      fprintf('%s\n',line);
      return;
    end
  else
    if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
    fxcfg = fast_fxcfg('getfxcfg',flacfg);
    if(isempty(fxcfg)) return; end
    model = fxcfg.model;
  end

  switch(model)
   case 'fir'
    rt = fast_fxcfg_fir(DoWhat,thing);
   case 'gamma'
    rt = fast_fxcfg_gamma(DoWhat,thing);
   case 'polynomial'
    rt = fast_fxcfg_poly(DoWhat,thing);
   otherwise
    fprintf('ERROR: model %s unrecognized\n',model);
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
evlabellist = '';
for nthfx = 1:nfx
  flacfg.nthfx = nthfx;
  if(fast_fxcfg('iserm',flacfg))
    fxcfg = fast_fxcfg('getfxcfg',flacfg);
    if(isempty(fxcfg)) return; end
    evidlist = [evidlist fxcfg.params(1)];
    evlabellist = strvcat(evlabellist,fxcfg.label);
  end
end

evidlistunique = unique(evidlist);
if(length(evidlist) ~= length(evidlistunique))
  fprintf('ERROR: replications found in event id list\n');
  fprintf('%d ',evidlist);
  fprintf('\n');
  ok = 0;
  return;
end

evlabellistunique = unique(evlabellist,'rows');
if(size(evlabellist,1) ~= size(evlabellistunique,1))
  fprintf('ERROR: replications found in event label list\n');
  fprintf('%s ',evlabellist');
  fprintf('\n');
  ok = 0;
  return;
end

ok = 1;
return;
%---------------------------------------------------------------%
function sesscfg = load_sesscfg(flacfg)

sesscfg = [];

if(isempty(flacfg.funcstem))
  fprintf('ERROR: flacfg.funcstem is empty\n');
  return;
end
if(isempty(flacfg.sesspath))
  fprintf('ERROR: flacfg.sesspath is empty\n');
  return;
end
d = dir(flacfg.sesspath);
if(isempty(d))
  fprintf('ERROR: flacfg.sesspath %s does not exist\n',flacfg.sesspath);
  return;
end
afsd = sprintf('%s/%s',flacfg.sesspath,flacfg.fsd);
d = dir(afsd);
if(isempty(d))
  fprintf('ERROR: %s does not exist\n',afsd);
  return;
end

if(isempty(flacfg.runlistfile))
  runlist = fast_runlist(afsd);
  if(isempty(runlist))
    fprintf('ERROR: cannot find any runs in %s\n',afsd);
    return;
  end
else
  rlf = sprintf('%s/%s',afsd,flacfg.runlistfile);
  runlist = fast_runlist(rlf);
  if(isempty(runlist))
    fprintf('ERROR: reading %s\n',rlf);
    return;
  end
end

sesscfg = fast_sesscfg_struct;
sesscfg.runlist = runlist;
sesscfg.fstemlist = '';

% Need to read run weight file 
% and check that each run is represented.

for nthrun = 1:length(runlist)
  rd = sprintf('%s/%s',afsd,runlist(nthrun,:));
  d = dir(rd);
  if(isempty(d))
    fprintf('ERROR: %s does not exist\n',rd);
    sesscfg = [];
    return;
  end

  fstem = sprintf('%s/%s',rd,flacfg.funcstem);
  [nslices nrows ncols ntp] = fmri_bvoldim(fstem);
  if(isempty(nslices))
    fprintf('ERROR: reading %s\n',fstem);
    sesscfg = [];
    return;
  end
  sesscfg.fstemlist = strvcat(sesscfg.fstemlist,fstem);
  sesscfg.ntp(nthrun) = ntp;
  % Could check spat dims here %

  if(~isempty(flacfg.evschfname))
    evschfile = sprintf('%s/%s',rd,flacfg.evschfname);
    evsch = fmri_ldpar(evschfile);
    if(isempty(evsch))
      fprintf('ERROR: reading %s\n',evschfile);
      sesscfg = [];
    end
    sesscfg.evschlist(nthrun).evsch = evsch;
  end

end
  
return;