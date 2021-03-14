function rt = fast_fxcfg(DoWhat,thing)
% rt = fast_fxcfg(DoWhat,thing)
%
% DoWhat - string instruction
% thing - either an fla or fx line. An fx line is only needed
%         for parseline.
%
% The fla has three fields that may or may not need to be set
% prior to calling fast_fxcfg:
%   nthfx 
%   nthrun
%   sesspath
%
% The return value depends upon the DoWhat
%
% DoWhats requiring neither nthfx nor nthrun
%
% loadsesscfg - loads session config (must set fla.sesspath, 
%   see also fast_sesscfg_struct). Returns the new fla.
% nfixedreg   - total number of fixed-effect regressors
% checkermid  - checks that all the ERMs are ok (ok = 1)
% parseline   - gets params from an fx line (thing = line);
%   see fast_fxcfg_FXMODEL.m where FXMODEL is name of the model.
%
% DoWhats requiring nthfx only
%
% getfxcfg    - gets nthfx fxcfg 
% createline  - creates an fx line 
% iserm       - returns 1 if nthfx is an ERM
% nparams     - number of params associated with nthfx
% nregressors - number of regressors associated with nthfx
%
% DoWhats requiring only the nthfx where the nthfx must be an ERM
%
% autopsd     - auto psd of nthfx
% irfmatrix   - impulse response matrix 
% irftaxis    - time associated with each row of irfmatrix 
% erfmatrix   - event response matrix 
% erftaxis    - time associated with each row of erfmatrix
% 
% DoWhats requiring nthrun only
%
% getevsch - gets nthrun evsch
% getntp   - gets number of timepoints in nthrun
% getitpx  - gets time point exclude indices for nthrun
%
% DoWhats requiring nthfx and (possibly) nthrun
%
% matrix - returns design matrix for nthfx (nthrun needed for rfx)
%
%
% See also: fast_flacfg_struct, fast_flacfg_load, fast_flacfg_write,
% fast_fxcfg_struct, fast_sesscfg_struct, fast_fla_desmat,
%
% See also FX Models:
% fast_fxcfg_poly, fast_fxcfg_fir, fast_fxcfg_extreg,
% fast_fxcfg_gamma, fast_fxcfg_spmhrf, fast_fxcfg_fourier
%
%
%


%
% fast_fxcfg.m
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

% Things to do:
%   nregressorstot - total number of regressors
%   nrandomreg     - total number of random fx regressors
%   indfxreg - indices of nthfx (nthrun needed for random fx) 
%   pmfmatrix - partial model fit for nthfx (nthrun needed for random fx)
%   main effect matrix for nthfx
%   omnibus matrix for nthfx
%   beta-to-fir matrix

% To add new FX Models, write the fast_fxcfg_NewModel.m file. Make
% sure that the model name in fx line is NewModel. Then add a case
% to switch(model). Note: if the switch model is replaced with
% eval, then no change to this file is necessary.

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
 case 'getitpx' % Time Point Exclude Indicies
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_itpx(flacfg);
 case 'getrunweight'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_runweight(flacfg);
 case 'loadsesscfg'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = load_sesscfg(flacfg);
 case 'nfixedreg'
  % number of fixed-effect regressors
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  nfx = length(flacfg.fxlist);
  svnthfx = flacfg.nthfx;
  nfixedreg = 0;
  for nthfx = 1:nfx
    if(strcmp(flacfg.fxlist(nthfx).fx.fxtype,'fixed'))
      flacfg.nthfx = nthfx;
      nr = fast_fxcfg('nregressors',flacfg);
      nfixedreg = nfixedreg + nr;
    end
  end
  flacfg.nthfx = svnthfx;
  rt = nfixedreg;
 
 case 'checkermid'
  % returns 1 if ok, 0 or [] otherwise
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = check_ermid(flacfg);
 
 case {'iserm','nparams','nregressors','autopsd',...
       'irfmatrix','erfmatrix','matrix','parseline','createline',...
      'irftaxis','erftaxis'}
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

  % Should do this with an eval %
  switch(model)
   case 'fir'
    rt = fast_fxcfg_fir(DoWhat,thing);
   case 'gamma'
    rt = fast_fxcfg_gamma(DoWhat,thing);
   case 'spmhrf'
    rt = fast_fxcfg_spmhrf(DoWhat,thing);
   case 'polynomial'
    rt = fast_fxcfg_poly(DoWhat,thing);
   case 'fourier'
    rt = fast_fxcfg_fourier(DoWhat,thing);
   case 'extreg'
    rt = fast_fxcfg_extreg(DoWhat,thing);
   otherwise
    fprintf('ERROR: model %s unrecognized\n',model);
    return;
  end
 
  if(strcmp(DoWhat,'matrix'))
    itpx = fast_fxcfg('getitpx',flacfg);
    if(~isempty(itpx) & flacfg.usetpexclude)
      rt(itpx,:) = 0;
    end
    if(~isempty(flacfg.nskip) & flacfg.nskip > 0)
      itpx = [1:flacfg.nskip];
      rt(itpx,:) = 0;
    end
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
  fprintf('ERROR: get_ntp: nthrun=%d, out of range\n');
  return;
end

ntp = flacfg.sesscfg.ntp(flacfg.nthrun);

return;

%------------------------------------------------------------%
function evsch = get_evsch(flacfg)
evsch = [];

% Check here to make sure model is erm?

if(isempty(flacfg.nthrun))
  fprintf('ERROR: get_evsch: nthrun is empty\n');
  return;
end
if(isempty(flacfg.sesscfg))
  fprintf('ERROR: get_evsch: sesscfg is empty\n');
  return;
end
if(isempty(flacfg.sesscfg.evschlist))
  fprintf('ERROR: get_evsch: flacfg.sesscfg.evschlist is empty\n');
  return;
end
nruns = length(flacfg.sesscfg.evschlist);
if(flacfg.nthrun < 1 | flacfg.nthrun > nruns)
  fprintf('ERROR: get_evsch: nthrun=%d, out of range\n',flacfg.nthrun);
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
  fprintf('ERROR: get_runweight: nthrun=%d, out of range\n',flacfg.nthrun);
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
function flacfg = load_sesscfg(flacfg)

if(isempty(flacfg.funcstem))
  fprintf('ERROR: flacfg.funcstem is empty\n');
  flacfg = []; return;
end
if(isempty(flacfg.sesspath))
  fprintf('ERROR: flacfg.sesspath is empty\n');
  flacfg = []; return;
end
d = dir(flacfg.sesspath);
if(isempty(d))
  fprintf('ERROR: flacfg.sesspath %s does not exist\n',flacfg.sesspath);
  flacfg = []; return;
end
afsd = sprintf('%s/%s',flacfg.sesspath,flacfg.fsd);
d = dir(afsd);
if(isempty(d))
  fprintf('ERROR: %s does not exist\n',afsd);
  flacfg = []; return;
end

runlist = fast_runlist(afsd,flacfg.runlistfile);
if(isempty(runlist))
  fprintf('ERROR: with runs or runlistfile in %s\n',afsd);
  flacfg = []; return;
end

sesscfg = fast_sesscfg_struct;
sesscfg.runlist = runlist;
sesscfg.fstemlist = '';

nfx = length(flacfg.fxlist);

% If one of the fxcfg needs run weight, read file and check 
% that each run is represented.


for nthrun = 1:length(runlist)
  rd = sprintf('%s/%s',afsd,runlist(nthrun,:));
  d = dir(rd);
  if(isempty(d))
    fprintf('ERROR: %s does not exist\n',rd);
    flacfg = []; return;
  end

  fstem = sprintf('%s/%s',rd,flacfg.funcstem);
  [nslices nrows ncols ntp] = fmri_bvoldim(fstem);
  if(isempty(nslices))
    fprintf('ERROR: reading %s\n',fstem);
    flacfg = []; return;
  end
  sesscfg.fstemlist = strvcat(sesscfg.fstemlist,fstem);
  sesscfg.ntp(nthrun) = ntp;
  sesscfg.volsize = [nrows ncols nslices];
  % Could check spat dims here %

  if(~isempty(flacfg.evschfname))
    evschfile = sprintf('%s/%s',rd,flacfg.evschfname);
    evsch = fmri_ldpar(evschfile);
    if(isempty(evsch))
      fprintf('ERROR: reading %s\n',evschfile);
      flacfg = []; return;
    end
    sesscfg.evschlist(nthrun).evsch = evsch;
  end

  % Go through each effect, determine if any are extreg %
  for nthfx = 1:nfx
    fxcfg = flacfg.fxlist(nthfx).fx;
    if(strcmp(fxcfg.model,'extreg'))
      extregstem = sprintf('%s/%s',rd,deblank(fxcfg.sparams(1,:)));
      extreg = squeeze(fast_ldbslice(extregstem))';
      %keyboard
      if(isempty(extreg))
	fprintf('ERROR: reading %s\n',extregstem);
	flacfg = []; return;
      end
      if(size(extreg,1) ~= ntp)
	fprintf('ERROR: dimension mismatch %s\n',extregstem);
	flacfg = []; return;
      end
      if(size(extreg,1) < fxcfg.params(1))
	fprintf('ERROR: nextreg=%d, exceeds size of extreg\n',fxcfg.params(1));
	flacfg = []; return;	
      end
      flacfg.fxlist(nthfx).fx.npmlist(nthrun).M  = extreg;
    end
  end
end
  
flacfg.sesscfg = sesscfg;

return;
%---------------------------------------------------------%
function itpx = get_itpx(flacfg)

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
  fprintf('ERROR: get_itpx: nthrun=%d, out of range\n',flacfg.nthrun);
  return;
end

evsch = flacfg.sesscfg.evschlist(flacfg.nthrun).evsch;
itpx = find(evsch(:,2) == -1);

return;
