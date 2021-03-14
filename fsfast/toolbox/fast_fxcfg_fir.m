function rt = fast_fxcfg_fir(DoWhat,thing)
% rt = fast_fxcfg_fir(DoWhat,thing)
%
% DoWhat can be:
%  iserm  - returns 1 if erm, 0 otherwise. thing not needed.
%  nparams - returns number of parameters in model (5). thing not needed.
%  nregressors - number of regressors in current X matrix. thing=flacfg
%  matrix - X matrix. thing=flacfg
%  amatrix - identity. thing=flacfg
%  parseline - parses the line, thing = line
%  createline - create a model line. thing=flacfg
%  (autopsd not available for fir)
%
% thing - line to be parsed or flacfg (see fast_flacfg_struct).
%
% FIR Parameters:
%  1. EventId
%  2. PSDMin
%  3. dPSD
%  4. PSDMax
%  5. BCW
%
%


%
% fast_fxcfg_fir.m
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
 
 case {'irfmatrix'}
  % BUG: breaks for ERF
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  nreg = fast_fxcfg_fir('nregressors',flacfg);
  if(isempty(nreg)) return; end
  rt = eye(nreg);
 
 case {'irftaxis','erftaxis'}
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_taxis(flacfg);
 
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
  
 case 'createline'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = createline(flacfg);
 
 otherwise
  fprintf('ERROR: fast_fxcfg_fir: DoWhat = %s, unrecognized\n',DoWhat);

end

return;

%------------------------------------------------------------%
function pr_fla_needed(DoWhat)
fprintf('ERROR: flacfg needed with %s\n',DoWhat);
return;

%------------------------------------------------------------%
function fxcfg = parseline(line)
% Read and check input line
% InputLine: Effect Fixed Label FIR EvId PSDMin dPSD PSDMax BCW
% Parameters:
%  1. EventId
%  2. PSDMin
%  3. dPSD
%  4. PSDMax
%  5. BCW

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
if(rem(evid,1) ~= 0)
  fprintf('ERROR: fir: evid=%g, must be integer \n',evid);
  fxcfg = [];
  return;
end

psdwin = fxcfg.params(2:5);
if(fast_psdwin(psdwin) ~= 1)
  fxcfg = [];
  return;
end

return;

%-----------------------------------------------------------%
function line = createline(flacfg)
line = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

line = sprintf('Effect %s %s %s %d %g %g %g %g',...
	       fxcfg.fxtype,fxcfg.label,fxcfg.model,...
	       fxcfg.params(1),fxcfg.params(2),...
	       fxcfg.params(3),fxcfg.params(4),...
	       fxcfg.params(5));
return;


%------------------------------------------------------------%
function nr = get_nregressors(flacfg)

nr = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

psdwin = fxcfg.params(2:5);
nr = fast_psdwin(psdwin,'npsdwin');

return;

%------------------------------------------------------------%
function  X = get_matrix(flacfg)
X = [];

nr = get_nregressors(flacfg);
if(isempty(nr)) return; end
ntp = fast_fxcfg('getntp',flacfg);
if(isempty(ntp)) return; end
fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

evsch = fast_fxcfg('getevsch',flacfg);
if(isempty(evsch)) return; end

tPres = evsch(:,1);
EvIdPres  = evsch(:,2);
if(size(evsch,2) == 3) wPres = evsch(:,3);
else wPres = [];
end

evid   = fxcfg.params(1);
psdmin = fxcfg.params(2);
dpsd   = fxcfg.params(3);
psdmax = fxcfg.params(4);
bcw    = fxcfg.params(5);
psd = [psdmin dpsd psdmax bcw];
if(fast_psdwin(psd) ~= 1) return; end

indEvId = find(EvIdPres == evid);
tPresEvId = tPres(indEvId);
if(~isempty(wPres)) wPresEvId = wPres(indEvId); 
else wPresEvId = [];
end

if(isempty(flacfg.useevschweight) | flacfg.useevschweight == 0)
  wPresEvId = [];
end

wRun = fast_fxcfg('getrunweight',flacfg);
if(isempty(wRun)) return; end

X = fast_sched2Xfir(tPresEvId,ntp,flacfg.TR,psd,flacfg.tDelay,wPresEvId); 

if(wRun ~= 1) X = wRun * X; end

% tpx and nskip are handled in fast_fxcfg('matrix',flacfg)

return;

%------------------------------------------------------------%
function taxis = get_taxis(flacfg,axistype)
% Ignore axis type as erf and irf are the same
taxis = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

psdmin = fxcfg.params(2);
dpsd   = fxcfg.params(3);
psdmax = fxcfg.params(4);
bcw    = fxcfg.params(5);

taxis = fast_psdwin([psdmin dpsd psdmax bcw],'irftaxis');

return;




