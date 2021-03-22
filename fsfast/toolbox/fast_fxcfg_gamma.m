function rt = fast_fxcfg_gamma(DoWhat,thing)
% rt = fast_fxcfg_gamma(DoWhat,thing)
%
% DoWhat can be:
%  iserm  - returns 1 because this is an erm. thing not needed.
%  nparams - returns number of parameters in model (9). thing not needed.
%  nregressors - number of regressors in current X matrix. thing=flacfg
%  matrix - X matrix. thing=flacfg
%  irfmatrix - thing=flacfg
%  erfmatrix - thing=flacfg
%  irftaxis, erftaxis
%  parseline - parses the line, thing = line
%  createline - create a model line. thing=flacfg
%  autopsd - returns psdwindow
%
% thing - line to be parsed or flacfg (see fast_flacfg_struct).
%
% Gamma Parameters:
%  1. EventId
%  2. PSDMin
%  3. dPSD
%  4. PSDMax
%  5. BoxCarWidth
%  6. AutoPSD
%  7. Delta
%  8. Tau
%  9. Number of Derivatives to add
%
%


%
% fast_fxcfg_gamma.m
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
  fprintf('rt = fast_fxcfg_gamma(DoWhat,<thing>)\n');
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
  rt = 9;
 
 case 'nregressors'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_nregressors(flacfg);
 
 case 'autopsd'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = [0 .1 40]; % need to do better here
 
 case {'irfmatrix','erfmatrix'}
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_amatrix(flacfg,DoWhat);
 
 case 'irftaxis'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_taxis(flacfg,'irftaxis');

 case 'erftaxis'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_taxis(flacfg,'erftaxis');

 case 'matrix'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  X = get_matrix(flacfg);
  if(isempty(X)) return; end
  rt = X;
 
 case 'parseline'
  if(isempty(line))
    fprintf('ERROR: gamma: line needed with parseline\n');
    return;
  end
  rt = parseline(line);
  
 case 'createline'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = createline(flacfg);
 
 otherwise
  fprintf('ERROR: fast_fxcfg_gamma: DoWhat = %s, unrecognized\n',DoWhat);

end

return;

%------------------------------------------------------------%
function pr_fla_needed(DoWhat)
fprintf('ERROR: flacfg needed with %s\n',DoWhat);
return;

%------------------------------------------------------------%
function fxcfg = parseline(line)
% Read and check input line
% InputLine: Effect Fixed Label Gamma EvId PSDMin dPSD PSDMax ...
% Parameters:
%  1. EventId
%  2. PSDMin
%  3. dPSD
%  4. PSDMax
%  5. BoxCarWidth
%  6. AutoPSD
%  7. Delta
%  8. Tau
%  9. Number of Derivatives to add

fxcfg = [];

nparams = fast_fxcfg_gamma('nparams');
[tmp nitems] = sscanf(line,'%s',inf);
if(nitems ~= nparams + 4)
  fprintf('ERROR: gamma: line has wrong number of items (%d)\n',nitems);
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
  fprintf('ERROR: gamma: must be fixed effect type\n');
  fprintf('%s\n',line);
  fxcfg = [];
  return;
end
if(~strcmp(fxcfg.model,'gamma'))
  fprintf('ERROR: gamma: not gamma model\n');
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
if(rem(evid,1) ~= 0)
  fprintf('ERROR: gamma: evid=%g, must be integer \n',evid);
  fxcfg = [];
  return;
end
if(evid < 1)
  fprintf('ERROR: gamma: evid=%d < 1\n',evid);
  fxcfg = [];
  return;
end

autopsdwin = fxcfg.params(6);
if(autopsdwin ~= 1)
  psdwin = [fxcfg.params(2:5)];
  if(fast_psdwin(psdwin) ~= 1)
    fxcfg = [];
    return;
  end
end

return;

%-----------------------------------------------------------%
function line = createline(flacfg)
line = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

line = sprintf('Effect %s %s %s',fxcfg.fxtype,fxcfg.label,fxcfg.model);
for n = 1:length(fxcfg.params);
  line = sprintf('%s %g',line,fxcfg.params(n));
end

return;

%------------------------------------------------------------%
function nr = get_nregressors(flacfg)

nr = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

nderiv = fxcfg.params(9);
nr = nderiv + 1;

return;

%------------------------------------------------------------%
function taxis = get_taxis(flacfg,axistype)
taxis = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

autopsdwin = fxcfg.params(6);
if(autopsdwin)
  psdwin = fast_fxcfg('autopsd',flacfg);
  if(isempty(psdwin)) return; end
  psdmin = psdwin(1);
  dpsd   = psdwin(2);
  psdmax = psdwin(3);
else
  psdmin = fxcfg.params(2);
  dpsd   = fxcfg.params(3);
  psdmax = fxcfg.params(4);
end
bcw = fxcfg.params(5);

taxis = fast_psdwin([psdmin dpsd psdmax bcw],axistype);

return;
%------------------------------------------------------------%
function A = get_amatrix(flacfg,matrixtype)
% matrixtype = 'irfmatrix' or 'erfmatrix'

A = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

% Get IRF time axis here, apply BCW below
t = get_taxis(flacfg,'irftaxis');
if(isempty(t)) return; end

bcw    = fxcfg.params(5);
delta  = fxcfg.params(7);
tau    = fxcfg.params(8);
nderiv = fxcfg.params(9);

r = (t-delta)./tau;
h = (r.^2) .* exp(-r);
indz = find(t<=delta);
h(indz) = 0;
h = h/max(h); % Do I really want to do this? Makes compat with
              % fmri_hemodyn()
hd = h;
for n = 1:nderiv
  % Should probably do real derivative here %
  h = [h gradient(hd)];
  hd = gradient(hd);
end

A = h;

return;

%------------------------------------------------------------%
function  X = get_matrix(flacfg)
X = [];

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

autopsdwin = fxcfg.params(6);
if(autopsdwin)
  psdwin = fast_fxcfg('autopsd',flacfg);
  if(isempty(psdwin)) return; end
  psdmin = psdwin(1);
  dpsd   = psdwin(2);
  psdmax = psdwin(3);
else
  psdmin = fxcfg.params(2);
  dpsd   = fxcfg.params(3);
  psdmax = fxcfg.params(4);
end
bcw = fxcfg.params(5);
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

% Build any BCW into Xfir
Xfir = fast_sched2Xfir(tPresEvId,ntp,flacfg.TR,psd,flacfg.tDelay,wPresEvId); 
if(isempty(Xfir)) return; end

% Get IRF
A = fast_fxcfg('irfmatrix',flacfg);
if(isempty(A)) return; end

if(size(Xfir,2) ~= size(A,1))
  fprintf('ERROR: gamma: dimension mismatch between Xfir and A\n');
  return;
end

X = Xfir*A;

% tpx and nskip are handled in fast_fxcfg('matrix',flacfg)

return;
%------------------------------------------------------------%





