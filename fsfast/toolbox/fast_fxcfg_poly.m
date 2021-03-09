function rt = fast_fxcfg_poly(DoWhat,thing)
% rt = fast_fxcfg_poly(DoWhat,thing)
%
% DoWhat can be:
%  iserm  - returns 0 (this is not an erm) thing not needed.
%  nparams - returns number of parameters in model (1). thing not needed.
%  nregressors - number of regressors in current X matrix. thing=flacfg
%  matrix - X matrix. thing=flacfg
%  parseline - parses the line, thing=line
%  createline - create a model line. thing=flacfg
%  amatrix - not appliciable
%  autopsd - not appliciable
%
% thing - line to be parsed or flacfg (see fast_flacfg_struct).
%
% Polynomial Parameters:
%  1. Order
%


%
% fast_fxcfg_poly.m
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
  fprintf('rt = fast_fxcfg_poly(DoWhat,<thing>)\n');
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
  rt = 0;
 
 case 'nparams'
  rt = 1;
 
 case 'nregressors'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_nregressors(flacfg);
 
 case 'autopsd'
  fprintf('ERROR: poly: cannot autopsd\n');
  rt = []; 
 case 'amatrix'
  fprintf('ERROR: poly: cannot amatrix\n');
  rt = []; 
 
 case 'matrix'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  X = get_matrix(flacfg);
  if(isempty(X)) return; end
  rt = X;
 
 case 'parseline'
  if(isempty(line))
    fprintf('ERROR: poly: line needed with parseline\n');
    return;
  end
  rt = parseline(line);
  
 case 'createline'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = createline(flacfg);
 
 otherwise
  fprintf('ERROR: fast_fxcfg_poly: getwhat = %s, unrecognized\n',getwhat);

end

return;

%------------------------------------------------------------%
function pr_fla_needed(DoWhat)
fprintf('ERROR: flacfg needed with %s\n',DoWhat);
return;

%------------------------------------------------------------%
function fxcfg = parseline(line)
% Read and check input line
% InputLine: Effect FxType Label polynomial order
% Parameters:
%  1. Order

fxcfg = [];

nparams = fast_fxcfg_poly('nparams');
[tmp nitems] = sscanf(line,'%s',inf);
if(nitems ~= nparams + 4)
  fprintf('ERROR: poly: line has wrong number of items (%d)\n',nitems);
  fprintf('%s\n',line);
  return;
end

fxcfg = fast_fxcfg_struct;

[fxcfg.fxtype n]  = sscanf(line,'%*s %s ',1);
[fxcfg.label  n]  = sscanf(line,'%*s %*s %s',1);
[fxcfg.model  n]  = sscanf(line,'%*s %*s %*s %s',1);

fxcfg.fxtype = lower(fxcfg.fxtype);
fxcfg.model = lower(fxcfg.model);

if(~strcmp(fxcfg.model,'polynomial'))
  fprintf('ERROR: poly: not poly model (%s)\n',fxcfg.model);
  fprintf('%s\n',line);
  fxcfg = [];
  return;
end

order = sscanf(line,'%*s %*s %*s %*s %d',1);

if(order < 0) 
  fprintf('ERROR: poly: order = %d, must be >= 0',order);
  fxcfg = [];
  return;
end

fxcfg.params(1) = order;

return;

%-----------------------------------------------------------%
function line = createline(flacfg)
line = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

line = sprintf('Effect %s %s %s %d %g %g %g\n',...
	       fxcfg.fxtype,fxcfg.label,fxcfg.model,fxcfg.params(1));
return;

%------------------------------------------------------------%
function nr = get_nregressors(flacfg)

nr = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

order = fxcfg.params(1);
nr = order + 1;

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

order = fxcfg.params(1);
X = fast_polytrendmtx(1,ntp,1,order);

% tpx and nskip are handled in fast_fxcfg('matrix',flacfg)

return;
%------------------------------------------------------------%





