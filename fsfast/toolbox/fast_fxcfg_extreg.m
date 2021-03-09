function rt = fast_fxcfg_extreg(DoWhat,thing)
% rt = fast_fxcfg_extreg(DoWhat,thing)
%
% DoWhat can be:
%  iserm  - returns 0 because this is not an erm. thing not needed.
%  nparams - returns number of parameters in model. thing not needed.
%  nregressors - number of regressors in current X matrix. thing=flacfg
%  matrix - X matrix. thing=flacfg
%  amatrix - not applicable.
%  parseline - parses the line, thing = line
%  createline - create a model line. thing=flacfg
%  autopsd - not applicable
%
% thing - line to be parsed or flacfg (see fast_flacfg_struct).
%
% ExtReg Parameters:
%  1. fname - name relative to fsd/RRR (stored in sparams not params)
%  2. nextreg  - number of extreg to use (-1 for all) (stored in params)
%


%
% fast_fxcfg_extreg.m
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
  fprintf('rt = fast_fxcfg_extreg(DoWhat,<thing>)\n');
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

%----------------------------%
DoWhat = lower(DoWhat);
switch(DoWhat)
 
 case 'iserm'
  rt = 0;
 
 case 'nparams'
  rt = 2;
 
 case 'nregressors'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_nregressors(flacfg);
 
 case 'autopsd'
  fprintf('ERROR: extreg: cannot autopsd\n');
  return;
 
 case 'amatrix'
  fprintf('ERROR: extreg: cannot amatrix\n');
  return;
 
 case 'matrix'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  X = get_matrix(flacfg);
  if(isempty(X)) return; end
  rt = X;
 
 case 'parseline'
  if(isempty(line))
    fprintf('ERROR: extreg: line needed with parseline\n');
    return;
  end
  rt = parseline(line);
  
 case 'createline'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = createline(flacfg);
 
 otherwise
  fprintf('ERROR: fast_fxcfg_extreg: getwhat = %s, unrecognized\n',getwhat);

end

return;

%------------------------------------------------------------%
function pr_fla_needed(DoWhat)
fprintf('ERROR: flacfg needed with %s\n',DoWhat);
return;

%------------------------------------------------------------%
function fxcfg = parseline(line)
% Read and check input line
% InputLine: Effect F/R Label ExtReg fname 
% Parameters: 
%  1. fname - name relative to fsd/RRR (stored in sparams)
%  2. nextreg (-1 for all) (stored in params)

fxcfg = [];

nparams = fast_fxcfg_extreg('nparams');
[tmp nitems] = sscanf(line,'%s',inf);
if(nitems ~= nparams + 4)
  fprintf('ERROR: extreg: line has wrong number of items (%d)\n',nitems);
  fprintf('%s\n',line);
  return;
end

fxcfg = fast_fxcfg_struct;

[fxcfg.fxtype n]  = sscanf(line,'%*s %s ',1);
[fxcfg.label  n]  = sscanf(line,'%*s %*s %s',1);
[fxcfg.model  n]  = sscanf(line,'%*s %*s %*s %s',1);

fxcfg.fxtype = lower(fxcfg.fxtype);
fxcfg.model = lower(fxcfg.model);

if(~strcmp(fxcfg.model,'extreg'))
  fprintf('ERROR: extreg: model not extreg (%s)\n',fxcfg.model);
  fprintf('%s\n',line);
  fxcfg = [];
  return;
end

fxcfg.sparams = sscanf(line,'%*s %*s %*s %*s %s',1);
fxcfg.params  = sscanf(line,'%*s %*s %*s %*s %*s %d',1);
return;

%-----------------------------------------------------------%
function line = createline(flacfg)
line = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

line = sprintf('Effect %s %s %s %s %d',fxcfg.fxtype,...
	       fxcfg.label,fxcfg.model,...
	       deblank(fxcfg.sparams(1,:)),fxcfg.params(1));

return;

%------------------------------------------------------------%
function nr = get_nregressors(flacfg)

nr = [];

fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

if(isempty(fxcfg.npmlist(flacfg.nthrun).M))
  fprintf('ERROR: extreg: extreg not loaded\n');
  return;
end

nextreg = fxcfg.params(1);
if(nextreg == -1)
  nr = size(fxcfg.npmlist(flacfg.nthrun).M,2);
else
  nr = nextreg;
end

return;

%------------------------------------------------------------%
function  X = get_matrix(flacfg)
X = [];

ntp = fast_fxcfg('getntp',flacfg);
if(isempty(ntp)) return; end
fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end
if(isempty(fxcfg.npmlist(flacfg.nthrun).M))
  fprintf('ERROR: extreg: extreg not loaded\n');
  return;
end
M = fxcfg.npmlist(flacfg.nthrun).M;
if(size(M,1) ~= ntp)
  fprintf('ERROR: extreg: extreg dimension mismatch\n');
  return;
end

% Really do this?
M = M - repmat(mean(M), [ntp 1]);
M = M./repmat(std(M), [ntp 1]);

nextreg = fxcfg.params(1);
if(nextreg == -1)
  X = M;
else
  if(nextreg > size(M,2))
    fprintf('ERROR: extreg: extreg does not have enough cols\n');
    return;
  end
  X = M(:,1:nextreg);
end

% tpx and nskip are handled in fast_fxcfg('matrix',flacfg)

return;
%------------------------------------------------------------%





