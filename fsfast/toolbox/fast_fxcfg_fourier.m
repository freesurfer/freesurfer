function rt = fast_fxcfg_fourier(DoWhat,thing)
% Handle configuration for a fourier regressor. This is
% not an Event Response Model.
%
% rt = fast_fxcfg_spmhrf(DoWhat,thing)
%
% DoWhat can be:
%  iserm  - returns 0 because this is not an erm. thing not needed.
%  nparams - returns number of parameters in model (3). thing not needed.
%  nregressors - number of regressors in current X matrix. thing=flacfg
%  matrix - X matrix. thing=flacfg
%  irfmatrix - not applicable
%  erfmatrix - not applicable
%  irftaxis, erftaxis - not applicable
%  parseline - parses the line, thing = line
%  createline - create a model line. thing=flacfg
%  autopsd - not applicable
%
% thing - line to be parsed or flacfg (see fast_flacfg_struct).
%
% Fourier Parameters:
%  1. Period of fundamental (sec)
%  2. Number of Harmonics
%  3. Delay (sec) - example: slice delay
%
% Need a way to do reversals.
% Need mechanism to compute phase.
%
%


%
% fast_fxcfg_fourier.m
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
  fprintf('rt = fast_fxcfg_fourier(DoWhat,<thing>)\n');
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
  rt = 3;
 
 case 'nregressors'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = get_nregressors(flacfg);
 
 case {'autopsd','irfmatrix','erfmatrix','irftaxis','erftaxis','amatrix'}
  fprintf('ERROR: fourier: flag %s not appropriate\n',DoWhat);
  rt = []; 

 case 'matrix'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  X = get_matrix(flacfg);
  if(isempty(X)) return; end
  rt = X;
 
 case 'parseline'
  if(isempty(line))
    fprintf('ERROR: spmhrf: line needed with parseline\n');
    return;
  end
  rt = parseline(line);
  
 case 'createline'
  if(isempty(flacfg)) pr_fla_needed(DoWhat); return; end
  rt = createline(flacfg);
 
 otherwise
  fprintf('ERROR: fast_fxcfg_fourier: DoWhat = %s, unrecognized\n',DoWhat);

end

return;

%------------------------------------------------------------%
function pr_fla_needed(DoWhat)
fprintf('ERROR: flacfg needed with %s\n',DoWhat);
return;

%------------------------------------------------------------%
function fxcfg = parseline(line)
% Read and check input line
% InputLine: Effect Fixed/Random Label fourier Period NHarmonics Delay
% Parameters:
%  1. Period (sec)
%  2. NHarmonics
%  3. Delay (sec)

fxcfg = [];

nparams = fast_fxcfg_fourier('nparams');
[tmp nitems] = sscanf(line,'%s',inf);
if(nitems ~= nparams + 4)
  fprintf('ERROR: fourier: line has wrong number of items (%d)\n',nitems);
  fprintf('%s\n',line);
  return;
end

fxcfg = fast_fxcfg_struct;

[fxcfg.fxtype n]  = sscanf(line,'%*s %s ',1);
[fxcfg.label  n]  = sscanf(line,'%*s %*s %s',1);
[fxcfg.model  n]  = sscanf(line,'%*s %*s %*s %s',1);

fxcfg.fxtype = lower(fxcfg.fxtype);
fxcfg.model = lower(fxcfg.model);

if(~strcmp(fxcfg.model,'fourier'))
  fprintf('ERROR: fourier: not fourier model\n');
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

period     = fxcfg.params(1);
nharmonics = fxcfg.params(2);
delay      = fxcfg.params(3);

if(period <= 0) 
  fprintf('ERROR: period = %g, cannot be <= 0\n',period);
  fxcfg = [];
  return;
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

nharmonics = fxcfg.params(2);
nr = 2 + 2*nharmonics;

return;

%------------------------------------------------------------%
function  X = get_matrix(flacfg)
X = [];

ntp = fast_fxcfg('getntp',flacfg);
if(isempty(ntp)) return; end
fxcfg = fast_fxcfg('getfxcfg',flacfg);
if(isempty(fxcfg)) return; end

period     = fxcfg.params(1);
nharmonics = fxcfg.params(2);
delay      = fxcfg.params(3);
t = flacfg.TR * [0:ntp-1]' + delay;
X = [];
for n = 1:nharmonics+1
  xc = cos(n*2*pi*t/period);
  xs = sin(n*2*pi*t/period);
  X = [X xc xs];
end

% tpx and nskip are handled in fast_fxcfg('matrix',flacfg)

return;
%------------------------------------------------------------%





