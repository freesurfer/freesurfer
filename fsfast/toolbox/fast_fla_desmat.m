function [X, flacfg] = fast_fla_desmat(flacfg,mthrun,mthruntype)
% [X flacfg] = fast_fla_desmat(flacfg,<mthrun,mthruntype>)
% mthruntype = 0 or [] (all runs), 1 (perrun), 2 (jkrun)
% Also assigns fxcfg.regind for each fx/run but only in the
% returned flacfg.


%
% fast_fla_desmat.m
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

X = [];

if(nargin < 1 | nargin > 3)
  fprintf('X = fast_fla_desmat(flacfg,<mthrun,mthruntype>)\n');
  return;
end
if(nargin == 2)
  fprintf('ERROR: must specify mthruntype with mthrun\n');
  return;
end
if(exist('mthrun') ~= 1) mthrun = []; end
if(isempty(flacfg.sesscfg.runlist))
  fprintf('ERROR: no runs in sesscfg run list\n');
  return;
end
if(isempty(flacfg.fxlist))
  fprintf('ERROR: no fx defined in fla.\n');
  return;
end

nruns = size(flacfg.sesscfg.runlist,1);
nfx = length(flacfg.fxlist);

% Get the list of run indices to use %
if(~isempty(mthrun))
  if(mthrun > nruns) 
    fprintf('ERROR: mthrun = %d > nruns = %d\n',mthrun,nruns);
    return;
  end
  if(mthruntype == 0) mthrunlist = 1:nruns; end % do all runs
  if(mthruntype == 1) mthrunlist = mthrun; end  % per run 
  if(mthruntype == 2)                           % jk run 
    mthrunlist = 1:nruns;
    ind = find(mthrunlist~=mthrun);
    mthrunlist = mthrunlist(ind);
  end
else
  mthrunlist = 1:nruns; 
end  

% Get the total number of fixed effects regressors
nfixedreg = fast_fxcfg('nfixedreg',flacfg);

% save in order to restore before return
nthrunsv = flacfg.nthrun;
nthfxsv = flacfg.nthfx;

% Build the Design Matrix
Xfe = []; % fixed  part
Xre = []; % random part
for nthrun = mthrunlist

  % build each run separately
  flacfg.nthrun = nthrun;
  Xferun = [];
  Xrerun = [];

  for nthfx = 1:nfx
    % build each effect 
    flacfg.nthfx = nthfx;
    fx = fast_fxcfg('getfxcfg',flacfg);
    if(isempty(fx)) return; end

    % Get matrix for this effect
    Xtmp = fast_fxcfg('matrix',flacfg);
    if(isempty(Xtmp)) return; end
    nreg = size(Xtmp,2);
    
    % Append it to the proper type
    if(strcmp(fx.fxtype,'fixed'))  
      flacfg.fxlist(nthfx).fx.regind{1} = [1:nreg]' + size(Xferun,2);
      Xferun = [Xferun Xtmp];  
    end
    if(strcmp(fx.fxtype,'random')) 
      flacfg.fxlist(nthfx).fx.regind{nthrun} = ...
	  [1:nreg]' + size(Xrerun,2) + nfixedreg;
      Xrerun = [Xrerun Xtmp]; 
    end
  end

  % Accumulated fixed effects
  Xfe = [Xfe; Xferun]; % Should check size here

  % Accumulated random effects
  Za = zeros(size(Xre,1),size(Xrerun,2)); % zero padding
  Zb = zeros(size(Xrerun,1),size(Xre,2)); % zero padding
  Xre = [Xre Za; Zb Xrerun];
  
end

% Finally!
X = [Xfe Xre];

% restore 
flacfg.nthrun = nthrunsv;
flacfg.nthfx = nthfxsv;

return;
