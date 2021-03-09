function flac = fast_ldflac(flacfile,flac)
% flac = fast_ldflac(flacfile,<flac>)
%
% Loads an fsfast flac file.
% If no args, returns an empty flac structure.
%
%
%


%
% fast_ldflac.m
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

if(nargin < 0 | nargin > 2)
  fprintf('flac = fast_ldflac(flacfile,<flac>)\n');
  flac = [];
  return;
end

if(~exist('flac','var')) flac = []; end
if(isempty(flac))
  flac.creator = 'unknown';
  flac.name = '';
  flac.fsd = '';
  flac.runlistfile = '';
  flac.basestem = '';
  flac.funcstem = '';
  flac.runlistfile = '';
  flac.TR = [];
  flac.mask = '';
  flac.inorm = 0;
  flac.tfilter = []; % Spec for temporal filtering matrix
  flac.TFmtx = []; % Actual matrix
  flac.stimulusdelay = 0;
  flac.fixacf = 1;
  flac.acfsegstem = '';
  flac.format = '';
  flac.formatext = '';
  flac.tpexcfile = '';
  flac.nskip = 0;
  flac.parfile  = '';
  flac.schdir = ''; % dir where schedules or pars can be found
  flac.par = [];
  flac.tpexc = [];
  flac.AllowMissingCond = 0;
  flac.autostimdur = 0; % Compute stim duration from par, TR otherwise
  % flag indicating the presence of a variable regressor EV
  flac.varregev = 0; % Does not actually do anything yet
  % VarRegEVs must be nuissance!
  %flac.ev  = []; % Leave commented for inherit
  %flac.con = []; % Leave commented for inherit
  flac.inheritlevel = 0;
  flac.acfbins = 10; % Set to 0 for no whitening
  flac.acffwhm = 0; % Set to 0 for no smoothing
  flac.acfsvd = 2; % Preproc residuals with svd before ar1 calc
  flac.hpfCutoffHz = [];
  flac.HeteroGCor = 0;
  % Allows flac TR and data TR to be different (will use flac TR).
  flac.OverrideTR = 0; 
  flac.fsv3_st2fir = 0;
  flac.fsv3_whiten = 0;  
  flac.RefEventDur = 1;
  flac.rawfwhm = 0;
  flac.volsurffwhm = [];
  flac.rawstc = '';
  flac.sdf = '';
  flac.RegDOF = 0;
  flac.RawSpaceType = '';
  flac.RawSpaceRes = [];
  flac.RawSpace = '';
  flac.PerSession = []; % Binary, 1 if persession
  flac.ExpKey = '';
  % subject surface analysis is being performed on. Can be 'self'.
  flac.subject = ''; 
  flac.hemi = '';
  % Actual subject name of subject being analyzed. This is the name
  % of the subject even if 'self' or fsaverage
  flacnew.sourcesubject = '';
  flac.mc = []; % For motion correction regressors
  flac.globalmean = [];
  flac.IsRetinotopy = 0;
  flac.IsABBlocked = 0;
  flac.stimtype = ''; % for retinotopy
  flac.direction = ''; % for retinotopy
  inherit = 0;
  ana.analysis     = '';
  ana.designtype   = '';
  ana.nconditions  = [];
  ana.PolyOrder    = [];
  ana.extreg       = '';
  ana.nextreg      = [];
  ana.firfit       = [];
  ana.timewindow   = [];
  ana.prestim      = [];
  ana.TER          = [];
  ana.gammafit     = [];
  ana.gamdelay     = [];
  ana.gamtau       = [];
  ana.gamexp       = [];
  ana.spmhrffit    = [];
  ana.nspmhrfderiv = [];
  ana.ncycles      = [];
  ana.nregressors  = [];
  ana.ConditionNames = '';
  ana.info = '';
  flac.ana = ana;
else
  flac.inheritlevel = flac.inheritlevel + 1;
  if(flac.inheritlevel > 10)
    fprintf('ERROR: maximum INHERIT level exceeded.\n');
    fprintf('       Check for circular INHERIT statements.\n');
    flac = [];
    return;
  end
  inherit = 1;
end
if(nargin == 0) return; end
fp = fopen(flacfile,'r');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',flacfile);
  flac = [];
  return;
end

% Derive name from flac file without .flac extension.
flac.name = basename(flacfile,'.flac');

% Get the first line, should look like:
% FSFAST-FLAC 1
tline = fgetl(fp);
flacid = sscanf(tline,'%s',1);
if(~strcmp(flacid,'FSFAST-FLAC'))
  fprintf('ERROR: format error (id string) in flac file %s\n',flacfile);
  fclose(fid);
  flac = [];
  return;
end
flacversion = sscanf(tline,'%*s %d',1);
if(flacversion ~= 1)
  fprintf('ERROR: flac file %s, version = %d\n',flacfile,flacversion);
  fclose(fid);
  flac = [];
  return;
end

nthline = 1;
nthev = 1;  % Keep track of the number of EVs
while(1)
  % scroll through any blank lines or comments
  while(1)
    tline = fgetl(fp);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  key = sscanf(tline,'%s',1);
  %fprintf('key = %s\n',key);
  
  switch(key)
   case 'flacname',    junk             = sscanf(tline,'%*s %s',1);
   case 'fsd',         flac.fsd         = sscanf(tline,'%*s %s',1);
   case 'TR',          flac.TR          = sscanf(tline,'%*s %f',1);
   case 'funcstem',    flac.funcstem    = sscanf(tline,'%*s %s',1);
   case 'mask',        flac.mask        = sscanf(tline,'%*s %s',1);
   case 'inorm',       flac.inorm       = sscanf(tline,'%*s %f',1);
   case 'runlistfile', flac.runlistfile = sscanf(tline,'%*s %s',1);
   case 'stimulusdelay', flac.stimulusdelay  = sscanf(tline,'%*s %f',1);
   case 'StimulusDelay', flac.stimulusdelay  = sscanf(tline,'%*s %f',1);
   case 'acfbins',     flac.acfbins     = sscanf(tline,'%*s %d',1);
   case 'acffwhm',     flac.acffwhm     = sscanf(tline,'%*s %f',1);
   case 'acfsvd',      flac.acfsvd      = sscanf(tline,'%*s %d',1);
   case 'acffix',      flac.fixacf      = sscanf(tline,'%*s %d',1);
   case 'fixacf',      flac.fixacf      = sscanf(tline,'%*s %d',1);
   case 'tpexclude',   flac.tpexcfile   = sscanf(tline,'%*s %s',1);
   case 'parfile',     flac.parfile     = sscanf(tline,'%*s %s',1);
   case 'schdir',      flac.schdir      = sscanf(tline,'%*s %s',1);
   case 'acfseg',      flac.acfsegstem  = sscanf(tline,'%*s %s',1);
   case 'RefEventDur', flac.RefEventDur = sscanf(tline,'%*s %f',1);
   case 'surface', 
    flac.subject = sscanf(tline,'%*s %s',1);
    flac.hemi = sscanf(tline,'%*s %*s %s',1);
   case 'UseTalairach',flac.UseTalairach = 1;
   case 'INHERIT',     
    inheritflacname  = sscanf(tline,'%*s %s',1);
    flacdir = fast_dirname(flacfile);
    inheritfile = sprintf('%s/%s',flacdir,inheritflacname);
    if(strcmp(flacfile,inheritfile))
      fprintf('ERROR: flac file %s cannot INHERIT itself.\n',flacfile);
      flac = [];
      return;
    end
    flac = fast_ldflac(inheritfile,flac);
    if(isempty(flac)) return; end
   case 'FORMAT',      
    flac.format  = sscanf(tline,'%*s %s',1);
    if(strcmp(flac.format,'bvolume')) flac.formatext = ''; 
    elseif(strcmp(flac.format,'mgh')) flac.formatext = '.mgh'; 
    elseif(strcmp(flac.format,'mgz')) flac.formatext = '.mgz'; 
    elseif(strcmp(flac.format,'nii')) flac.formatext = '.nii'; 
    elseif(strcmp(flac.format,'nii.gz')) flac.formatext = '.nii.gz'; 
    else
      fprintf('ERROR: format %s unrecognized\n',flac.format);
      fprintf('line %d\n',nthline);
      flac=[]; return; 
    end
   case 'TFILTER', 
    tfilter = flac_tfilter_parse(tline);
    if(isempty(tfilter)) 
      flac=[]; 
      fprintf('line %d\n',nthline);
      return; 
    end
    flac.tfilter = tfilter;
   case 'EV', 
    ev = flac_ev_parse(tline);
    if(isempty(ev)) 
      flac=[]; 
      fprintf('line %d\n',nthline);
      return; 
    end
    flac.ev(nthev) = ev;
    nthev = nthev + 1;
   case 'CONTRAST-BEGIN', 
    flac = load_contrast(fp,flac);
    if(isempty(flac)) 
      fprintf('line %d\n',nthline);
      return; 
    end
   otherwise
    fprintf('INFO: key %s unrecognized, line %d, skipping\n',key,nthline);
  end
  nthline = nthline + 1;
end % while (1)

fclose(fp);

if(isempty(flac.funcstem)) 
  fprintf('ERROR: no funcstem specified in %s\n',flacfile);
  flac = [];
end
nevs = length(flac.ev);
if(nevs == 0)
  fprintf('ERROR: no EVs specified in %s\n',flacfile);
  flac = [];
end


if(isempty(flac.fsd)) flac.fsd = 'bold'; end 
if(isempty(flac.acfsegstem)) flac.acfsegstem = 'acfseg'; end 

% Check each contrast
if(~isfield(flac,'con')) flac.con = []; end
ncon = length(flac.con);
if(ncon == 0)
  fprintf('WARNING: no contrasts in FLAC file %s\n',flacfile);
  return;
end
for nthcon = 1:ncon

  % Make sure that each EV in the contrast is an EV in the FLAC
  % Handle variable regressor EVs too
  for nthev = 1:length(flac.con(nthcon).ev)
    % This might not work with VarRegEV
    evindex = flac_evindex(flac,flac.con(nthcon).ev(nthev).name);
    if(isempty(evindex))
      fprintf('ERROR: in FLAC file %s, Contrast %s\n',...
	      flacfile,flac.con(nthcon).name);
      fprintf('Contrast EV %s is not in the model\n',...
	      flac.con(nthcon).ev(nthev).name);
      flac = [];
      return;
    end
  end
  
  % Compute the contrast matrices, unless variable reg EV
  if(~flac.varregev)
    flactmp = flac_conmat(flac,nthcon);
    if(isempty(flactmp))
      fprintf('ERROR: with contrast %s in %s\n',...
	      flac.con(nthcon).name,flacfile);
      flac = [];
      return;
    end
    flac = flactmp;
  end
end

% Should do some tests here to make sure names are not rep,
% contrasts, etc



return;

%------------------------------------------------%
% If one per-evrw is spec, then all must be spec (why?)
% If global evrw is spec, then per-evrw not allowed (why?)
% If global evrw is spec, then all non-zero EVs must have same
%  number of regressors equal to the number of weights
% Per-EVRW is not allowed for F-tests
% If Per-EVRW is not spec, then defaults to all ones.
% If Global-EVRW are not spec, then defaults to all ones.
% The number of weights in EVRW for an EV be same as nreg
% Check for variable regressor EVs
% Functional forms for EVRW?
% ANOVA?
function flac = load_contrast(fp,flac)
  % Just loads the contrast spec, does not compute
  % contrast matrix.
  if(~isfield(flac,'con')) ncon = 0; 
  else ncon = length(flac.con);
  end
  ncon = ncon + 1;

  flac.con(ncon).name     = '';
  flac.con(ncon).varsm    = 0;
  flac.con(ncon).sumev    = 0;
  flac.con(ncon).sumevreg = 0;
  flac.con(ncon).sumevrw  = [];
  flac.con(ncon).rmprestim  = 0; % 4/8/09
  
  nthev = 1;
  while(1)
    % scroll through any blank lines or comments
    while(1)
      tline = fgetl(fp);
      if(~isempty(tline) & tline(1) ~= '#') break; end
    end
    if(tline(1) == -1) break; end

    key = sscanf(tline,'%s',1);
    switch(key)
     case 'NAME'
      % Contrast Name
      [tmp c] = sscanfitem(tline,2);
      if(c ~= 1) fprintf('FLAC-CON format error\n');flac=[];return;end
      flac.con(ncon).name = tmp;     
     case 'VARSM'
      % Variance smoothing (?)
      [tmp c] = sscanfitem(tline,2);
      if(c ~= 1) fprintf('FLAC-CON format error: VARSM\n');flac=[];return;end
      flac.con(ncon).varsm = sscanf(tmp,'%f',1);
     case 'SUMEV'
      % SUMEV sumevflag
      % sumevflag=0 : do not sum EVs (give each EV a separate set
      % of rows in C)
      [tmp c] = sscanfitem(tline,2);
      if(c ~= 1) fprintf('FLAC-CON format error\n');flac=[];return;end
      flac.con(ncon).sumev = sscanf(tmp,'%d',1);
     case 'SUMEVREG'
      % SUMEVREG sumevregflag
      % sumevregflag=0 : do not sum regressors within an EV (give
      % each regressor a separate row in C)
      [tmp c] = sscanfitem(tline,2);
      if(c ~= 1) fprintf('FLAC-CON format error\n');flac=[];return;end
      flac.con(ncon).sumevreg = sscanf(tmp,'%d',1);
     case 'EV'
      % EV evname evweight <evreg1w evreg2w ... evregNREGw>
      % evweight - overall weight for a regressor
      % evregNw - weight for Nth regressor (overrides EVRW)
      % Cannot spec weights for EV Regressors when there are a
      % variable number of regressors in the EV
      [tmp c] = sscanfitem(tline,2);
      if(c ~= 1) fprintf('FLAC-CON format error\n');flac=[];return;end
      flac.con(ncon).ev(nthev).name = tmp;
      [tmp c] = sscanfitem(tline,3);
      if(c ~= 1) fprintf('FLAC-CON format error\n');flac=[];return;end
      flac.con(ncon).ev(nthev).evw = sscanf(tmp,'%f',1);
      nthevrw = 0;
      while(1)
	[tmp c] = sscanfitem(tline,nthevrw+4);
	if(~c) break; end
	nthevrw = nthevrw + 1;
	flac.con(ncon).ev(nthev).evrw(nthevrw) = sscanf(tmp,'%f',1);
      end
      if(nthevrw == 0) flac.con(ncon).ev(nthev).evrw = []; end
      nthev = nthev+1;
     case 'EVRW' % Default EV Regressor Weights
      % EVRW evreg1w evreg2w ... evregNREGw
      % Weights to use for all EVs in contrast. 
      % Can be overriden if EV Reg weights are specified with the EV.
      % EVs must have nreg = NREG
      % Cannot spec weights for EV Regressors when there are a
      % variable number of regressors in the EV
      nthevrw = 0;
      while(1)
	[tmp c] = sscanfitem(tline,nthevrw+2);
	if(~c) break; end
	nthevrw = nthevrw + 1;
	flac.con(ncon).evrw(nthevrw) = sscanf(tmp,'%f',1);
      end
      if(nthevrw == 0)
	fprintf('ERROR: EVRW has no weights\n');
	flac = [];
	return;
      end
     case 'CONTRAST-END'
      if(~isfield(flac.con(ncon),'evrw')) flac.con(ncon).evrw = []; end
      return;
     otherwise
      fprintf('ERROR: key %s in contrast %d not recognized\n',key,ncon);
      flac = [];
      return;
    end
  end
  

return

