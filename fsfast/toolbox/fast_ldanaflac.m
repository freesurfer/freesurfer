function flac = fast_ldanaflac(anadir)
% flac = fast_ldanaflac(anadir)
%
%


%
% fast_ldanaflac.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2010/04/15 16:56:16 $
%    $Revision: 1.46 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(nargin ~= 1)
  fprintf('flac = fast_ldanaflac(anadir)\n');
  flac = [];
  return;
end

flac = fast_ldflac; % creates empty struct
flac.name = basename(anadir);
flac.AllowMissingCond = 1;
flac.autostimdur = [];
flac.acfbins = 10;
% Default in mkanalysis is to fix
flac.fixacf = 1; 
  
flac.mask = '';
flac.con = [];
% format is handled diff than in fast_ldflac.m
flac.format = getenv('FSF_OUTPUT_FORMAT');
if(isempty(flac.format)) flac.format = 'nii'; end
flac.formatext = sprintf('.%s',flac.format);
%flac.ev = []; % Leave commented?


%----------- Read in the analysis.info -------------------
info = sprintf('%s/analysis.info',anadir);
designtype = 'event-related';
nconditions = [];
ConditionNames = '';
fp = fopen(info,'r');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',info);
  flac = [];
  return;
end
nthline = 1;
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
   case 'analysis',    analysis         = sscanf(tline,'%*s %s',1);
   case 'TR',          flac.TR          = sscanf(tline,'%*s %f',1);
   case 'OverrideTR',  flac.OverrideTR  = sscanf(tline,'%*s %d',1);
   case 'RefEventDur', flac.RefEventDur = sscanf(tline,'%*s %f',1);
   case 'fsd',         flac.fsd         = sscanf(tline,'%*s %s',1);
   case 'funcstem',    flac.funcstem    = sscanf(tline,'%*s %s',1);
   case 'maskstem',    flac.mask        = sscanf(tline,'%*s %s',1);
   case 'inorm',       flac.inorm       = sscanf(tline,'%*s %f',1);
   case 'runlistfile', flac.runlistfile = sscanf(tline,'%*s %s',1);
   case 'tpexclude',   flac.tpexcfile   = sscanf(tline,'%*s %s',1);
   case 'parname',     flac.parfile     = sscanf(tline,'%*s %s',1);
   case 'designtype',  designtype       = sscanf(tline,'%*s %s',1);
   case 'nconditions', nconditions      = sscanf(tline,'%*s %d',1);
   case 'surface',     
    flac.subject = sscanf(tline,'%*s %s',1);
    flac.hemi = sscanf(tline,'%*s %*s %s',1);
   case 'UseTalairach',flac.UseTalairach = 1;
   case 'Condition', 
    conditionid   = sscanf(tline,'%*s %d',1);
    ConditionName = sscanfitem(tline,3);
    ConditionNames = strvcat(ConditionNames,ConditionName);
   otherwise
    fprintf('INFO: key %s unrecognized, line %d, skipping\n',key,nthline);
  end
  nthline = nthline + 1;
end % while (1)
fclose(fp);

fprintf('RED %g\n',flac.RefEventDur);

%----------- Read in the analysis.cfg -------------------
TER = flac.TR;
PolyOrder = 0;
extregList = '';
nextregList = [];
nskip = 0;
ncycles = [];
delay = 0;
timeoffset = 0;
gammafit = 0;
gamdelay = 0;
gamtau = 0;
gamexp = 2;
spmhrffit = 0;
nspmhrfderiv = 0;
timewindow = 0;
prestim = 0;
cfg  = sprintf('%s/analysis.cfg',anadir);
fp = fopen(cfg,'r');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',info);
  flac = [];
  return;
end
nthline = 1;
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
   case '-gammafit',   
    gammafit = 1;
    gamdelay = sscanf(tline,'%*s %f',1);
    gamtau   = sscanf(tline,'%*s %*f %f',1);
   case '-gammaexp',   gamexp       = sscanf(tline,'%*s %f',1);
   case '-spmhrf',     
    nspmhrfderiv = sscanf(tline,'%*s %d',1);
    spmhrffit = 1;
   case '-polyfit',    PolyOrder   = sscanf(tline,'%*s %f',1);
   case '-TER',        TER         = sscanf(tline,'%*s %f',1);
   case '-acfbins',    flac.acfbins = sscanf(tline,'%*s %d',1);
   case '-autostimdur',flac.autostimdur = 1;
   case '-noautostimdur',flac.autostimdur = 0;
   case '-extreg',     
    extreg      = sscanf(tline,'%*s %s',1);
    extregList = strvcat(extregList,extreg);
    nextreg     = sscanf(tline,'%*s %*s %d',1);
    nextregList = [nextregList nextreg];
   case '-nextreg',    nextreg     = sscanf(tline,'%*s %d',1);
   case '-rescale',    flac.inorm  = sscanf(tline,'%*s %f',1);
   case '-nskip',      nskip       = sscanf(tline,'%*s %d',1);
   case '-prestim',    prestim     = sscanf(tline,'%*s %f',1);
   case '-timewindow', timewindow  = sscanf(tline,'%*s %f',1);
   case '-ncycles',    ncycles     = sscanf(tline,'%*s %f',1);
   case '-delay',      delay       = sscanf(tline,'%*s %f',1);
   case '-timeoffset', timeoffset  = sscanf(tline,'%*s %f',1);
   case '-fwhm',       sscanf(tline,'%*s %f',1); % dont worry about it
   case '-fix-acf',    flac.fixacf = 1; 
   case '-no-fix-acf', flac.fixacf = 0;
   case '-fsv3-st2fir',    flac.fsv3_st2fir = 1;
   case '-no-fsv3-st2fir', flac.fsv3_st2fir = 0;
   case '-fsv3-whiten',     flac.fsv3_whiten = 1;
   case '-no-fsv3-whiten',  flac.fsv3_whiten = 0;
   otherwise
    fprintf('INFO: key %s unrecognized, line %d, skipping\n',key,nthline);
  end
  nthline = nthline + 1;
end % while (1)
fclose(fp);

if(isempty(flac.funcstem)) 
  fprintf('ERROR: no funcstem specified in %s\n',flacfile);
  flac = [];
  return;
end

if(isempty(flac.mask)) 
  flac.mask = 'brain';
  fprintf('INFO: mask is not set, setting to brain\n');
end
if(isempty(flac.fsd)) flac.fsd = 'bold'; end 
if(isempty(flac.acfsegstem)) flac.acfsegstem = 'acfseg'; end 

if(timeoffset ~= 0) flac.stimulusdelay = timeoffset; end
if(delay ~= 0)      flac.stimulusdelay = delay; end

ana.analysis     = analysis;
ana.designtype   = designtype;
ana.PolyOrder    = PolyOrder;
ana.extregList   = extregList;
ana.nextregList  = nextregList;
ana.ncycles      = ncycles ;
ana.nconditions  = nconditions;
if(~isempty(ConditionNames))
  ana.ConditionNames = ConditionNames;
else
  ana.ConditionNames = '';
  for n = 1:ana.nconditions
    tmp = sprintf('Condition%02d',n);
    ana.ConditionNames = strvcat(ana.ConditionNames,tmp);
  end
end


if(~spmhrffit & ~gammafit) ana.firfit = 1; 
else                       ana.firfit = 0; 
end
ana.timewindow   = timewindow;
ana.prestim      = prestim;
ana.TER          = TER;

ana.gammafit     = gammafit;
ana.gamdelay     = gamdelay;
ana.gamtau       = gamtau;
ana.gamexp       = gamexp;

ana.spmhrffit    = spmhrffit;
ana.nspmhrfderiv = nspmhrfderiv;

if(ana.gammafit) 
  ana.nregressors = 1; 
end
if(ana.spmhrffit) 
  ana.nregressors = ana.nspmhrfderiv + 1;
end
if(ana.firfit)
  ana.nregressors = round(ana.timewindow/ana.TER);
end

flac.ana = ana;
flac.ana.con = [];

nthev = 1;
tline = sprintf('EV Baseline baseline nuis');
flac.ev(nthev) = flac_ev_parse(tline);
nthev = nthev+1;

if(strcmp(designtype,'event-related') | strcmp(designtype,'blocked'))
  for n=1:nconditions
    if(gammafit)
      tline = sprintf('EV Condition%02d gamma task %d %g %g %g 0 %g',...
		      n,n,gamdelay,gamtau,gamexp,TER);
    elseif(spmhrffit)
      tline = sprintf('EV Condition%02d spmhrf task %d %d %g',...
		      n,n,nspmhrfderiv,TER);
    else
      tline = sprintf('EV Condition%02d fir task %d %g %g %g %g',...
		      n,n,-prestim,timewindow-prestim,TER);
      
    end
    flac.ev(nthev) = flac_ev_parse(tline);
    nthev = nthev+1;
  end
end

if(strcmp(designtype,'abblocked') | strcmp(designtype,'retinotopy'))
  % par file will have:
  %   stimtype  eccen (retinotopy)
  %   direction neg
  % rtopy output will go in ananame/{eccen,polar}
  % Need to add nuis on either side of fund and harm to be
  % compatible with selfreqavg
  period = ncycles * flac.TR;
  nharmonics = 1;
  tline = sprintf('EV Fourier fourier task %g %g %g',...
		  period,nharmonics,delay);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;

  nthcon = 1;
  flac.con(nthcon).name     = 'omnibus';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [1 1 1 1];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.ana.con = flac.con(nthcon);

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'fund';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [1 1 0 0];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.ana.con(nthcon) = flac.con(nthcon);

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'fund-sin';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [1 0 0 0];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.ana.con(nthcon) = flac.con(nthcon);

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'fund-cos';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [0 1 0 0];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.ana.con(nthcon) = flac.con(nthcon);

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'harm';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [0 0 1 1];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.ana.con(nthcon) = flac.con(nthcon);
  
  ncontrasts = length(flac.con);
end

if(PolyOrder > 0)
  tline = sprintf('EV Poly polynomial nuis %d',PolyOrder);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;
end

if(~isempty(extregList))
  nlist = size(extregList,1);
  for n = 1:nlist
    extreg = deblank(extregList(n,:));
    nextreg = nextregList(n);
    extregevname = basename(extreg,'dat'); % remove .dat
    tline = sprintf('EV %s nonpar nuis %s %d',extregevname,extreg,nextreg);
    flac.ev(nthev) = flac_ev_parse(tline);
    nthev = nthev+1;
  end
end

if(~isempty(flac.tpexcfile))
  tline = sprintf('EV TExclude texclude nuis %s',flac.tpexcfile);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;
end


if(strcmp(designtype,'event-related') | strcmp(designtype,'blocked'))
  %-------------- contrasts --------------------------
  tmpstr = sprintf('%s/*.mat',anadir);
  clist = dir(tmpstr);
  ncontrasts = length(clist);
  for nthcon = 1:ncontrasts
    fprintf('%2d %s\n',nthcon,clist(nthcon).name);
    tmpstr = sprintf('%s/%s',anadir,clist(nthcon).name);
    cspec = load(tmpstr);
    bug = fast_gui_bug(cspec);
    if(bug)
      fprintf('\n\n');
      fprintf('ERROR: detected the FS-FAST GUI Bug in ');
      fprintf('contrast %s\n', clist(nthcon).name);
      fprintf('Please see https://surfer.nmr.mgh.harvard.edu/fswiki/FsFastGuiBug\n');
      fprintf('\n\n');
      ProcAnyway = getenv('FSF_PROC_GUI_BUG');
      if(isempty(ProcAnyway)) ProcAnyway = '0'; end
      ProcAnyway = sscanf(ProcAnyway,'%d');
      if(ProcAnyway ~= 1)
	fprintf('FSF_PROC_GUI_BUG not set to 1, so exiting\n');
	flac = [];
	return;
      end
      fprintf(' ... but FSF_PROC_GUI_BUG is set to 1, so continuing\n');
      fprintf('\n\n');
    end
    if(~isfield(cspec,'setwcond')) cspec.setwcond = 1; end
    if(~isfield(cspec,'sumconds')) cspec.sumconds = 1; end
    if(~isfield(cspec,'sumdelays')) cspec.sumdelays = 0; end
    if(~isfield(cspec,'setwdelay')) 
      if(ana.nregressors > 1) cspec.setwdelay = 1; 
      else                    cspec.setwdelay = 0; 
      end
    end
    cspec.name = clist(nthcon).name(1:end-4);
    if(~isfield(cspec,'CondState'))
      cspec.CondState = zeros(1,nconditions);
    end
    % This assures that weights are either +1, -1, or 0 when the
    % conditions are not summed. 4/9/09
    if(cspec.sumconds == 0) cspec.WCond = sign(cspec.WCond); end
    flac.ana.con(nthcon).cspec = cspec;
    flac.con(nthcon).name     = clist(nthcon).name(1:end-4);
    flac.con(nthcon).varsm    = 0;
    flac.con(nthcon).sumev    = cspec.sumconds;
    flac.con(nthcon).sumevreg = cspec.sumdelays;
    flac.con(nthcon).sumevrw  = [];
    flac.con(nthcon).rmprestim = cspec.RmPreStim; % 4/8/09
    flac.con(nthcon).ContrastMtx_0 = cspec.ContrastMtx_0; % 4/3/09
    con_nthev = 0;
    % Actual contrast matrix is computed below
    for nthcondition = 1:cspec.NCond
      if(cspec.WCond(nthcondition)==0) continue; end
      con_nthev = con_nthev + 1;
      flac.con(nthcon).ev(con_nthev).name = sprintf('Condition%02d',nthcondition);
      flac.con(nthcon).ev(con_nthev).evw  = cspec.WCond(nthcondition);
      flac.con(nthcon).ev(con_nthev).evrw = cspec.WDelay;
    end
  end % contrasts
  
  if(~isempty(extregList))
    if(isempty(nthcon)) nthcon = 0; end
    nlist = size(extregList,1);
    for n = 1:nlist
      extreg = deblank(extregList(n,:));
      nextreg = nextregList(n);
      extregevname = basename(extreg,'dat'); % remove .dat
      if(strcmp(extregevname,'mcextreg')) continue; end
      if(strcmp(extregevname,'mcprextreg')) continue; end
      nthcon = nthcon + 1;
      flac.con(nthcon).name     = extregevname;
      flac.con(nthcon).varsm    = 0;
      flac.con(nthcon).sumev    = 0;
      flac.con(nthcon).sumevreg = 0;
      flac.con(nthcon).sumevrw  = [];
      flac.con(nthcon).rmprestim = 0;
      flac.con(nthcon).ContrastMtx_0 = [];
      flac.con(nthcon).ev(1).name = extregevname;
      flac.con(nthcon).ev(1).evw  = 1;
      flac.con(nthcon).ev(1).evrw = ones(1,nextreg);
    end
  end
  
end

% Check each contrast
for nthcon = 1:ncontrasts

  % Make sure that each EV in the contrast is an EV in the FLAC
  for nthev = 1:length(flac.con(nthcon).ev)
    evindex = flac_evindex(flac,flac.con(nthcon).ev(nthev).name);
    if(isempty(evindex))
      fprintf('Contrast EV %s is not in the model\n',...
	      flac.con(nthcon).ev(nthev).name);
      flac = [];
      return;
    end
  end
  
  % Compute the contrast matrices
  % May want to defer this until customize to account for EVs with
  % a variable number of regressors
  flactmp = flac_conmat(flac,nthcon);
  if(isempty(flactmp))
    fprintf('ERROR: with contrast %s in anadir %s\n',...
	    flac.con(nthcon).name,anadir);
    flac = [];
    return;
  end
  flac = flactmp;
end

return;
