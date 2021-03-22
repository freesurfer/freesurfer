function flac = fast_ldanaflac(anadir)
% flac = fast_ldanaflac(anadir)
%
%


%
% fast_ldanaflac.m
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
nconditions = [];
ConditionNames = '';
PolyOrder = 0;
period = [];
gammafit = 0;
gamdelay = 0;
gamtau = 0;
gamexp = 2;
ngamderiv = 0;
spmhrffit = 0;
nspmhrfderiv = 0;
timewindow = 0;
prestim = 0;
TER = [];
taskregList = [];
ntaskregList = [];
nuisregList = [];
nnuisregList = [];
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
    if(~isempty(tline)) 
      if(tline(1) ~= '#') break; end; 
    end
  end
  if(tline(1) == -1) break; end

  key = sscanf(tline,'%s',1);
  %fprintf('key = %s\n',key);
  
  switch(key)
   case 'analysis',    analysis         = sscanf(tline,'%*s %s',1);
   case 'funcstem',    flac.funcstem    = sscanf(tline,'%*s %s',1);
   case 'mcstem',      flac.mcstem      = sscanf(tline,'%*s %s',1);
   case 'mask',        flac.mask        = sscanf(tline,'%*s %s',1);
   case 'fsd',         flac.fsd         = sscanf(tline,'%*s %s',1);
   case 'runlistfile', flac.runlistfile = sscanf(tline,'%*s %s',1);
   case 'inorm',       flac.inorm       = sscanf(tline,'%*s %f',1);
   case 'TR',          flac.TR          = sscanf(tline,'%*s %f',1);
   case 'OverrideTR',  flac.OverrideTR  = sscanf(tline,'%*s %d',1);
   case 'RegDOF',      flac.RegDOF      = sscanf(tline,'%*s %d',1);
   case 'RawFWHM',     flac.rawfwhm     = sscanf(tline,'%*s %f',1);
   case 'VolSurfFWHM', flac.volsurffwhm = sscanf(tline,'%*s %f',1);
   case 'RawSTC',      flac.stc         = sscanf(tline,'%*s %s',1);
   case 'SliceDelayFile', flac.sdf         = sscanf(tline,'%*s %s',1);
   case 'acfbins',     flac.acfbins     = sscanf(tline,'%*s %d',1);
   case 'acffwhm',     flac.acffwhm     = sscanf(tline,'%*s %f',1);
   case 'acfsvd',      flac.acfsvd      = sscanf(tline,'%*s %d',1);
   case 'fixacf',      flac.fixacf      = sscanf(tline,'%*s %d',1);
   case 'fsv3-whiten', flac.fsv3_whiten = sscanf(tline,'%*s %d',1);
   case 'HPFCutoffHz', flac.hpfCutoffHz = sscanf(tline,'%*s %f',1);
   case 'HeteroGCor',  flac.HeteroGCor  = sscanf(tline,'%*s %f',1);
   case 'polyfit',     PolyOrder        = sscanf(tline,'%*s %d',1);
   case 'tpexclude',   flac.tpexcfile   = sscanf(tline,'%*s %s',1);
   case 'nskip',       flac.nskip       = sscanf(tline,'%*s %d',1);
   case 'parname',     flac.parfile     = sscanf(tline,'%*s %s',1);
   case 'designtype',  designtype       = sscanf(tline,'%*s %s',1);
   case 'RefEventDur', flac.RefEventDur = sscanf(tline,'%*s %f',1);
   case 'prestim',     prestim          = sscanf(tline,'%*s %f',1);
   case 'timewindow',  timewindow       = sscanf(tline,'%*s %f',1);
   case 'TER',         TER              = sscanf(tline,'%*s %f',1);
   case 'stimulusdelay', flac.stimulusdelay = sscanf(tline,'%*s %f',1);
   case 'nconditions', nconditions      = sscanf(tline,'%*s %d',1);
   case 'TFILTER', 
    tfilter = flac_tfilter_parse(tline);
    if(isempty(tfilter)) 
      flac=[]; 
      fprintf('line %d\n',nthline);
      return; 
    end
    flac.tfilter = tfilter;
   case 'gamma', 
    gammafit = 1;
    gamdelay = sscanf(tline,'%*s %f',1);
    gamtau   = sscanf(tline,'%*s %*f %f',1);
    gamexp   = sscanf(tline,'%*s %*f %*f %f',1);
    ngamderiv = sscanf(tline,'%*s %*f %*f %*f %d',1);
   case 'spmhrf', 
    spmhrffit = 1; 
    nspmhrfderiv = sscanf(tline,'%*s %d',1);
   case 'period', period  = sscanf(tline,'%*s %f',1);
    
   case 'RawSpace',     
    flac.RawSpaceType = sscanf(tline,'%*s %s',1);
    flac.RawSpace = sscanf(tline,'%*s %*s %s',1);
    if(strcmp(flac.RawSpace,'mni305'))
      flac.RawSpaceRes = sscanf(tline,'%*s %*s %*s %f',1); 
    end
    if(strcmp(flac.RawSpace,'cvs_avg35_inMNI152'))
      flac.RawSpaceRes = sscanf(tline,'%*s %*s %*s %f',1); 
    end
    if(strcmp(flac.RawSpaceType,'surface'))
      flac.subject = sscanf(tline,'%*s %*s %s',1);
      flac.hemi = sscanf(tline,'%*s %*s %*s %s',1);
    end
   %case 'UseTalairach',flac.UseTalairach = 1;
   case 'ExpKey', flac.ExpKey = sscanf(tline,'%*s %s',1);
   case 'Condition', 
    conditionid   = sscanf(tline,'%*s %d',1);
    ConditionName = sscanfitem(tline,3);
    ConditionNames = strvcat(ConditionNames,ConditionName);
   case 'taskreg'
    taskreg = sscanf(tline,'%*s %s',1);
    ntaskreg = sscanf(tline,'%*s %*s %d',1);
    taskregList = strvcat(taskregList,taskreg);
    ntaskregList = [ntaskregList ntaskreg];
   case 'nuisreg'
    nuisreg = sscanf(tline,'%*s %s',1);
    nnuisreg = sscanf(tline,'%*s %*s %d',1);
    nuisregList = strvcat(nuisregList,nuisreg);
    nnuisregList = [nnuisregList nnuisreg];
   case 'UseB0DC'
    UseB0DC = sscanf(tline,'%*s %d',1);
   case 'ApplySubCortMask'
    ApplySubCortMask = sscanf(tline,'%*s %d',1);
   case 'PerSession'
    flac.PerSession = sscanf(tline,'%*s %d',1);
   otherwise
    fprintf('INFO: key %s unrecognized, line %d, skipping\n',key,nthline);
  end
  nthline = nthline + 1;
end % while (1)
fclose(fp);
if(isempty(flac.funcstem))
  flac.funcstem = flac_funcstem(flac,0);
end

if(isempty(flac.PerSession))
  if(strcmp(flac.RawSpace,'native'))
    flac.PerSession = 1;
  else
    flac.PerSession = 0;
  end
end

%fprintf('RED %g\n',flac.RefEventDur);

if(isempty(flac.mcstem)) 
  fprintf('ERROR: no mcstem specified in %s\n',infofile);
  flac = [];
  return;
end

if(isempty(flac.mask)) 
  flac.mask = 'brain';
  fprintf('INFO: mask is not set, setting to brain\n');
end
if(strcmp(flac.mask,'nomask')) flac.mask = ''; end

if(isempty(flac.fsd)) flac.fsd = 'bold'; end 
if(isempty(flac.acfsegstem)) flac.acfsegstem = 'acfseg'; end 

flac.designtype = designtype;

ana.analysis     = analysis;
ana.info = char(textread(info,'%s','delimiter','\n'));
ana.designtype   = designtype;
ana.PolyOrder    = PolyOrder;
%ana.extregList   = extregList;
%ana.nextregList  = nextregList;
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
ana.ngamderiv    = ngamderiv;

ana.spmhrffit    = spmhrffit;
ana.nspmhrfderiv = nspmhrfderiv;

if(strcmp(designtype,'event-related') | strcmp(designtype,'blocked'))
if(ana.gammafit) 
  ana.nregressors = 1; 
end
if(ana.spmhrffit) 
  ana.nregressors = ana.nspmhrfderiv + 1;
end
if(ana.firfit)
  ana.nregressors = round(ana.timewindow/ana.TER);
end
end

flac.ana = ana;
%flac.ana.con = [];

nthev = 1;
tline = sprintf('EV Baseline baseline nuis');
flac.ev(nthev) = flac_ev_parse(tline);
nthev = nthev+1;

if(strcmp(designtype,'event-related') | strcmp(designtype,'blocked'))
  for n=1:nconditions
    if(gammafit)
      tline = sprintf('EV Condition%02d gamma task %d %g %g %g %d %g',...
		      n,n,gamdelay,gamtau,gamexp,ngamderiv,TER);
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

nthcon = 0;
if(strcmp(designtype,'retinotopy')) 
  flac.IsRetinotopy = 1; 
  % par file will have:
  %   stimtype  eccen (retinotopy)
  %   direction neg
  % rtopy output will go in ananame/{eccen,polar}
  % Need to add nuis on either side of fund and harm to be
  % compatible with selfreqavg
  if(isempty(period))
    fprintf('ERROR: Must specify -period in mkanalysis-sess');
    flac = [];
    return;
  end

  tline = sprintf('EV eccen abret task %g',period);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;
  tline = sprintf('EV polar abret task %g',period);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;

  for nthstimtype = 1:2
    nthcon = nthcon + 1;
    if(nthstimtype == 1) stimtype = 'eccen'; end
    if(nthstimtype == 2) stimtype = 'polar'; end
    flac.con(nthcon).name     = stimtype;
    flac.con(nthcon).varsm    = 0;
    flac.con(nthcon).sumev    = 0;
    flac.con(nthcon).sumevreg = 0;
    flac.con(nthcon).sumevrw  = [];
    flac.con(nthcon).evrw  = [];
    flac.con(nthcon).ev(1).name = stimtype;
    flac.con(nthcon).ev(1).evw  = 1;
    flac.con(nthcon).ev(1).evrw = [1 1 zeros(1,10)];
    flac.con(nthcon).rmprestim = 0;
    flac.con(nthcon).cspec.name = flac.con(nthcon).name;
    flac.con(nthcon).UseExtC = 0;
    flac.ana.con(nthcon) = flac.con(nthcon);
  end
  ncontrasts = length(flac.con);
end

if(strcmp(designtype,'abblocked'))
  flac.IsABBlocked = 1;
  nharmonics = 1;
  tline = sprintf('EV Fourier abret task %g',period);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'omnibus';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).evrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [1 1 1 1 zeros(1,8)];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.con(nthcon).UseExtC = 0;
  flac.ana.con = flac.con(nthcon);

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'fund';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).evrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [1 1 0 0 zeros(1,8)];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.con(nthcon).UseExtC = 0;
  flac.ana.con(nthcon) = flac.con(nthcon);

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'fund-sin';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).evrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [1 0 0 0 zeros(1,8)];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.con(nthcon).UseExtC = 0;
  flac.ana.con(nthcon) = flac.con(nthcon);

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'fund-cos';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).evrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [0 1 0 0 zeros(1,8)];
    flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.con(nthcon).UseExtC = 0;
  flac.ana.con(nthcon) = flac.con(nthcon);

  nthcon = nthcon + 1;
  flac.con(nthcon).name     = 'harm';
  flac.con(nthcon).varsm    = 0;
  flac.con(nthcon).sumev    = 0;
  flac.con(nthcon).sumevreg = 0;
  flac.con(nthcon).sumevrw  = [];
  flac.con(nthcon).evrw  = [];
  flac.con(nthcon).ev(1).name = 'Fourier';
  flac.con(nthcon).ev(1).evw  = 1;
  flac.con(nthcon).ev(1).evrw = [0 0 1 1 zeros(1,8)];
  flac.con(nthcon).rmprestim = 0;
  flac.con(nthcon).cspec.name = flac.con(nthcon).name;
  flac.con(nthcon).UseExtC = 0;
  flac.ana.con(nthcon) = flac.con(nthcon);
  
  ncontrasts = length(flac.con);
end

% Polynomial and highpass filtering
if(PolyOrder > 0 & isempty(flac.hpfCutoffHz))
  tline = sprintf('EV Poly polynomial nuis %d',PolyOrder);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;
end
if(PolyOrder == 0 & ~isempty(flac.hpfCutoffHz))
  tline = sprintf('EV HPF hpf nuis %f',flac.hpfCutoffHz);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;
end
if(PolyOrder > 0 & ~isempty(flac.hpfCutoffHz))
  tline = sprintf('EV HPFPoly hpf+poly nuis %f %d',flac.hpfCutoffHz,PolyOrder);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;
end

% Non-paramatric (external) Task regressors 
if(~isempty(taskregList))
  nlist = size(taskregList,1);
  for n = 1:nlist
    taskreg = deblank(taskregList(n,:));
    ntaskreg = ntaskregList(n);
    taskregevname = basename(taskreg,'dat'); % remove .dat
    tline = sprintf('EV %s nonpar task %s %d',taskregevname,taskreg,ntaskreg);
    flac.ev(nthev) = flac_ev_parse(tline);
    nthev = nthev+1;
  end
end

% Non-paramatric (external) Nuissance regressors 
if(~isempty(nuisregList))
  nlist = size(nuisregList,1);
  for n = 1:nlist
    nuisreg = deblank(nuisregList(n,:));
    nnuisreg = nnuisregList(n);
    nuisregevname = basename(nuisreg,'dat'); % remove .dat
    tline = sprintf('EV %s nonpar nuis %s %d',nuisregevname,nuisreg,nnuisreg);
    flac.ev(nthev) = flac_ev_parse(tline);
    nthev = nthev+1;
  end
end

% Skipping and time-point exclude file
if(~isempty(flac.tpexcfile) | flac.nskip > 0)
  if(flac.nskip > 0 & isempty(flac.tpexcfile)) tpexcfile = 'junk.nskip.junk';
  else tpexcfile = flac.tpexcfile;
  end
  tline = sprintf('EV TPExclude texclude nuis %s',tpexcfile);
  flac.ev(nthev) = flac_ev_parse(tline);
  nthev = nthev+1;
end

ncontrasts = 0;
if(strcmp(designtype,'event-related') | strcmp(designtype,'blocked'))
  %-------------- contrasts --------------------------
  tmpstr = sprintf('%s/*.mat',anadir);
  clist = dir(tmpstr);
  ncontrasts = length(clist);
  for nthc = 1:ncontrasts
    nthcon = nthcon+1;
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
    flac.con(nthcon).evrw  = [];
    flac.con(nthcon).rmprestim = cspec.RmPreStim; % 4/8/09
    flac.con(nthcon).ContrastMtx_0 = cspec.ContrastMtx_0; % 4/3/09
    flac.con(nthcon).UseExtC = 0;
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

  % Get "external" contrasts
  tmpstr = sprintf('%s/*.mtx',anadir);
  clist = dir(tmpstr);
  ncontrasts2 = length(clist);
  for nthc = 1:ncontrasts2
    nthcon = nthcon+1;
    cname = clist(nthc).name;
    fprintf('%2d %s\n',nthcon,cname);
    tmpstr = sprintf('%s/%s',anadir,cname);
    cspec = [];
    cspec.name = cname(1:length(cname)-4);
    cspec.ContrastMtx_0 = load(tmpstr);
    flac.con(nthcon).name = cspec.name;
    flac.con(nthcon).cspec = cspec;
    flac.con(nthcon).ContrastMtx_0 = cspec.ContrastMtx_0;
    flac.con(nthcon).UseExtC = 1;
  end
end

% Contrast for task regressors
if(~isempty(taskregList))
  nlist = size(taskregList,1);
  for n = 1:nlist
    taskreg = deblank(taskregList(n,:));
    ntaskreg = ntaskregList(n);
    taskregevname = basename(taskreg,'dat'); % remove .dat
    nthcon = nthcon + 1;
    flac.con(nthcon).name     = taskregevname;
    flac.con(nthcon).varsm    = 0;
    flac.con(nthcon).sumev    = 0;
    flac.con(nthcon).sumevreg = 0;
    flac.con(nthcon).sumevrw  = [];
    flac.con(nthcon).evrw  = [];
    flac.con(nthcon).rmprestim = 0;
    flac.con(nthcon).ContrastMtx_0 = [];
    flac.con(nthcon).ev(1).name = taskregevname;
    flac.con(nthcon).ev(1).evw  = 1;
    flac.con(nthcon).UseExtC = 0;
    if(ntaskreg > 0)
      flac.con(nthcon).ev(1).evrw = ones(1,ntaskreg);
    else
      flac.con(nthcon).ev(1).evrw = [];
    end
  end
end
ncontrasts = nthcon;

if(ncontrasts == 0) flac.ana.con = []; end

% Check each contrast
for nthcon = 1:ncontrasts
  con = flac.con(nthcon);
  if(con.UseExtC)
    ntaskevs = length(flac_evtaskind(flac));
    if(ntaskevs ~= size(con.ContrastMtx_0,2))
      fprintf('ERROR: Contrast %s. Number of columns %d != number of task evs %d\n',...
	      con.name, size(con.ContrastMtx_0,2), ntaskevs);
      flac = [];
      return;
    end
    continue;
  end
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
  if(0)
  flactmp = flac_conmat(flac,nthcon);
  if(isempty(flactmp))
    fprintf('ERROR: with contrast %s in anadir %s\n',...
	    flac.con(nthcon).name,anadir);
    flac = [];
    return;
  end
  flac = flactmp;
  end
end

return;

