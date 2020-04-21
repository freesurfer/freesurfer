function [err,msg] = fast_svana(ananame,flac)
% [err,msg] = fast_svana(ananame,ana)

err = 1;
if(nargin ~= 2)
  fprintf('[err,msg] = fast_svana(ananame,ana)\n');
  return;
end

err = mkdirp(ananame);
if(err)
  msg = sprintf('ERROR: creating dir %s in %s',ananame,pwd);
  fprtinf('%s\n',msg);
  return;
end

flac.ana.inorm = flac.inorm;
flac.ana.delay = flac.stimulusdelay;
ana = flac.ana;

if(strcmp(ana.designtype,'event-related') | ...
   strcmp(ana.designtype,'blocked')) IsERBlock = 1;
else                                 IsERBlock = 0;
end

% Save the full flac.
flac.CreationDate = date;
flacmat = sprintf('%s/fsfast.flac',ananame);
save(flacmat,'flac');

% Create creator.txt
tmpfile = sprintf('%s/creator.txt',ananame);
fp = fopen(tmpfile,'w');
fprintf(fp,'%s',flac.creator);
fclose(fp);

%------------- analysis.cfg ---------------------
anacfg = sprintf('%s/analysis.cfg',ananame);
[fp msg] = fopen(anacfg,'w');
if(fp == -1)
  fprtinf('%s\n',msg);
  return;
end

if(IsERBlock)
  if(ana.gammafit)
    fprintf(fp,'-gammafit %f %f\n',ana.gamdelay,ana.gamtau);
    fprintf(fp,'-gammaexp %f\n',ana.gamexp);
  end
  if(ana.spmhrffit)
    fprintf(fp,'-spmhrf %d\n',ana.nspmhrfderiv);
  end
  fprintf(fp,'-timewindow %f\n',ana.timewindow);
  fprintf(fp,'-prestim %f\n',ana.prestim);
  fprintf(fp,'-acfbins %d\n',flac.acfbins);
  if(flac.fixacf)
    fprintf(fp,'-fix-acf\n');
  end
  if(flac.autostimdur) fprintf(fp,'-autostimdur\n');
  else fprintf(fp,'-noautostimdur\n');
  end  
  fprintf(fp,'-TER %f\n',ana.TER);
  fprintf(fp,'-polyfit %d\n',ana.PolyOrder);
else
  if(~isempty(ana.ncycles))
    fprintf(fp,'-ncycles %d\n',ana.ncycles);
  end  
end

if(ana.inorm)
  fprintf(fp,'-rescale %g\n',ana.inorm);
end  
  
if(~isempty(ana.extreg))
  fprintf(fp,'-extreg %s\n',ana.extreg);
  fprintf(fp,'-nextreg %d\n',ana.nextreg);
end  
if(~isempty(ana.delay))
  fprintf(fp,'-delay %f\n',ana.delay);
end  

fclose(fp);

%------------- analysis.info ---------------------
anainfo = sprintf('%s/analysis.info',ananame);
[fp msg] = fopen(anainfo,'w');
if(fp == -1)
  fprtinf('%s\n',msg);
  return;
end

fprintf(fp,'analysis %s\n',ana.analysis);
fprintf(fp,'TR %f\n',flac.TR);
fprintf(fp,'RefEventDur %f\n',flac.RefEventDur);
fprintf(fp,'fsd %s\n',flac.fsd);
fprintf(fp,'funcstem %s\n',flac.funcstem);
fprintf(fp,'runlistfile %s\n',flac.runlistfile);
fprintf(fp,'tpexclude %s\n',flac.tpexcfile);
if(~isempty(flac.parfile))
  fprintf(fp,'parname %s\n',flac.parfile);
end
fprintf(fp,'designtype %s\n',ana.designtype);
if(~isempty(ana.nconditions))
  fprintf(fp,'nconditions %d\n',ana.nconditions);
end

if(~isempty(flac.mask))
  fprintf(fp,'maskstem %s\n',flac.mask);
end
if(~isempty(flac.subject))
  fprintf(fp,'surface %s %s\n',flac.subject,flac.hemi);
end
if(flac.UseTalairach)
  fprintf(fp,'UseTalairach\n');
end

if(IsERBlock)
  % Write out the condition names
  if(~isfield(ana,'ConditionNames'))ana.ConditionNames = '';end
  if(isempty(ana.ConditionNames))
    for n = 1:ana.nconditions
      tmp = sprintf('Condition%02d',n);
      ana.ConditionNames = strvcat(ana.ConditionNames,tmp);
    end
  end
  for n = 1:ana.nconditions
    fprintf(fp,'Condition %d %s\n',n,deblank(ana.ConditionNames(n,:)));
  end
end
fclose(fp);

%---------------------------------------------------------
delete_old_contrasts = 1;
if(delete_old_contrasts)
  tmp = sprintf('%s/*.mat',ananame);
  d = dir(tmp);
  if(~isempty(d))
    for nthcon = 1:length(d)
      fname = sprintf('%s/%s',ananame,d(nthcon).name);
      fprintf('Deleting old contrast %s\n',fname);
      delete(fname);
    end
  end
end

if(~IsERBlock) return; end

%---------------------------------------------------------
ncon = length(ana.con);
for nthcon = 1:ncon
  cmtxfile = sprintf('%s/%s.mat',ananame,flac.ana.con(nthcon).cspec.name);
  CNorm = ana.con(nthcon).cspec.CNorm;
  ContrastMtx_0 = ana.con(nthcon).cspec.ContrastMtx_0;
  NCond = ana.con(nthcon).cspec.NCond;
  WCond = ana.con(nthcon).cspec.WCond;
  % This line fixes the FSFAST GUI bug
  if(CNorm) WCond = fast_norm_con(WCond); end
  WDelay = ana.con(nthcon).cspec.WDelay;
  TER = ana.con(nthcon).cspec.TER;
  TimeWindow = ana.con(nthcon).cspec.TimeWindow;
  TPreStim = ana.con(nthcon).cspec.TPreStim;
  RmPreStim = ana.con(nthcon).cspec.RmPreStim;
  sumconds = ana.con(nthcon).cspec.sumconds;
  sumdelays = ana.con(nthcon).cspec.sumdelays;
  nircorr = ana.con(nthcon).cspec.nircorr;
  rdelta = ana.con(nthcon).cspec.rdelta;
  rtau = ana.con(nthcon).cspec.rtau;
  setwcond  = ana.con(nthcon).cspec.setwcond;
  setwdelay = ana.con(nthcon).cspec.setwdelay;
  CondState = ana.con(nthcon).cspec.CondState;
  save(cmtxfile,'ContrastMtx_0','NCond','WCond','WDelay','CNorm','TER',...
    'TimeWindow','TPreStim','RmPreStim','sumconds','sumdelays',...
    'nircorr','rdelta','rtau','setwdelay','setwcond','CondState','-V4');
end

err = 0;
return;
