function err = mkcontrast2m(configfile)
% err = mkcontrast2(configfile)
% Config file is as created by the mkcontrast2 script
% This matlab script was written to be compiled
% $Id: mkcontrast2m.m,v 1.1 2015/05/08 19:02:12 greve Exp $

fprintf('starting mkcontrast2m\n');

err = 1;
fid = fopen(configfile,'r');
if(fid == -1) 
  fprintf('ERROR: cannot open %s\n',configfile);
  return; 
end

monly = 0;
NCond = 0;
NDelay = 0;
WCond = [];
WDelay = [];
rdelta = [];
rtau = [];

while(1)

  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fid);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  % Get count of number of items in the line %
  [items count] = sscanf(tline,'%s');
  
  % Read the key %  
  key = sscanf(tline,'%s',1);
  tlinesplit = splitstring(tline);
  
  switch(key)
   case 'cmtxfile',  cmtxfile = sscanf(tline,'%*s %s',1);
   case 'monly',  monly = sscanf(tline,'%*s %d',1);
   case 'NCond',  NCond = sscanf(tline,'%*s %d',1);
   case 'WCond', 
    for n=1:(size(tlinesplit,1)-1)
      WCond(n) = sscanf(tlinesplit(n+1,:),'%f',1);  
    end
   case 'NDelay', NDelay = sscanf(tline,'%*s %d',1);
   case 'WDelay', 
    for n=1:(size(tlinesplit,1)-1)
      WDelay(n) = sscanf(tlinesplit(n+1,:),'%f',1);  
    end
   case 'TER',       TER = sscanf(tline,'%*s %f',1);    
   case 'sumconds',  sumconds = sscanf(tline,'%*s %d',1);    
   case 'sumdelays', sumdelays = sscanf(tline,'%*s %d',1);    
   case 'nircorr',   nircorr = sscanf(tline,'%*s %d',1);    
   case 'TPreStim',  TPreStim = sscanf(tline,'%*s %f',1);    
   case 'RmPreStim', RmPreStim = sscanf(tline,'%*s %d',1);    
   case 'rdelta', 
    for n=1:(size(tlinesplit,1)-1)
      rdelta(n) = sscanf(tlinesplit(n+1,:),'%f',1);  
    end
   case 'rtau', 
    for n=1:(size(tlinesplit,1)-1)
      rtau(n) = sscanf(tlinesplit(n+1,:),'%f',1);  
    end
   case 'ndelays',   ndelays = sscanf(tline,'%*s %d',1);    
   case 'CNorm',     CNorm = sscanf(tline,'%*s %d',1);    
   case 'setwdelay', setwdelay = sscanf(tline,'%*s %d',1);    
   case 'setwcond',  setwcond = sscanf(tline,'%*s %d',1);    
  end
end % while (1)
fclose(fid);

TimeWindow = NDelay*TER; % Hack so that it works with gammafit

if(isempty(WCond))
  WCond = ones(1,NCond);
end
if(isempty(WDelay))
  WDelay = ones(1,NDelay);
end
if(TimeWindow == 0) TimeWindow = NDelay*TER; end

if(nircorr == 0)
  % ---- Do not correlate with assumed response ----- %
  ContrastMtx_0 = fast_contrastmtx(TER,TimeWindow,TPreStim,...
				   NCond,sumconds,WCond,sumdelays,WDelay,RmPreStim,CNorm);
else
  % ---- Correlate with assumed response ------- %
  if(nircorr == 1)
    delta = (rdelta(2) + rdelta(1))/2;
    tau   = (rtau(2) + rtau(1))/2;
  else
    ddelta = (rdelta(2) - rdelta(1))/(nircorr-1);
    delta  = rdelta(1) + ddelta*[0:nircorr-1];
    dtau   = (rtau(2) - rtau(1))/(nircorr-1);
    tau    = rtau(1) + dtau*[0:nircorr-1];
  end
  
  t = -TPreStim + TER * [0:ndelays-1];
  h = fmri_hemodyn(t,delta,tau);
  
  ContrastMtx_0 = [];
  for n = 1 : nircorr
    wdelay = h(:,n)'; %'
    wdelay = wdelay/sum(wdelay);
    cm = fmri_mrestriction2([WCond],wdelay,sumconds,sumdelays);
    ContrastMtx_0 = [  ContrastMtx_0; cm ];
  end
  
end

save(cmtxfile,'ContrastMtx_0','NCond','WCond','WDelay','CNorm','TER',...
     'TimeWindow','TPreStim','RmPreStim','sumconds','sumdelays',...
     'nircorr','rdelta','rtau','setwdelay','setwcond','-V4');
fprintf('mkcontrast2m done\n');

return;
  
