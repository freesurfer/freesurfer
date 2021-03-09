function flacfg = fast_flacfg_load(cfgfile)
% flacfg = fast_flacfg_load(cfgfile)


%
% fast_flacfg_load.m
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


flacfg = [];

if(nargin ~= 1)
  fprintf('flacfg = fast_flacfg_load(cfgfile)\n');
  return;
end

fid = fopen(cfgfile,'r');
if(fid == -1)
  fprintf('ERROR: cannot open cfg file %s\n',cfgfile);
  return;
end

tline = fgetl(fid);

flacfgid = sscanf(tline,'%s',1);
if(~strcmp(flacfgid,'FSFAST-FLACFG'))
  fprintf('ERROR: format error (id string) in cfg file %s\n',cfgfile);
  fclose(fid);
  return;
end

flacfgversion = sscanf(tline,'%*s %d',1);
if(flacfgversion ~= 2)
  fprintf('ERROR: cfg file %s, version = %d\n',cfgfile,flacfgversion);
  fclose(fid);
  return;
end

flacfg = fast_flacfg_struct;
nthfx = 1;  
while(1)

  % scroll through any blank lines or comments
  while(1)
    tline = fgetl(fid);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  key = lower(sscanf(tline,'%s',1));
  %fprintf('key = %s\n',key);
  
  switch(key)
   case 'flaname',     flacfg.flaname     = sscanf(tline,'%*s %s',1);
   case 'nskip',       flacfg.nskip       = sscanf(tline,'%*s %d',1);
   case 'tr',          flacfg.TR          = sscanf(tline,'%*s %f',1);
   case 'fsd',         flacfg.fsd         = sscanf(tline,'%*s %s',1);
   case 'funcstem',    flacfg.funcstem    = sscanf(tline,'%*s %s',1);
   case 'maskstem',    flacfg.maskstem    = sscanf(tline,'%*s %s',1);
   case 'inorm',       flacfg.inorm       = sscanf(tline,'%*s %f',1);
   case 'runlistfile', flacfg.runlistfile = sscanf(tline,'%*s %s',1);
   case 'evschfname',  flacfg.evschfname  = sscanf(tline,'%*s %s',1);
   case 'slicetiming', flacfg.slicetiming = sscanf(tline,'%*s %s',1);
   case 'usetpexclude', flacfg.usetpexclude = sscanf(tline,'%*s %d',1);
   case 'effect', 
    flacfg.fxlist(nthfx).fx = fast_fxcfg('parseline',tline);
    nthfx = nthfx + 1;
   otherwise
    fprintf('INFO: key %s unrecognized, skipping\n',key);
  end
end % while (1)

fclose(fid);

flacfg.sesscfg = fast_sesscfg_struct;


return;
