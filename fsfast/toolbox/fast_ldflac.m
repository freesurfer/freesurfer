function flac = fast_ldflac(flacfile)
% flac = fast_ldflac(flacfile)
%
% $Id: fast_ldflac.m,v 1.1 2004/10/16 05:16:07 greve Exp $

flac = [];
if(nargin > 1)
  fprintf('flac = fast_ldflac(flacfile)\n');
  return;
end

flac.name = '';
flac.fsd = '';
flac.funcstem = '';
flac.runlistfile = '';
flac.TR = [];
flac.mask = '';
flac.inorm = [];
flac.whiten = 0;
flac.fixacf = 0;
%flac.ev = []; % Leave commented
if(nargin == 0) return; end

fp = fopen(flacfile,'r');
if(fp == -1)
  fpritf('ERROR: could not open %s\n',flacfile);
  flac = [];
  return;
end

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
   case 'flacname',    flac.name        = sscanf(tline,'%*s %s',1);
   case 'fsd',         flac.fsd         = sscanf(tline,'%*s %s',1);
   case 'TR',          flac.TR          = sscanf(tline,'%*s %f',1);
   case 'funcstem',    flac.funcstem    = sscanf(tline,'%*s %s',1);
   case 'mask',        flac.mask        = sscanf(tline,'%*s %s',1);
   case 'inorm',       flac.inorm       = sscanf(tline,'%*s %f',1);
   case 'runlistfile', flac.runlistfile = sscanf(tline,'%*s %s',1);
   case 'whiten',      flac.whiten      = sscanf(tline,'%*s %d',1);
   case 'acfbins',     flac.acfbins     = sscanf(tline,'%*s %d',1);
   case 'fixacf',      flac.fixacf      = sscanf(tline,'%*s %d',1);
   case 'EV', 
    flac.ev(nthev) = fast_ev_parse(tline);
    nthev = nthev + 1;
   otherwise
    fprintf('INFO: key %s unrecognized, skipping\n',key);
  end
end % while (1)

fclose(fp);

return;
