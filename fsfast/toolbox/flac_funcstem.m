function stem = flac_funcstem(flac,maskflag)
% stem = flac_funcstem(flac,maskflag)
% Determines the funcstem or maskstem. Should be
% consistent with the getana script.
% 
% $Id: flac_funcstem.m,v 1.2 2011/03/11 23:48:16 greve Exp $

stem = [];
if(nargin < 1 | nargin > 2)
  fprintf('stem = flac_funcstem(flac,maskflag)\n');
  return;
end
if(~exist('maskflag','var')) maskflag = 0; end

ExpKey = '';
if(~isempty(flac.ExpKey)) ExpKey = sprintf('.%s',flac.ExpKey); end

if(~maskflag)
  % Return funcstem
  if(strcmp(flac.stc,'none')) stc = ''; 
  else stc = sprintf('.%s',flac.stc);
  end
  if(strcmp(flac.RawSpace,'native'))
    if(flac.rawfwhm == 0) fwhmstr = '';
    else fwhmstr = sprintf('.sm%g',flac.rawfwhm);
    end
    stem = sprintf('%s%s%s',flac.mcstem,stc,fwhmstr);
  elseif(strcmp(flac.RawSpace,'mni305'))
    stem = sprintf('%s%s.sm%g.mni305.%dmm%s',flac.mcstem,stc,flac.rawfwhm,flac.RawSpaceRes,ExpKey);
  else % surface
    stem = sprintf('%s%s.sm%g.%s.%s%s',flac.mcstem,stc,flac.rawfwhm,flac.subject,flac.hemi,ExpKey);
  end
else
  % Return maskstem
  if(strcmp(flac.RawSpace,'native'))
    stem = sprintf('%s',flac.mask);
  elseif(strcmp(flac.RawSpace,'mni305'))
    stem = sprintf('%s.mni305.%dmm%s',flac.mask,flac.RawSpaceRes,ExpKey);
  else
    stem = sprintf('%s.%s.%s%s',flac.mask,flac.subject,flac.hemi,ExpKey);
  end
end  





return;


