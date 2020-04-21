function stem = flac_funcstem(flac,maskflag)
% stem = flac_funcstem(flac,maskflag)
% Determines the funcstem or maskstem. Should be
% consistent with the getana script.
% 

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
    if(~isempty(flac.sdf)) stc = '.sdf'; end
  else 
    stc = sprintf('.%s',flac.stc);
  end
  if(strcmp(flac.RawSpace,'native'))
    if(flac.rawfwhm == 0) fwhmstr = '';
    else fwhmstr = sprintf('.sm%g',flac.rawfwhm);
    end
    stem = sprintf('%s%s%s',flac.mcstem,stc,fwhmstr);
  elseif(strcmp(flac.RawSpace,'mni305'))
    stem = sprintf('%s%s.sm%g.mni305.%dmm%s',flac.mcstem,stc,flac.rawfwhm,flac.RawSpaceRes,ExpKey);
  elseif(strcmp(flac.RawSpace,'cvs_avg35_inMNI152'))
    stem = sprintf('%s%s.sm%g.cvs_avg35_inMNI152.%dmm%s',flac.mcstem,stc,flac.rawfwhm,flac.RawSpaceRes,ExpKey);
  else % surface
    if(isempty(flac.volsurffwhm))
      stem = sprintf('%s%s.sm%g.%s.%s%s',flac.mcstem,stc,flac.rawfwhm,flac.subject,flac.hemi,ExpKey);
    else
      % Smooth in the volume before sampling on the surface
      stem = sprintf('%s%s.vsm%g.sm%g.%s.%s%s',flac.mcstem,stc,flac.volsurffwhm,flac.rawfwhm,flac.subject,flac.hemi,ExpKey);
    end
  end
else
  % Return maskstem
  if(strcmp(flac.RawSpace,'native'))
    stem = sprintf('%s',flac.mask);
  elseif(strcmp(flac.RawSpace,'mni305'))
    stem = sprintf('%s.mni305.%dmm%s',flac.mask,flac.RawSpaceRes,ExpKey);
  elseif(strcmp(flac.RawSpace,'cvs_avg35_inMNI152'))
    stem = sprintf('%s.cvs_avg35_inMNI152.%dmm%s',flac.mask,flac.RawSpaceRes,ExpKey);
  else
    stem = sprintf('%s.%s.%s%s',flac.mask,flac.subject,flac.hemi,ExpKey);
  end
end  

return;


