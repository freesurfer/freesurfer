function val = tdr_measasc(measasc,varname)
% val = tdr_measasc(measasc,varname)
%
% $Id: tdr_measasc.m,v 1.3 2005/03/19 00:18:57 greve Exp $
%
% EPI ------------------------------------------------------
% Number of echoes: sWiPMemBlock.alFree[2] (if it exists.
%  if it does not exist, then it's 1 echo).
% tDwell    = tdr_measasc(measasc,'sRXSPEC.alDwellTime[0]'); % nsec
% tRampUp   = tdr_measasc(measasc,'m_alRegridRampupTime');  % us
% tFlat     = tdr_measasc(measasc,'m_alRegridFlattopTime'); % us
% tRampDown = tdr_measasc(measasc,'m_alRegridRampdownTime'); % us
% tDelSamp  = tdr_measasc(measasc,'m_alRegridDelaySamplesTime'); % us
% echospacing = tdr_measasc(measasc,'m_lEchoSpacing'); % us
%   Should be the same as tRampUp+tFlat+tRampDown
% Can also get these from MDH
%   nlines (nrows) = sKSpace.lPhaseEncodingLines
%   nkcols = 'm_iNoOfFourierColumns = [N]';
%
% FID ---------------------------------------------------
%  Applies only after 11/15
% time to first echo (us) - alTE[0]  % us
% number of echoes - sWiPMemBlock.alFree[10]
% echo spacing     - sWiPMemBlock.adFree[2] % ms
%

val = [];
if(nargin ~= 2)
  fprintf('val = tdr_measasc(measasc,varname)\n');
  return;
end

fp = fopen(measasc,'r');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',measasc);
  return;
end

while(1)

  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fp);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  % Get count of number of items in the line %
  [items count] = sscanf(tline,'%s');

  % Read the keyword %  
  key = sscanf(tline,'%s',1);
  %fprintf('key = %s\n',key);

  if(strcmp(key,varname))
    valstr = sscanf(tline,'%*s %*s %s',1);
    ind = find(valstr ~= '[' & valstr ~= ']');
    valstr = valstr(ind);
    val = sscanf(valstr,'%f');
    %keyboard
    return;
  end

end % while (1)

fprintf('ERROR: could not fine %s in %s\n',varname,measasc);


return;