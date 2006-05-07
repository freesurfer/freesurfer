function val = tdr_measasc(measasc,varname,flag)
% val = tdr_measasc(measasc,varname,<flag>)
%
% $Id: tdr_measasc.m,v 1.4 2006/05/07 23:18:13 greve Exp $
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
% Can also val = tdr_measasc(measasc,[],'epi'). Then val
%   is an epipar structure.
%
% FID ---------------------------------------------------
%  Applies only after 11/15
% time to first echo (us) - alTE[0]  % us
% number of echoes - sWiPMemBlock.alFree[10]
% echo spacing     - sWiPMemBlock.adFree[2] % ms
%

val = [];
if(nargin < 2 | nargin > 3)
  fprintf('val = tdr_measasc(measasc,varname,<flag>)\n');
  return;
end

if(~exist('flag','var')) flag = []; end
if(~isempty(flag))
  if(~strcmp(flag,'epi'))
    fprintf('ERROR: flag (%s) must be epi or empty\n',flag);
    return;
  end
  % Read in as epipar struct.
  val.tDwell = tdr_measasc(measasc,'sRXSPEC.alDwellTime[0]'); % nsec
  val.tDwell = val.tDwell/1000; % convert to usec
  val.tRampUp   = tdr_measasc(measasc,'m_alRegridRampupTime');  % us
  val.tFlat     = tdr_measasc(measasc,'m_alRegridFlattopTime'); % us
  val.tRampDown = tdr_measasc(measasc,'m_alRegridRampdownTime'); % us
  val.tDelSamp  = tdr_measasc(measasc,'m_alRegridDelaySamplesTime'); % us
  val.echospacing = tdr_measasc(measasc,'m_lEchoSpacing'); % us
  val.TE        = tdr_measasc(measasc,'alTE[0]'); % us
  val.nkcols    =  tdr_measasc(measasc,'m_iNoOfFourierColumns');
  val.nkrows    =  tdr_measasc(measasc,'m_iNoOfFourierLines');
  val.timeunits = 'usec';
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
