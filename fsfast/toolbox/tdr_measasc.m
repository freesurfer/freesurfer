function val = tdr_measasc(measasc,varname,flag)
% val = tdr_measasc(measasc,varname,<flag>)
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
% If flag = 'epi', reads in as an epipar struct
%  val.tDwell      % converted to usec
%  val.tRampUp     % us
%  val.tFlat       % us
%  val.tRampDown   % us
%  val.tDelSamp    % us
%  val.echospacing % us
%  val.TE          % us
%  val.nkcols      
%  val.nkrows    
%  val.timeunits = 'usec';



%
% tdr_measasc.m
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
  if(isempty(val.tDwell)) return; end
  val.tDwell = val.tDwell/1000; % convert to usec
  val.tRampUp   = tdr_measasc(measasc,'alRegridRampupTime[0]');  % us
  val.tFlat     = tdr_measasc(measasc,'alRegridFlattopTime[0]'); % us
  val.tRampDown = tdr_measasc(measasc,'alRegridRampdownTime[0]'); % us
  val.tDelSamp  = tdr_measasc(measasc,'alRegridDelaySamplesTime[0]'); % us

  val.echospacing  = tdr_measasc(measasc,'sFastImaging.lEchoSpacing'); % us
  val.TE        = tdr_measasc(measasc,'alTE[0]'); % us
  val.nkcols    =  tdr_measasc(measasc,'iNoOfFourierColumns');
  val.nkrows    =  tdr_measasc(measasc,'iNoOfFourierLines');
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
