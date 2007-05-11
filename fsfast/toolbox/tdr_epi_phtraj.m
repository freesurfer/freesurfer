function [phx, phy, epipar] = tdr_epi_phtraj(arg1,delay)
% [phx, phy, epipar] = tdr_epi_phtraj(arg1,<delay>)
% 
% Assumes phase space trajectory starts in upper right (phx>0 and
% phy>0), then goes to upper left on the first line (phx>0 and
% phy<0), then drops down one line, then goes to the right, etc.
%
% If delay is specified, then the waveform delayed by delay usec.
%
% $Id: tdr_epi_phtraj.m,v 1.4 2007/05/11 20:52:07 greve Exp $

%
% tdr_epi_phtraj.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/05/11 20:52:07 $
%    $Revision: 1.4 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

phx = [];
phy = [];
epipar = [];
if(nargin ~= 1 & nargin ~= 2)
  fprintf('[phx phy epipar] = tdr_epi_phtraj(measasc)\n');
  fprintf('[phx phy epipar] = tdr_epi_phtraj(epipar)\n');
  return;
end

% The 1st arg is either the path to a meas.asc file or an epipar struct
if(isstr(arg1))
  measasc = arg1;
  epipar.tDwell = tdr_measasc(measasc,'sRXSPEC.alDwellTime[0]'); % nsec
  if(isempty(epipar.tDwell)) return; end
  epipar.tDwell = epipar.tDwell/1000; % convert to usec
  epipar.tRampUp   = tdr_measasc(measasc,'alRegridRampupTime[0]');  % us
  epipar.tFlat     = tdr_measasc(measasc,'alRegridFlattopTime[0]'); % us
  epipar.tRampDown = tdr_measasc(measasc,'alRegridRampdownTime[0]'); % us
  epipar.tDelSamp  = tdr_measasc(measasc,'alRegridDelaySamplesTime[0]'); % us
  epipar.echospacing = tdr_measasc(measasc,'sFastImaging.lEchoSpacing'); % us
  epipar.nkcols    = tdr_measasc(measasc,'iNoOfFourierColumns'); 
  epipar.nkrows    = tdr_measasc(measasc,'iNoOfFourierLines'); 
else
  epipar = arg1;
end

if(~exist('delay','var')) delay = []; end
if(isempty(delay)) delay = 0; end

% This is the phase trajectory for a column. It starts at -pi and
% goes to pi-delta, passing thru 0 between (nkcols/2) and
% (nkcols/2)+1. 
phx0 = kspacevector2(epipar.nkcols,epipar.tDwell,epipar.tRampUp,...
		     epipar.tFlat,epipar.tRampDown,epipar.tDelSamp,delay);

phx = repmat(phx0,[epipar.nkrows 1]);
% Assume first phx line starts pos and goes neg
phx(1:2:end,:) = fliplr(phx(1:2:end,:)); 

tmp = [0:epipar.nkrows-1]'/epipar.nkrows;
tmp = tmp - tmp(round(epipar.nkrows/2));
% Assume first phy col starts pos and goes neg
phy0 = -2*pi*tmp;
phy = repmat(phy0, [1 epipar.nkcols]);

phx = reshape1d(phx');
phy = reshape1d(phy');


return
