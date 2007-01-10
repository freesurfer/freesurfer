function [evw, evrw, hit] = flac_evconw(flac,nthev,nthcon)
% [evw, evrw, hit] = flac_evconw(flac,nthev,nthcon)
%
% Retuns the weights of the nth EV as found in the nth contrast. If
% the EV is not specified in the contrast, hit = 0;
%
%


%
% flac_evconw.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.3 $
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

evw = [];
evrw = [];
hit = 0;

if(nargin ~= 3)
  fprintf('[evw evrw] = flac_evconw(flac,nthev,nthcon)\n');
  return;
end

% Nth contrast
con = flac.con(nthcon);

% Number of evs as found in the nth contrast.
nevscon = length(con.ev);

% Find the index of the EV in the contrast corresponding to the nth EV
for nthconev = 1:nevscon

  if(strcmp(con.ev(nthconev).name,flac.ev(nthev).name))
    evw = con.ev(nthconev).evw; % Weight for this EV
    if(~isempty(con.ev(nthconev).evrw))
      % Reg Weights specified for this EV
      evrw = con.ev(nthconev).evrw;
    elseif(~isempty(con.evrw))
      % Default Reg Weights for contrast
      evrw = con.evrw;
    else
      % Final default all ones
      evrw = ones(1,flac.ev(nthev).nreg);
    end
    if(length(evrw) ~= flac.ev(nthev).nreg)
      fprintf('ERROR: flac_evconw(): %s, %s evrw dim mismatch\n',...
	      con.name,flac.ev(nthev).name);
      fprintf('This condition has %d regressors, but evrm has %d\n',...
	      flac.ev(nthev).nreg,length(evrw));
      evw = []; evrw = [];
      return;
    end
    hit = 1;
    return;
  end

end

% The nth EV was not specified in the contrast, so return zeros
hit = 0;
evw = 0;
evrw = zeros(1,flac.ev(nthev).nreg);

return;


















