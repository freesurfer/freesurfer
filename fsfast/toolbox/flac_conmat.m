function flacnew = flac_conmat(flac,nthcon)
% flacnew = flac_conmat(flac,nthcon)
%
% Computes the contrast matrix for the nth contrast. 
% Adds to flac.con(nthcon).C
%
%


%
% flac_conmat.m
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

flacnew = [];

if(nargin ~= 2) 
  fprintf('flacnew = flac_conmat(flac,nthcon)\n');
  return;
end

if(nthcon > length(flac.con))
  fprintf('ERROR: flac_conmat: nthcon = %d\n',nthcon);
  return;
end

con = flac.con(nthcon);
nevs = length(flac.ev);

if(con.sumev == 1 & con.sumevreg == 1)
  % t-test, C has only one row
  C = [];
  for nthev = 1:nevs
    [evw evrw hit] = flac_evconw(flac,nthev,nthcon);    
    if(isempty(evw)) return; end
    C = [C evw*evrw];
  end
end

if(con.sumev == 1 & con.sumevreg == 0)
  % Sum the EVs but not the EV Regs
  % Determine J, the number of rows in C
  J = 0;
  for nthev = 1:nevs
    [evw evrw hit] = flac_evconw(flac,nthev,nthcon);    
    if(isempty(evw)) return; end    
    if(hit)
      if(J==0) J = flac.ev(nthev).nreg;
      else
	if(J ~= flac.ev(nthev).nreg)
	  fprintf('ERROR: flac_conmat: %s: (SumEv=1,SumEVReg=0)\n',con.name);
	  fprintf('Two or more EVs have different numbers of regressors\n');
	  return;
	end
      end
    end % if hit
  end % for nthev

  % Now determine C
  C = [];
  for nthev = 1:nevs
    [evw evrw hit] = flac_evconw(flac,nthev,nthcon);    
    if(hit) Cev = evw*diag(evrw);
    else    Cev = zeros(J,flac.ev(nthev).nreg);
    end
    C = [C Cev];
  end

end

if(con.sumev == 0 & con.sumevreg == 1)
  % Determine J, the number of rows in C
  J = 0;
  for nthev = 1:nevs
    [evw evrw hit] = flac_evconw(flac,nthev,nthcon);    
    if(isempty(evw)) return; end
    if(hit) J = J + 1; end
  end 
  % Now determine C
  C = [];
  nthhit = 0;
  for nthev = 1:nevs
    [evw evrw hit] = flac_evconw(flac,nthev,nthcon);    
    Cev = zeros(J,flac.ev(nthev).nreg);
    if(hit)
      nthhit = nthhit + 1;
      Cev(nthhit,:) = evw*evrw;
    end
    C = [C Cev];
  end
end

if(con.sumev == 0 & con.sumevreg == 0)
  % Just put everything on the diagonal and remove the zeros later
  c = [];
  for nthev = 1:nevs
    [evw evrw hit] = flac_evconw(flac,nthev,nthcon);    
    c = [c evw*evrw];
  end 
  C = diag(c);
end

% Now remove rows of C that are all 0
J = size(C,1);
keeprows = [];
for row = 1:J
  ind = find(C(row,:) ~= 0);
  if(~isempty(ind)) keeprows = [keeprows row]; end
end
C = C(keeprows,:);

if(isempty(C))
  fprintf('ERROR: flac_conmat: %s: all weights are 0\n',con.name);
  return;
end

flacnew = flac;
flacnew.con(nthcon).C = C;


return














