function flacnew = flac_desmat(flac)
% flacnew = flac_desmat(flac)
%
% Builds design matrices for each EV and performs the horizontal
% concatenation. Requires that flac.ntp and flac.ev(n).st already be
% set. If a nonpar is used, the matrix must already be set. These are
% done by flac_customize but could be done in some other way (allows
% for optseq-type optimization).
%
%


%
% flac_desmat.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/12/04 18:02:58 $
%    $Revision: 1.12 $
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

if(nargin ~= 1)
 fprintf('flacnew = flac_desmat(flac)\n');
 return;
end

flacnew = flac;
flacnew.X = [];

nev = length(flac.ev);
for nthev = 1:nev
  ev = flac.ev(nthev);
  
  if(ev.ishrf)  
    % HRF Regressors
    st = ev.st; % Delay has already been added
    flacnew.ev(nthev).Xfir = fast_st2fir(st,flac.ntp,flac.TR,ev.psdwin,1);
    if(isempty(flacnew.ev(nthev).Xfir)) 
      fprintf('ERROR: creating FIR design matrix for %s\n',...
	      flacnew.ev(nthev).name);
      flacnew = [];
      return; 
    end
    flacnew.ev(nthev).Xirf = flac_ev2irf(flac,nthev);
    if(isempty(flacnew.ev(nthev).Xirf)) 
      fprintf('ERROR: creating IRF design matrix for %s\n',...
	      flacnew.ev(nthev).name);
      flacnew = [];
      return; 
    end
    flacnew.ev(nthev).X = flacnew.ev(nthev).Xfir * flacnew.ev(nthev).Xirf;
  else
    switch(ev.model)
     case {'baseline'}
      flacnew.ev(nthev).X = ones(flac.ntp,1);
     case {'polynomial'}
      polyorder = ev.params(1);
      X = fast_polytrendmtx(1,flac.ntp,1,polyorder);
      flacnew.ev(nthev).X = X(:,2:end);
     case {'fourier'}  
      period     = ev.params(1);
      nharmonics = ev.params(2);
      tdelay     = ev.params(3); % Need to add
      X = fast_fourier_reg(period,flac.ntp,flac.TR,nharmonics);
      flacnew.ev(nthev).X = X;
     case {'nyquist'}  
      X = 2*(rem([1:flac.ntp]',2)-.5);
      flacnew.ev(nthev).X = X;
     case {'asl'}  
      X = zeros(flac.ntp,1);
      X(1+ev.params(1):2:end) = 1;
      flacnew.ev(nthev).X = X;
     case {'nonpar','selfregseg'}  
      % Nonpar X must be already loaded with flac_customize.
      if(isempty(flacnew.ev(nthev).X))
	fprintf('ERROR: empty %s matrix for %s\n',ev.model,ev.name);
	flacnew = [];
	return;
      end
     otherwise
      fprintf('ERROR: model %s unrecognized\n');
      flacnew = [];
      return;
    end
  end
%keyboard
  
  flacnew.X = [flacnew.X flacnew.ev(nthev).X];
  
end % for ev


return;














