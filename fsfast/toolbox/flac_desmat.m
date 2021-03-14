function flacnew = flac_desmat(flac,IRFOnly)
% flacnew = flac_desmat(flac,<IRFOnly>)
%
% Builds design matrices for each EV and performs the horizontal
% concatenation. Requires that flac.ntp and flac.ev(n).st already be
% set. If a nonpar is used, the matrix must already be set. These are
% done by flac_customize but could be done in some other way (allows
% for optseq-type optimization). If IRFOnly=1 and the FIR matrix
% exists, then the FIR matrix is not rebuilt (for speed).
%


%
% flac_desmat.m
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

flacnew = [];

if(nargin < 1 || nargin > 2)
 fprintf('flacnew = flac_desmat(flac,IRFOnly)\n');
 return;
end

if(nargin == 1) IRFOnly = 0; end

flacnew = flac;
flacnew.X = [];

nev = length(flac.ev);
for nthev = 1:nev
  ev = flac.ev(nthev);
  
  if(ev.ishrf)  
    % HRF Regressors
    st = ev.st; % Delay has already been added
    if(~isfield(ev,'Xfir')) ev.Xfir = []; end
    if(~IRFOnly | isempty(ev.Xfir))
      TER = ev.psdwin(3);
      ssr = round(flac.TR/TER);
      ntpssr = flac.ntp * ssr;
      useweight = 1;
      flacnew.ev(nthev).Xfir = fast_st2fir(st,ntpssr,TER,ev.psdwin,useweight);
      if(isempty(flacnew.ev(nthev).Xfir)) 
	fprintf('ERROR: creating FIR design matrix for %s\n',...
		flacnew.ev(nthev).name);
	flacnew = [];
	return; 
      end
    end
    [Xirf tirf] = flac_ev2irf(ev,TER,flac.RefEventDur);
    if(isempty(Xirf)) 
      fprintf('ERROR: creating IRF design matrix for %s\n',...
	      flacnew.ev(nthev).name);
      flacnew = [];
      return; 
    end
    flacnew.ev(nthev).Xirf = Xirf;
    flacnew.ev(nthev).tirf = tirf;
    Xtmp = flacnew.ev(nthev).Xfir * flacnew.ev(nthev).Xirf;
    flacnew.ev(nthev).X = Xtmp(1:ssr:end,:);
    flacnew.ev(nthev).Xfir = []; % To save memory
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
      if(strcmp(flacnew.direction,'neg')) 
	X(:,[2:2:end]) = -X(:,[2:2:end]); % negate imaginary part
      end
      if(flac.IsRetinotopy & ~strcmp(flac.stimtype,ev.name)) 
	X = zeros(size(X));
      end
      flacnew.ev(nthev).X = X;
     case {'abret'}  
      % Impelemnts old ABBlocked and retinotopy
      % 12 Regressors: fundR fundI harmR harmI 
      % fundR- fundI- fundR+ fundI+ harmR- harmI- harmR+ harmI+
      period  = ev.params(1);
      fund = 1/period;
      [fftaxis delta] = fast_fftaxis(flac.ntp,flac.TR);
      freqs = [fund 2*fund fund-delta fund+delta 2*fund-delta 2*fund+delta];
      t = flac.TR*[0:flac.ntp-1]';
      ph = 2*pi*t*freqs;
      X = zeros(flac.ntp,12);
      X(:,1:2:end) = cos(ph); % real
      X(:,2:2:end) = sin(ph); % imag
      if(strcmp(flacnew.direction,'neg')) 
	X(:,[2:2:end]) = -X(:,[2:2:end]); % negate imaginary part
      end
      X = X - repmat(mean(X),[flac.ntp 1]); % remove mean
      if(flac.IsRetinotopy & ~strcmp(flac.stimtype,ev.name)) 
	X = zeros(size(X));
      end
      flacnew.ev(nthev).X = X;
     case {'hpf'}  
      CutoffHz  = ev.params(1);
      fftaxis = fast_fftaxis(flac.ntp,flac.TR);
      ind = find(fftaxis > 0 & fftaxis <= .75*CutoffHz);
      t = flac.TR*[0:flac.ntp-1]';
      ph = 2*pi*t*fftaxis(ind);
      %X = [cos(ph) sin(ph)]; % DFT
      X = [cos(ph)]; % DCT
      X = X - repmat(mean(X),[flac.ntp 1]);
      [u s] = fast_svd(X);
      pvs = 100*(diag(s))/sum(diag(s)); % percent var explained
      indkeep = find(pvs > 1); % keep components over 1%
      X = u(:,indkeep);
      flacnew.ev(nthev).X = X;
      flacnew.ev(nthev).nreg = length(indkeep); % Needed!
     case {'hpf+poly'}  
      CutoffHz  = ev.params(1);
      polyorder = ev.params(2);
      fftaxis = fast_fftaxis(flac.ntp,flac.TR);
      ind = find(fftaxis > 0 & fftaxis <= .75*CutoffHz);
      t = flac.TR*[0:flac.ntp-1]';
      ph = 2*pi*t*fftaxis(ind);
      %X = [cos(ph) sin(ph)]; % DFT
      X = [cos(ph)]; %DCT
      X = X - repmat(mean(X),[flac.ntp 1]);
      Xp = fast_polytrendmtx(1,flac.ntp,1,polyorder);
      X = [X Xp(:,2:end)];
      [u s] = fast_svd(X);
      pvs = 100*(diag(s))/sum(diag(s));
      %cpvs = 100*cumsum(diag(s))/sum(diag(s));
      indkeep = find(pvs > 1);
      X = u(:,indkeep);
      flacnew.ev(nthev).X = X;
      flacnew.ev(nthev).nreg = length(indkeep); % Needed!
     case {'nyquist'}  
      X = 2*(rem([1:flac.ntp]',2)-.5);
      flacnew.ev(nthev).X = X;
     case {'asl'}  
      X = zeros(flac.ntp,1);
      X(1+ev.params(1):2:end) = 1;
      flacnew.ev(nthev).X = X;
     case {'aareg'}  
      nkeep = ev.params(1);
      fmax = ev.params(2);
      if(fmax < 0) fmax = 60/90; end % 90 beats per min
      fdelta = ev.params(3);
      if(fdelta < 0) fdelta = []; end % Use fftdelta
      X = fast_aareg(flac.ntp,flac.TR,fmax,fdelta);
      flacnew.ev(nthev).X = X(:,1:nkeep);
     case {'nonpar','selfregseg'}
      % Nonpar X must be already loaded with flac_customize.
      if(isempty(flacnew.ev(nthev).X))
	fprintf('ERROR: empty %s matrix for %s\n',ev.model,ev.name);
	flacnew = [];
	return;
      end
     case {'texclude'}  
      % texclude X must be already loaded with flac_customize,
      % OR it can also be empty.
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














