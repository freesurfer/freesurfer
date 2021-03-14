function h = yvplt_acffit(varargin)
% h = yvplt_acffit(t,y)
%
% Computes, fits, and plots autocorrelation function. This
% can be used with yakview by adding -rawfunc yvplt_acffit
%


%
% yvplt_acffit.m
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

DTNONE = 0;
DTPF   = 1;
DTHPF  = 2;
DTEXT  = 3;

h = [];

if(nargin == 0)
  fprintf('USAGE: h = yvplt_acffit(string,options)\n');
  return;
end

cbflag = varargin{1};

if(~isstr(cbflag))
  fprintf('First argument must be a string');
  return;
end

switch(cbflag)
  case 'init',
    if(nargin ~= 3)
      msg = 'USAGE: hdrview(init,t,y)';
      qoe(msg);error(msg);
    end
    ud.yvplt_acffit = 1;
    ud.arorder = 3;
    ud.pforder = 2;
    ud.nskip = 0;
    ud.hpfcutoff = 2;
    ud.whiten = 0;
    ud.showres = 0; % show residuals (only with external)
    ud.showbeta = 0; % show beta (or fft)
    ud.showhist = 0; % show histogram
    ud.pct = 0; % show percent signal change
    ud.acfshow1 = 0; % show acf from first lag
    ud.XX = [];
    ud.XFile = [];
    ud.XDir = [];
    ud.ParDir = [];
    ud.ParFile = [];
    ud.par = [];
    ud.hX = [];
    ud.dtmethod = DTPF;
    ud.ThisHDRView = gcf;
    ud.t = varargin{2};
    ud.y = varargin{3};
    ud.hCtl = ''; % Currently not used
    set(gcf,'KeyPressFcn','yvplt_acffit(''kbd'');');

  case 'plot',
    if(nargin ~= 3)
      msg = 'USAGE: yvplt_acffit(plot,t,y);';
      qoe(msg);error(msg);
    end
    ud = get(gcf,'UserData'); 
    ud.t = varargin{2};
    ud.y = varargin{3};

  case 'kbd',
    c = get(gcf,'CurrentCharacter'); 
    %fprintf(1,'Key %s\n',c);
    ud = get(gcf,'UserData'); 

    switch(c)

      case {'a'},
        ttl  = 'Enter AR Order';
        prompt = {'Order:'};
        lines  = 1;
        def = {sprintf('%7.4f',ud.arorder)};
        answer   = inputdlg(prompt,ttl,lines,def);
        arorder = sscanf(answer{1},'%f');
        if(arorder <= 0)
          msg = sprintf('AROrder = %d <= 0',arorder);
          errordlg(msg);
          return;
        end
        ud.arorder = arorder;

      case {'d'},
        ttl  = 'Detrending Method';
        liststr{1} = 'None';
        liststr{2} = 'Poly';
        liststr{3} = 'HPF';
        liststr{4} = 'External';
        [s ok] = listdlg('Name',ttl,'ListString',liststr,...
                 'SelectionMode','single','InitialValue',ud.dtmethod+1);
        %fprintf('s=%d, ok = %d\n',s,ok);
        if(~ok) return; end
        if(s-1 == ud.dtmethod) return; end
        if(s-1 == DTEXT & isempty(ud.XX)) 
          msg = 'You must load an external matrix first. Press x.';
          errordlg(msg);
          return; 
        end
        ud.dtmethod = s-1;

      case {'p'},
        ttl  = 'Enter PolyFit Order';
        prompt = {'Order:'};
        lines  = 1;
        def = {sprintf('%7.4f',ud.pforder)};
        answer   = inputdlg(prompt,ttl,lines,def);
        ud.pforder = sscanf(answer{1},'%f');

      case {'n'},
        ttl  = 'Enter NSkip';
        prompt = {'NSkip:'};
        lines  = 1;
        def = {sprintf('%7.4f',ud.nskip)};
        answer   = inputdlg(prompt,ttl,lines,def);
        ud.nskip = sscanf(answer{1},'%f');

      case {'h'},
        ttl  = 'Enter HPF CutOff (TRs)';
        prompt = {'CutOff:'};
        lines  = 1;
        def = {sprintf('%7.4f',ud.hpfcutoff)};
        answer   = inputdlg(prompt,ttl,lines,def);
        hpfcutoff = sscanf(answer{1},'%f');
        if(hpfcutoff <= 0)
          msg = sprintf('HPF CutOff = %g <= 0',hpfcutoff);
          errordlg(msg);
          return;
        end
        ud.hpfcutoff = hpfcutoff;

      case {'i'},
        if(isempty(ud.XX)) 
          msg = 'You must load an external matrix first. Press x.';
          errordlg(msg);
          return; 
        end
        hCur = gcf;
        if(isempty(ud.hX)) ud.hX = figure; end
        figure(ud.hX);
        imagesc(ud.XX.Xfinal); colorbar;
        figure(hCur);

      case {'x'},
        curdir = pwd;
        if(~isempty(ud.XDir)) cd(ud.XDir); end
        [fname pname] = uigetfile('*.mat','Load External Matrix');
        cd(curdir);
        XFile = sprintf('%s/%s',pname,fname);
        fprintf('XFile: %s\n',XFile);
        XX = load(XFile);
        if(isempty(XX))
          msg = sprintf('Error loading %s',XFile);
          errordlg(msg);
          return;
        end
        ud.dtmethod = DTEXT;

        if(~isfield(XX,'Xfinal'))
          msg = sprintf('ERROR: matfile does not have Xfinal variable');
          errordlg(msg);
          return;
        end
        if(~isfield(XX,'Nnnc'))
          msg = sprintf('ERROR: matfile does not have Nnnc variable');
          errordlg(msg);
          return;
        end
        if(~isfield(XX,'Navgs_per_cond'))
          msg = sprintf('ERROR: matfile does not have Navgs variable');
          errordlg(msg);
          return;
        end
        if(size(XX.Xfinal,1) ~= length(ud.t))
          msg = sprintf('ERROR: X is wrong size (%d,%d)',...
             size(XX.Xfinal,1),length(ud.t))
          errordlg(msg);
          return;
        end
        if(isempty(pname)) ud.XDir  = pwd; 
        else               ud.XDir  = pname;
        end
        ud.XFile = XFile;
        ud.XX = XX;

      case {'m'},
        curdir = pwd;
        if(~isempty(ud.ParDir)) cd(ud.ParDir); end
        [fname pname] = uigetfile('*.pdg','Load Paradigm');
        cd(curdir);
        ParFile = sprintf('%s/%s',pname,fname);
        fprintf('ParFile: %s\n',ParFile);
        par = fmri_ldpar(ParFile);
        if(isempty(par))
          msg = sprintf('Error loading %s',ParFile);
          errordlg(msg);
          return;
        end
        ud.ParFile = ParFile ;
        ud.par = par;
        if(isempty(pname)) ud.ParDir  = pwd; 
        else               ud.ParDir  = pname;
        end

      case {'c'},
        ud.pct = ~ud.pct;

      case {'f'},
        ud.showhist = ~ud.showhist;

      case {'r'},
        ud.showres = ~ud.showres;

      case {'1'},
        ud.acfshow1 = ~ud.acfshow1;

      case {'b'},
        if(isempty(ud.XX)) return; end
        ud.showbeta = ~ud.showbeta;

      case {'w'},
        ud.whiten = ~ud.whiten;
        if(ud.whiten) fprintf('Whitening On\n');
        else fprintf('Whitening Off\n');
        end

    end % switch(c) %
end

%-- Update the user data so its there for the next time ----%
set(gcf,'UserData',ud); 
%--------------------------------------------------------%

arorder   = ud.arorder;
pforder   = ud.pforder;
hpfcutoff =  ud.hpfcutoff;
dtmethod  =  ud.dtmethod;
t = ud.t;
if(~isempty(ud.XX)) X = ud.XX.Xfinal;
else X = [];
end
y = ud.y;
par = ud.par;

if(ud.nskip > 0)
  mtmp = ud.nskip+1:length(y);
  y = y(mtmp);
  if(~isempty(X)) X = X(mtmp,:); end
  t = t(mtmp);
  if(~isempty(par))
     i = find(par(:,1) > ud.XX.TR*ud.nskip);
     par = ud.par(i,:);
  end
end


if(ud.pct) y = 100*y/mean(y); end

y0 = y;

TR = mean(diff(t));

ntrs = length(y);
nmax = round(ntrs/2);

freqmax = (1/TR)/2;  % Nyquist
deltafreq = freqmax/(ntrs/2);  % Measured from 0 to Nyquist
freqaxis = deltafreq * [0:nmax-1]'; %'

% Detrend data %
switch(dtmethod)
  case DTNONE,
    r = y;
  case DTEXT,
    T = X*inv(X'*X)*X';
    R = eye(ntrs) - T;
    s = T*y;
    r = R*y;
    beta = (inv(X'*X)*X')*y;
    if(ud.whiten)
      arparams = arburg(r,arorder);
      acf = fast_ar2acf(arparams,ntrs);
      L = fmri_acorr2covmtx(acf,ntrs,1);
      W = chol(inv(L));
      Z = W*X;
      T = X*inv(Z'*Z)*Z'*W;
      R = W*(eye(ntrs) - T);
      s = T*y;
      r = R*y;
      beta = (inv(Z'*Z)*Z'*W)*y;
    end
    NavgsTot = ud.XX.Nnnc * ud.XX.Navgs_per_cond;
    beta_task = beta(1:NavgsTot);
    fprintf('beta: ');
    fprintf('%g ',beta_task);
    fprintf('\n');
    ttle = sprintf('Raw Time Course and External Estimate');
  case DTPF,
    if(pforder >= 0)
      D = fast_polytrendmtx(1,ntrs,1,pforder);
      E = eye(ntrs) - D*inv(D'*D)*D';
      r = E*y;
    end
    ttle = sprintf('Raw Time Course after polynomial(%d) detrending',pforder);
  case DTHPF,
    G = fast_mkgausmtx(hpfcutoff,ntrs);
    H = eye(ntrs)-G;
    r = H*y;
    ttle = sprintf('Raw Time Course after HPF(%g) detrending',hpfcutoff);
end

% AR estimation with Burg Method
%[ayw yvaryw]         = aryule(y,arorder);
[ab  yvarb ]          = arburg(r,arorder);
%[a_cov  yvar_cov ]   = arcov(y,arorder);
%[a_mcov  yvar_mcov ] = armcov(y,arorder);

acf0 = xcorr(r,nmax,'unbiased');
acf1 = acf0/max(acf0);
acf = acf1(nmax+1:length(acf0));
acf = acf';%'

% Compute ACF based on AR coefficiencts
%acfyw    = fast_ar2acf(ayw,nmax+1);
acfb     = fast_ar2acf(ab,nmax+1);
%acfcov   = fast_ar2acf(a_cov,nmax+1);
%acfmcov  = fast_ar2acf(a_mcov,nmax+1);

h = subplot(3,1,1);
if(dtmethod ~= DTEXT | ud.showres) 
  if(~ud.showhist) 
    plot(t,r,t,zeros(size(t)),'-.');
  else
    [rh xrh] = hist(r,round(ntrs/10));
  end
else 
  if(isempty(ud.par)) plot(t,y0,t,s,'r-');
  else
    rng = max(s)-min(s);
    plot(t,y0,t,s,'r-',par(:,1),rng*par(:,2)+min(s));
  end
  legend('Measured','Estimate');
end
title(ttle);

subplot(3,1,2);
ff = abs(fft(r));
ff = ff(1:round(ntrs/2));
if(ud.showbeta == 0)
  plot(freqaxis,ff);
  ttle = sprintf('FFT of Raw Time Course after detrending',pforder);
  title(ttle);
else
  mm = 1:length(beta_task);
  plot(mm,beta_task,mm,zeros(size(mm)),'k-.');
  title('Beta');
end

subplot(3,1,3);
if(ud.acfshow1) nn = 0:length(acf)-1;
else            nn = 1:length(acf)-1;
end
plot(nn,acf(nn+1),nn,acfb(nn+1));
legend('Actual','Burg');
%plot(nn,acf,nn,acfyw,nn,acfb,nn,acfcov,nn,acfmcov);
%legend('Actual','YW','Burg','Cov','MCov');
ttle = sprintf('ACF with AR%d fit',arorder);
title(ttle);

subplot(3,1,1);

%fprintf('AROrder = %d, PFOrder = %d\n',arorder,pforder);
%fprintf('YW     %6.4f ',yvaryw);
%fprintf('%6.4f ',ayw);
%fprintf('\n');

fprintf('Burg   %6.4f ',yvarb);
fprintf('%6.4f ',ab);
fprintf('\n');

%fprintf('Cov    %6.4f ',yvar_cov);
%fprintf('%6.4f ',a_cov);
%fprintf('\n');

%fprintf('MCov   %6.4f ',yvar_mcov);
%fprintf('%6.4f ',a_mcov);
%fprintf('\n');

return;
