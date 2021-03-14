function h = yvpltraw(varargin)
% h = yvpltraw(varargin)
% 
%  ('-init',hdatfile)
%  ('-base',base)
%  ('-plot',c,r,s)
%
% Plots raw time courses from an analysis.
%

%
% yvpltraw.m
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


% Things to do:
% 1. allow user to spec analysis
% 2. per-run or full session


if(nargin == 0)
  fprintf('USAGE: h = yvpltraw(string,options)\n');
  return;
end

cbflag = varargin{1};

if(~isstr(cbflag))
  fprintf('First argument must be a string');
  return;
end

switch(cbflag)
  case 'init',
    if(nargin ~= 2 & nargin ~= 3)
      msg = 'USAGE: yvpltraw(init,hdatfile,<hparent>)';
      qoe(msg);error(msg);
    end
    ud.yvpltraw = 1;
    ud.percent  = 0;
    ud.hdatfile = varargin{2};
    if(nargin == 3) ud.hparent = varargin{3};
    else ud.hparent = [];
    end

    ud = load_datfile(ud);
    if(isempty(ud)) return; end
    ud.nthrun = 1;
    ud.showpmf = 0; % show partial model fit
    ud.perrun = 0;  % analyze on a per-run basis
    ud.anadir = fast_dirname(ud.hdatfile);
    xfile = sprintf('%s/X.mat',ud.anadir);
    ud.XX = load(xfile);
    if(isempty(ud.XX))
      fprintf('ERROR: could not load %s\n',xmatfile);
      return;
    end
    ud.fsd = fast_dirname(ud.anadir);
    ud = set_stem(ud);
    [ud err] = set_matrices(ud);
    ud.base = 0;
    ud.c = floor(ud.ncols/2) + ud.base;
    ud.r = floor(ud.nrows/2) + ud.base;
    ud.s = floor(ud.nslices/2) + ud.base;
    ud.showpar = 0; % show partial model fit
    ud.showyhat = 1;
    ud.showraw = 1;
    ud.showlegend = 1;
    ud.fixaxes = 0;
    ud.axis = [];
    ud.tpx = [];
    ud.rescale = []; % this should not be necessary
    [ud err] = load_parfile(ud);
    [ud err] = load_tpx(ud);
    ud = load_rescale(ud);
    ud = load_voxel(ud);
    h = figure;
    set(gcf,'KeyPressFcn',         'yvpltraw(''kbd'');');
    set(gcf,'WindowButtonDownFcn', 'yvpltraw(''wbd'')');
    ud.curpostxt = uicontrol('Style', 'text','Position',  [1 1 250 20]);

  case 'plot',
    if(nargin ~= 4)
      msg = 'USAGE: yvplt_acffit(plot,c,r,s);';
      qoe(msg);error(msg);
    end
    ud = get(gcf,'UserData'); 
    ud.c = varargin{2};
    ud.r = varargin{3};
    ud.s = varargin{4};
    % check crs %
    ud = load_voxel(ud);

  case 'wbd',
    figure(gcbo);
    ud = get(gcf,'UserData'); 
    tvz = get(gca,'CurrentPoint');
    t = tvz(1,1);
    v = tvz(1,2);
    [m tind] = min(abs(ud.t-t));
    [m tparind] = min(abs(ud.par(:,1)-t));
    tpar = ud.par(tparind,1);
    cond = ud.par(tparind,2);
    fprintf('t=%g, i=%d, v=%g, raw=%g, yhat=%g, cond=%d (%g)\n',...
	    t,tind,v,ud.raw(tind),ud.yhat(tind),cond,tpar);
    ud.curpos = [t v];
    set(gcf,'UserData',ud); 
    setcurpostxt(ud);
    if(~isempty(ud.hparent)) figure(ud.hparent); end
    return;

  case 'kbd',
    figure(gcbo);
    c = get(gcf,'CurrentCharacter'); 
    %fprintf(1,'Key %s\n',c);
    ud = get(gcf,'UserData'); 

    switch(c)
    case 'h' % Help
      s = sprintf('Keypress Commands: \n'); 
      s = sprintf(' %sa - toggle between per-run and all runs\n',s);
      s = sprintf(' %sf - toggle fixing of axis range\n',s);
      s = sprintf(' %sm - toggle partial-model fit\n',s);
      s = sprintf(' %sp - toggle display of paradigm\n',s);
      s = sprintf(' %sr - select new run to display \n',s);
      s = sprintf(' %ss - select new functional stem\n',s);
      s = sprintf(' %sv - save raw timecourse to ascii file\n',s);
      msgbox(s,'Help','help');
      if(~isempty(ud.hparent)) figure(ud.hparent); end
      return;

      case {'a'},
        ud.perrun = ~ud.perrun;
        [ud err] = load_voxel(ud);
        [ud err] = set_matrices(ud);
      case {'c'}, % percent signal change
        ud.percent = ~ud.percent;
	[ud err] = load_voxel(ud);
      case {'f'},
        ud.fixaxes = ~ud.fixaxes;
        ud.axis = axis;
      case {'l'},
        ud.showlegend = ~ud.showlegend;
      case {'m'},
        ud.showpmf = ~ud.showpmf;
        [ud err] = set_matrices(ud);
        if(err) return; end
      case {'='}, % advance to the next run
        ud.nthrun = ud.nthrun + 1;
	if(ud.nthrun > length(ud.ad.runlist))
	  ud.nthrun = 1;
	end
        [ud err] = set_stem(ud);
        [ud err] = set_matrices(ud);
        [ud err] = load_voxel(ud);
        [ud err] = load_parfile(ud);
      case {'-'}, % goto prev run
        ud.nthrun = ud.nthrun - 1;
	if(ud.nthrun < 1)
	  ud.nthrun = length(ud.ad.runlist);
	end
        [ud err] = set_stem(ud);
        [ud err] = set_matrices(ud);
        [ud err] = load_voxel(ud);
        [ud err] = load_parfile(ud);
      case {'p'},
        ud.showpar = ~ud.showpar;
      case {'r'},
        ttl  = 'Show Run';
        for n = 1:ud.ad.Nruns
  	  liststr{n} = sprintf('%03d',ud.ad.runlist(n));
        end
        [s ok] = listdlg('Name',ttl,'ListString',liststr,...
                 'SelectionMode','single','InitialValue',ud.nthrun);
        %fprintf('s=%d, ok = %d\n',s,ok);
        if(~ok) return; end
        if(s == ud.nthrun) return; end
        ud.nthrun = s;
        [ud err] = set_stem(ud);
        [ud err] = set_matrices(ud);
        [ud err] = load_voxel(ud);
        [ud err] = load_parfile(ud);

      case {'s'},
        ttl  = 'Enter Functional Stem';
        prompt = {'Stem:'};
        lines  = 1;
        def = {basename(ud.stem)};
        answer   = inputdlg(prompt,ttl,lines,def);
        oldstem = ud.ad.funcstem;
        ud.ad.funcstem = answer{1};
        [ud err] = set_stem(ud);
        if(err) 
           ud.ad.funcstem = oldstem;
           return; 
        end
        [ud err] = set_matrices(ud);
        [ud err] = load_voxel(ud);

      case {'v'},
        svfile = sprintf('%s-%03d-c%02d-r%02d-s%02d.dat',...
			 ud.ad.funcstem,ud.ad.runlist(ud.nthrun),...
			 ud.c, ud.r, ud.s);
        [fname pname] = uiputfile(svfile,'Save Raw Time Course');
        if(fname == 0) return; end
        svfile = sprintf('%s/%s',pname,fname);
        fid = fopen(svfile,'w');
        if(fid == -1)
          fprintf('ERROR: could not open %s for writing\n',svfile);
          return;
        end
        tmp = [ud.t ud.raw ud.yhat];
        fprintf(fid,'%g %g %g\n',tmp'); %'
	fclose(fid);

    end % switch(c) %

end
%+++++++++++++++++++++++++++++++++++%

nf = ud.nframes(ud.nthrun);
ud.t = ud.ad.TR*[0:nf-1]'; %'

if(ud.perrun)
  ud.yhat = ud.Trun*ud.v ;
  if(ud.showpmf) ud.raw = ud.Rrun*ud.v;
  else           ud.raw = ud.v;
  end
  ud.resvar = sum((ud.v - ud.yhat).^2)/trace(ud.Rrun);
else
  yhatall = ud.T*ud.vall;
  if(ud.showpmf) rawall = ud.R*ud.vall;
  else           rawall = ud.vall;
  end
  if(ud.nthrun == 1)   n1 = 1;
  else  n1 = sum(ud.nframes(1:ud.nthrun-1))+1;
  end
  n2 = n1 + ud.nframes(ud.nthrun)-1;
  ud.yhat = yhatall(n1:n2);
  ud.raw  = rawall(n1:n2);
  ud.resvar = sum((ud.vall - yhatall).^2)/trace(ud.R);
end

if(~isempty(ud.tpx) & ~isempty(ud.tpx(ud.nthrun).incl))
  ud.t = ud.t(ud.tpx(ud.nthrun).incl);
  ud.raw = ud.raw(ud.tpx(ud.nthrun).incl);
  ud.yhat = ud.yhat(ud.tpx(ud.nthrun).incl);
end


if(~ud.showpar)
  plot(ud.t,ud.raw,'r-',ud.t,ud.yhat,'b-');
  if(ud.showlegend) 
    if(ud.showpmf) legend('Residual','Task-Related'); 
    else           legend('Raw','Best-Fit'); 
    end
  end
else
  pltmin = min([ud.raw;ud.yhat]);
  pltmax = max([ud.raw;ud.yhat]);
  parmax = max(ud.par(:,2));
  partrace = .9*(pltmax-pltmin)*ud.par(:,2)/parmax + pltmin;
  plot(ud.t,ud.raw,'r-',ud.t,ud.yhat,'b-',...
       ud.par(:,1),partrace,'kx');
  %if(ud.ad.runlist(ud.nthrun)== 9) keyboard; end
  if(ud.showlegend) 
    if(ud.showpmf) legend('Residual','Task-Related','Paradigm'); 
    else           legend('Raw','Best-Fit','Pardigm');  
    end
  end
end

if(ud.showpmf)
  hold on;
  plot(ud.t,zeros(size(ud.t)),'k-.');
  hold off;
end

if(ud.fixaxes) axis(ud.axis); end

if(ud.perrun) tmpstr = '(per-run)';
else tmpstr = '';
end

tit = sprintf('%s %s %03d/%s c=%d, r=%d, s=%d (resvar=%g)\n',...
	      basename(ud.anadir),tmpstr,...
	      ud.ad.runlist(ud.nthrun),ud.ad.funcstem,...
	      ud.c,ud.r,ud.s,ud.resvar);
title(tit);
xlabel('time (sec)');
set(gcf,'Name','FS-FAST Raw Time-Course Viewer (RTCV)');

set(gcf,'UserData',ud); 
if(~isempty(ud.hparent)) figure(ud.hparent); end
return;
%---------------------- End main --------------------------%
%------------------------------------------------------%
%--------------------------------------------------%

%-------------------------------------------%
function  ud = load_datfile(ud);
ad = fmri_lddat3(ud.hdatfile);
if(isempty(ad))
  msg = sprintf('Error loading %s',ud.hdatfile);
  errordlg(msg);
  ud = [];
  return;
end
if(ad.Version ~= 3)
  msg = 'ERROR: analysis header file is out-of-date. This is not a problem. You just need to re-run selxavg-sess to create a new header. Your results will not change.';
  errordlg(msg);
  ud = [];
  return;
end
ud.ad = ad;
return;

%-------------------------------------------%
function [ud, err] = load_parfile(ud)
  err = 0;
  run = ud.ad.runlist(ud.nthrun);
  parfile = sprintf('%s/%03d/%s',ud.fsd,run,ud.ad.parname);
  par = fmri_ldpar(parfile);
  if(isempty(par))
    fprintf('ERROR: loading info from %s\n',parfile);
    err = 1;
    return;
  end
  % remove zeros from par%
  indnz = find(par(:,2) ~= 0);
  ud.par = par(indnz,:);
  fprintf('INFO: parfile is %s\n',parfile);
return;


%-------------------------------------------%
function [ud, err] = set_stem(ud)
  err = 0;
  for nthrun = 1:ud.ad.Nruns,
    run = ud.ad.runlist(nthrun);
    stem = sprintf('%s/%03d/%s',ud.fsd,run,ud.ad.funcstem);
    [nrows ncols nf fs nslices endian bext] = fmri_bfiledim(stem);
    if(isempty(nrows))
      fprintf('ERROR: loading info from %s\n',stem);
      err = 1;
      return;
    end
    nframes(nthrun) = nf;
  end
  run = ud.ad.runlist(ud.nthrun);
  ud.stem = sprintf('%s/%03d/%s',ud.fsd,run,ud.ad.funcstem);
  ud.ncols = ncols;
  ud.nrows = nrows;
  ud.nslices = nslices;
  ud.nframes = nframes;

  fprintf('INFO: yvpltraw: stem = %s\n',ud.stem);
return;


%-------------------------------------------%
function [ud, err] = load_voxel(ud)
  %tic;
  err = 0;
  v = fast_ldbvoxel(ud.stem,ud.c,ud.r,ud.s,ud.base);
  if(isempty(v))
    fprintf('ERROR: yvpltraw: loading voxel %d %d %d (base=%d)\n',...
	   ud.c,ud.r,ud.s,ud.base);
    err = 1;
    return;
  end
  if(~isempty(ud.rescale))
    v = v*ud.rescale(ud.nthrun);
  end
  if(ud.percent)
    mnstem = sprintf('%s/h-offset',ud.anadir);
    vmn = fast_ldbvoxel(mnstem,ud.c,ud.r,ud.s,ud.base);
    v = 100*(v-vmn)/vmn;
  end
  ud.v = v;
  if(~ud.perrun) [ud err] = load_voxel_allruns(ud); end

  %fprintf('Time to load: %g\n',toc);
return;

%-------------------------------------------%
function [ud, err] = load_voxel_allruns(ud)
  err = 0;
  vall = [];
  for nthrun = 1:ud.ad.Nruns,
    run = ud.ad.runlist(nthrun);
    stem = sprintf('%s/%03d/%s',ud.fsd,run,ud.ad.funcstem);
    v = fast_ldbvoxel(stem,ud.c,ud.r,ud.s,ud.base);
    if(isempty(v))
      fprintf('ERROR: yvpltraw: loading voxel %d %d %d, run %d\n',...
	      ud.c,ud.r,ud.s,run);
      err = 1;
      return;
    end
    if(~isempty(ud.tpx) & ~isempty(ud.tpx(nthrun).excl))
      v(ud.tpx(nthrun).excl) = 0;
    end
    if(~isempty(ud.rescale))
      v = v*ud.rescale(nthrun);
    end
    vall = [vall; v];
  end
  if(ud.percent)
    mnstem = sprintf('%s/h-offset',ud.anadir);
    vmn = fast_ldbvoxel(mnstem,ud.c,ud.r,ud.s,ud.base);
    vall = 100*(vall-vmn)/vmn;
  end
    
  ud.vall = vall;
return;

%-------------------------------------------%
function [ud, err] = set_matrices(ud)
  err = 0;
  if(ud.nthrun == 1)   n1 = 1;
  else  n1 = sum(ud.nframes(1:ud.nthrun-1))+1;
  end
  n2 = n1 + ud.nframes(ud.nthrun)-1;
  c2 = ud.XX.Nnnc*ud.ad.Nh;

  X = ud.XX.Xfinal;
  Xruntmp = X(n1:n2,:);
  Xrun = [];
  for c = 1:size(Xruntmp,2)
    n = length(find(Xruntmp(:,c)~=0));
    if(n > 0) Xrun = [Xrun Xruntmp(:,c)]; end
  end

  Xpmf = X;
  Xpmfrun = Xrun;
  if(ud.showpmf)  
    Xpmfrun(:,c2+1:end) = 0;  
    Xpmf(:,c2+1:end) = 0;  
  end

  if(ud.perrun)
    % Need to put trap here for when X is ill-conditioned
    % because not enough conditions. Or is this fixed by
    % removing all cols that are zero?
    ud.Trun = Xpmfrun * inv(Xrun'*Xrun)* Xrun'; 
    ud.Rrun = eye(size(Xrun,1)) - Xrun * inv(Xrun'*Xrun)* Xrun'; 
  end

  ud.T = Xpmf * inv(X'*X)* X'; 
  ud.R = eye(size(X,1)) - X*inv(X'*X)* X'; 

  %if(ud.nthrun ~= 1) keyboard; end

return;
%----------------------------------------------------------%
function [ud, err] = load_tpx(ud)

err = 0;
if(isempty(ud.XX.tpxlist) & ud.ad.Nskip == 0) return; end

for nthrun = 1:ud.ad.Nruns,
  run = ud.ad.runlist(nthrun);
  stem = sprintf('%s/%03d/%s',ud.fsd,run,ud.ad.funcstem);
  if(~isempty(ud.XX.tpxlist))
    tpxfile = deblank(ud.XX.tpxlist(nthrun,:));
    if(~strcmp(tpxfile,'noexcl'))
      [nslices nrows ncols ntrs] = fmri_bvoldim(stem);
      tpxfile = sprintf('%s/%s',ud.fsd,tpxfile);
      [indTPExcl indTPIncl] = fast_ldtpexcl(tpxfile,ud.ad.TR,...
					    ntrs,ud.ad.Nskip);
    end
  else
    [nslices nrows ncols ntrs] = fmri_bvoldim(stem);
    indTPExcl = 1:ud.ad.Nskip;
    indTPIncl = ud.ad.Nskip+1:ntrs;
  end
  
  ud.tpx(nthrun).incl = indTPIncl;
  ud.tpx(nthrun).excl = indTPExcl;
  
end

return;
%----------------------------------------------------------%
function [ud, err] = load_rescale(ud)

err = 1;

for nthrun = 1:ud.ad.Nruns,
  if(~isempty(ud.XX.RescaleTarget) & ud.XX.RescaleTarget ~= 0)
    run = ud.ad.runlist(nthrun);
    mvf= sprintf('%s/%03d/%s.meanval',ud.fsd,run,ud.ad.funcstem);
    fid = fopen(mvf,'r');
    if(fid == -1) 
      fprintf('ERROR: loading meanval %s \n',mvf);
      return;
    end
    ud.meanval(nthrun) = fscanf(fid,'%f');
    ud.rescale(nthrun) = ud.XX.RescaleTarget/ud.meanval(nthrun);
  else
    ud.rescale(nthrun) = 1;
  end
end

err = 0;
return;

%-------------------------------------------------%
function setcurpostxt(ud)

cpstring = sprintf('t = %g, v = %g ',ud.curpos(1),ud.curpos(2));
set(ud.curpostxt,'string',cpstring);

return;

