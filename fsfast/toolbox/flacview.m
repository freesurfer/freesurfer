function r = flacview(varargin)
% r = flacview(varargin)


%
% flacview.m
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

version = 'flacview.m @FS_VERSION@';
r = 1;

%% Print usage if there are no arguments %%
if(nargin == 0)
  print_usage;
  return;
end

hcurrentfig = get(0,'CurrentFigure');
init = 0;
flag = deblank(varargin{1});
%fprintf('flag1 = %s\n',flag);

%--------- Initialize ----------------------------%
if(~strcmp(flag,'-cb') | isempty(hcurrentfig))
  init = 1;

  % Parse the args
  s = parse_args(varargin);
  if(isempty(s)) return; end
  s = check_params(s);
  if(isempty(s)) return; end

  % Load the flac
  s.flac = fast_ldflac(s.flacfile);
  if(isempty(s.flac))
    fprintf('ERROR: loading %s\n',s.flacfile);
    return;
  end
  
  % Customize
  s.flac.sess = s.sesspath;
  s.flac.nthrun = s.nthrun;
  s = load_data(s);
  if(isempty(s)) return; end  
  s.nruns = size(s.flac.runlist,1);
  
  % Prepare the overlay
  s = make_overlay(s);
  
  if(s.mosview)
    s.displayimg = vol2mos(s.displayvol);
  else
    s.displayimg = squeeze(s.displayvol(:,:,s.curvox(3),:));
  end
  
  % Set up the image window
  figure;
  s.himage = image(s.displayimg);
  set(s.himage,'EraseMode','none');
  axis image;
  s.hfig = gcf;
  s.haxis = gca;
  set(gcf,'pointer','crosshair');

  set(gcf,'KeyPressFcn',          'flacview(''-cb'',''kbd'');');
  set(gcf,'WindowButtonDownFcn',  'flacview(''-cb'',''wbd'');');
  %set(gcf,'WindowButtonMotionFcn','flacview(''-cb'',''wbm'');');
  set(gcf,'Name',s.title);

  set(gcf,'Interruptible','off'); % smooth animations %
  set(gcf,'DoubleBuffer','on');   % smooth animations %
  set(gcf,'BusyAction','cancel'); % dont build up a lot of events %
  set(gcf,'renderer','painters'); % seems to be the best

  % Set up the plot window
  s.hplot = figure;
  set(gcf,'Interruptible','off'); % smooth animations %
  set(gcf,'DoubleBuffer','on');   % smooth animations %
  set(gcf,'BusyAction','cancel'); % dont build up a lot of events %
  set(gcf,'renderer','painters'); % seems to be the best
  s.t = s.flac.TR*[0:s.flac.ntp-1]';
  s = make_timecourse(s);

  % Return control to the image window
  figure(s.hfig);

  % Set up the run menu
  s.runmenu = uimenu('Label','Run');
  for nthrun = 1:s.nruns
    cbf = sprintf('flacview(''-cb'',''Run'',''%d'');',nthrun);
    fprintf('%s\n',cbf);
    uimenu(s.runmenu,'Label',s.flac.runlist(nthrun,:),'Callback',cbf);
  end
  hc = get(s.runmenu,'children');
  set(hc(end),'checked','on');
  
  % Set up the contrast menu
  s.conmenu = uimenu('Label','Contrast');
  ncon = length(s.flac.con);
  for nthcon = ncon:-1:1
    conname = s.flac.con(nthcon).name;
    cbf = sprintf('flacview(''-cb'',''contrast'',''%d'');',nthcon);
    %fprintf('%s\n',cbf);
    uimenu(s.conmenu,'Label',conname,'Callback',cbf);
  end
  hc = get(s.conmenu,'children');
  set(hc(1),'checked','on');
  
  s.plotmenu = uimenu('Label','Plot');
  uimenu(s.plotmenu,'Label','Raw','Callback',...
	 'flacview(''-cb'',''plot'',''Raw'');');  
  uimenu(s.plotmenu,'Label','Task','Callback',...
	 'flacview(''-cb'',''plot'',''Task'');');  
  uimenu(s.plotmenu,'Label','Resid','Callback',...
	 'flacview(''-cb'',''plot'',''Resid'');');  
  uimenu(s.plotmenu,'Label','Nuis','Callback',...
	 'flacview(''-cb'',''plot'',''Nuis'');');  
  uimenu(s.plotmenu,'Label','NuisResid','Callback',...
	 'flacview(''-cb'',''plot'',''NuisResid'');');  
  uimenu(s.plotmenu,'Label','PMF','Callback',...
	 'flacview(''-cb'',''plot'',''PMF'');');  
  uimenu(s.plotmenu,'Label','PMFResid','Callback',...
	 'flacview(''-cb'',''plot'',''PMFResid'');');  
  uimenu(s.plotmenu,'Label','ACF','Callback',...
	 'flacview(''-cb'',''plot'',''ACF'');');  
  plotmenu_check(s);

  % Set the user data
  set(gcf,'UserData',s);

  % force an event to get things rolling %
  flacview('-cb','redraw');

  return;
end

% ------ Only get's here if the first arg is -cb ---------%

%----------- Parse the call-back function ----------%
s = get(gcf,'UserData');

% Exit if this is not a flacview window
if(~isfield(s,'flacview')) return; end

figure(s.hfig);
redraw = 0;
flag = deblank(varargin{2});
%fprintf('flag2 = %s\n',flag);
switch(flag)
 
 case {'contrast'},
  nthcon = sscanf(deblank(varargin{3}),'%d');
  %fprintf('nthcon = %d, %d\n',nthcon,s.nthcon);
  if(s.nthcon == nthcon) return; end
  hc = get(s.conmenu,'children');
  set(hc(s.nthcon),'checked','off');
  set(hc(nthcon),'checked','on');
  s.nthcon = nthcon;
  s = make_overlay(s);
  s = make_timecourse(s);
  redraw = 1;  

 case {'plot'},
  plottype = deblank(varargin{3});
  %fprintf('plottype = %s, %s\n',plottype,s.plottype);
  if(strcmp(plottype,s.plottype)) return; end
  s.plottype = plottype;
  plotmenu_check(s);
  s = make_timecourse(s);

 case {'Run'},
  nthrun = sscanf(deblank(varargin{3}),'%d');
  if(s.flac.nthrun == nthrun) return; end
  s.flac.nthrun = nthrun;
  s = load_data(s);
  if(isempty(s)) return; end  
  s = make_overlay(s);
  s = make_timecourse(s);
  set(gcf,'Name',s.title);
  redraw = 1;    

 case {'redraw'},
  redraw = 1;

 %case {'wbd','wbm'}
 case {'wbd'}
  xyz = get(gca,'CurrentPoint');
  c = round(xyz(1,1));
  r = round(xyz(1,2));
  sz = size(s.displayimg);
  if(r < 1 | r > sz(1) | c < 1 | c > sz(2) )
    set(gcf,'pointer','arrow');
    return;
  end
  set(gcf,'pointer','crosshair');
  if(s.mosview)
    [rv cv sv] = mossub2volsub(r, c, s.beta.volsize);
  else
    rv = r;
    cv = c;
    sv = s.curvox(3);
  end
  if(strcmp(flag,'wbd'))
    s.curvox(1) = rv;
    s.curvox(2) = cv;
    s.curvox(3) = sv;
  end
  %fprintf('r = %d, c = %d, s = %d\n',rv,cv,sv);
  s = make_timecourse(s);
  
 case {'kbd'} %----------- Keyboard -------------%
  c = get(s.hfig,'CurrentCharacter'); 
  %fprintf('c = %s (%d)\n',c,c);
  switch(c)
   case {'m'}
    s.mosviewprev = s.mosview;
    s.mosview = ~s.mosview;
    redraw = 1;
   case {'o'}
    tit  = 'Adjust Overlay Threshold';
    prompt = {'Sat Threshold:','Min Threshold'};
    lines  = 1;
    def = {sprintf('%7.4f',s.thsat),sprintf('%7.4f',s.thmin)};
    answer   = inputdlg(prompt,tit,lines,def);
    if(~isempty(answer))
      asat = sscanf(answer{1},'%f');
      amin = sscanf(answer{2},'%f');
      if(asat <= amin)
	msg = sprintf('Sat threshold (%f) cannot be less than min (%f)',...
		      asat,amin);
	errordlg(msg);
      elseif(asat ~= s.thsat | amin ~= s.thmin)
	s.thmin = amin;
	s.thsat = asat;
	s = make_overlay(s);	 
	redraw  = 1;
      end 
    end
  end % switch c
  
end % switch flag

%fprintf('redraw = %d\n',redraw);
if(redraw > 0)
  figure(s.hfig);
  axes(s.haxis);
  if(s.mosview)
    s.displayimg = vol2mos(s.displayvol);
  else
    s.displayimg = squeeze(s.displayvol(:,:,s.curvox(3),:));
  end
  if(s.mosviewprev == s.mosview)
    set(s.himage,'CData',s.displayimg);
  else
    s.himage = image(s.displayimg);    
    axis image;
  end
  tit = sprintf('Contrast: %s',s.flac.con(s.nthcon).name);
  title(tit);
end

set(gcf,'UserData',s);

return;

%--------------->>>>>>><<<<<<<<<<<<<<--------------%
%--------------------------------------------------%
%-- end main --------------------------------------%
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = main_struct
  s.flacview       = 1;
  s.sesspath       = '';
  s.nthrun         = 1;
  s.flac           = [];
  s.func           = [];
  s.mask           = [];
  s.acfseg         = [];
  s.beta           = [];
  s.con            = [];
  s.thmin          = 2;
  s.thsat          = 5;
  s.nthcon         =  1;
  s.hfig           = [];
  s.haxis          = [];
  s.himage         = [];
  s.MarkerOn       = 1;
  s.hMarker        = [];
  s.displayimg     = [];
  s.displayimg1    = [];
  s.displayimg2    = [];
  s.ncmap          = 64;
  s.curpoint       = [1 1]; % current point in the display img
  s.curvox         = [1 1 1 1]; % current vox index [r c s f]
  s.prevvox        = [1 1 1 1]; % prev vox index [r c s f]
  s.curvol         = 1; % current view volume
  s.prevvol        = 1; % previous current view volume
  s.curview        = 1; % current plane view 1 = rc
  s.mosview        = 1; % 1 = view mos, 0 = view vol  
  s.mosviewprev    = 1; % previous mosview
  s.sliceslider    = [];
  s.zoomstate      = 0;
  s.verbose        = 0;
  s.title          = '';
  s.hplot          = [];
  s.plottype       = 'Raw';
return;

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = main_struct;

  inputargs = varargin{1};
  ninputargs = length(inputargs);

  narg = 1;
  while(narg <= ninputargs)

    flag = deblank(inputargs{narg});
    narg = narg + 1;
    %fprintf(1,'Argument: %s\n',flag);
    if(~isstr(flag))
      flag
      fprintf(1,'ERROR: All Arguments must be a string\n');
      error;
    end

    switch(flag)

      case {'-f'},
        arg1check(flag,narg,ninputargs);
        s.flacfile = inputargs{narg};
        narg = narg + 1;

      case {'-s'},
        arg1check(flag,narg,ninputargs);
        s.sesspath = inputargs{narg};
        narg = narg + 1;

      case {'-r'},
        arg1check(flag,narg,ninputargs);
        s.nthrun = inputargs{narg};
        narg = narg + 1;

      case {'-t'},
        arg1check(flag,narg,ninputargs);
        s.title = inputargs{narg};
        narg = narg + 1;

      case '-verbose',
        s.verbose = 1;
     
      case {'-debug','-echo'}, % ignore

      otherwise
        fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
        s = [];
        return;

    end % --- switch(flag) ----- %

  end % while(narg <= ninputargs)

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

if(isempty(s.flacfile))
  fprintf(2,'ERROR: No flac file specified\n');
  s=[]; return;
end

if(isempty(s.sesspath))
  fprintf(2,'ERROR: No session specified\n');
  s=[]; return;
end

return;

%--------------------------------------------------%
%% Print Usage 
function print_usage(dummy)
  fprintf(1,'USAGE:\n');
  fprintf(1,'  flacview\n');
  fprintf(1,'     -f flacfile \n');
  fprintf(1,'     -s sesspath  \n');
  fprintf(1,'     -r nthrun\n');
  fprintf(1,'     -t title\n');
return
%--------------------------------------------------%

%--------------------------------------------------%
function printstate(s)
  fprintf('Sess:   %s\n',s.flac.sess);
  fprintf('NthRun: %d\n',s.flac.nthrun);
  fprintf('Current Voxel: ');
  fprintf('%2d ',s.curvox);
  fprintf('\n');
  fprintf('Current Point: ');
  fprintf('%2d ',s.curpoint);
  fprintf('\n');
  fprintf('Current Plane View: %d\n',s.curview);
  fprintf('Mos View: %d\n',s.mosview);
  if(~isempty(s.sliceslider))
    fprintf('SliceSlider: %g\n',get(s.sliceslider,'value'));
  end
return;

%--------------------------------------------------%
%% Check that there is at least one more argument %%
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;

%--------------------------------------------------%
function s = make_overlay(s)

basemos = vol2mos(s.base);
conmos = vol2mos(s.con(s.nthcon).fsig.vol);
[displaymos s.cmap s.cscale ] = ...
    imgoverlaytc2(basemos,conmos,s.thmin,s.thsat,'abs',0);
s.displayvol = mos2vol(displaymos,s.func.volsize);

return;

%--------------------------------------------------%
function s = make_timecourse(s)

switch(s.plottype)
 case 'Raw'
  s.y = squeeze(s.func.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
 case 'Task'
  b = squeeze(s.beta.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));  
  taskregind = flac_taskregind(s.flac);
  s.y = s.flac.X(:,taskregind) * b(taskregind);
 case 'Resid'
  b = squeeze(s.beta.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
  y = squeeze(s.func.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
  s.y = y - s.flac.X * b;
 case 'PMF'
  b = squeeze(s.beta.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
  C = s.flac.con(s.nthcon).C;
  s.y = s.flac.X*(C'*C*b);
 case 'PMFResid'
  b = squeeze(s.beta.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
  C = s.flac.con(s.nthcon).C;
  y = squeeze(s.func.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
  s.y = y - s.flac.X*(C'*C*b);
 case 'ACF'
  b = squeeze(s.beta.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
  y = squeeze(s.func.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
  r = y - s.flac.X * b;
  s.y = fast_acorr(r);
 case 'Nuis'
  b = squeeze(s.beta.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));  
  nuisregind = flac_nuisregind(s.flac);
  s.y = s.flac.X(:,nuisregind) * b(nuisregind);
 case 'NuisResid'
  y = squeeze(s.func.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));
  b = squeeze(s.beta.vol(s.curvox(1),s.curvox(2),s.curvox(3),:));  
  nuisregind = flac_nuisregind(s.flac);
  s.y = y - s.flac.X(:,nuisregind) * b(nuisregind);
end
  
figure(s.hplot)
plot(s.t,s.y);
tit = sprintf('%s (%d,%d,%d)\n',s.plottype,...
	      s.curvox(1),s.curvox(2),s.curvox(3));
title(tit);
xlabel('time (sec)');
figure(s.hfig);

return;

%--------------------------------------------------%
function plotmenu_check(s)

hc = get(s.plotmenu,'children');
nhc = length(hc);
for nthhc = 1:nhc
  label = get(hc(nthhc),'Label');
  if(strcmp(label,s.plottype))
    set(hc(nthhc),'checked','on');
  else
    set(hc(nthhc),'checked','off');
  end
end

return;

%--------------------------------------------------%
function s = load_data(s)

fprintf('Loading nthrun = %d\n',s.flac.nthrun);

  s.flac = flac_customize(s.flac);
  if(isempty(s.flac)) s = []; return; end
  s.flac = flac_desmat(s.flac);
  if(isempty(s.flac)) s = []; return; end

  s.title = sprintf('%s %s %d (%s)',s.flac.sess,s.flac.name,...
		    s.flac.nthrun,s.flac.runlist(s.flac.nthrun,:));

  % Load the mask
  fprintf('  Loading mask\n');
  mstem = sprintf('%s/%s/masks/%s',s.flac.sess,s.flac.fsd,s.flac.mask);
  s.mask = MRIread(mstem);
  if(isempty(s.mask)) s = []; return; end

  % Load beta
  fprintf('  Loading beta\n');
  bstem = sprintf('%s/%s/%s/%s/beta',s.flac.sess,s.flac.fsd,s.flac.name,...
		  s.flac.runlist(s.flac.nthrun,:));
  s.beta = MRIread(bstem);
  if(isempty(s.beta)) s = []; return; end

  % Load contrasts
  ncon = length(s.flac.con);
  for nthcon = 1:ncon
    fprintf('  Loading contrast %s\n',s.flac.con(nthcon).name);
    cstem = sprintf('%s/%s/%s/%s/%s/fsig',s.flac.sess,s.flac.fsd,...
		    s.flac.name,s.flac.runlist(s.flac.nthrun,:),...
		    s.flac.con(nthcon).name);
    s.con(nthcon).fsig = MRIread(cstem);
    if(isempty(s.con(nthcon).fsig)) s = []; return; end
  end
  
  % Load functional
  fprintf('  Loading functional\n');
  fstem = sprintf('%s/%s/%s/%s',s.flac.sess,s.flac.fsd,...
		  s.flac.runlist(s.flac.nthrun,:),s.flac.funcstem);
  s.func = MRIread(fstem);
  if(isempty(s.func)) s = []; return; end

  % Prepare the base
  base = s.beta.vol(:,:,:,1);
  basemin = min(base(:));
  basemax = max(base(:));
  base = round((s.ncmap-1)*(base-basemin)/(basemax-basemin))+1;
  cmgray = gray(s.ncmap);
  s.base = reshape(cmgray(base,:),[size(base) 3]);

  fprintf('Done loading data\n');
  
return;
