function r = swapview(varargin)
% r = swapview(varargin)


%
% swapview.m
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

version = 'swapview.m @FS_VERSION@';
r = 1;

%% Print usage if there are no arguments %%
if(nargin == 0)
  print_usage;
  return;
end

hcurrentfig = get(0,'CurrentFigure');

%--------------- Initialize if needed -----------------%
init = 0;
flag = deblank(varargin{1});
if(~isempty(strmatch(flag,'-init')) | isempty(hcurrentfig))
  init = 1;
  s = parse_args(varargin);
  if(isempty(s)) return; end
  s = check_params(s);
  if(isempty(s)) return; end

  if(~isempty(s.stem1))
    fprintf('Loading %s\n',s.stem1);
    s.vol1 = fast_ldbslice(s.stem1,-1);
    if(isempty(s.vol1))
      fprintf('ERROR: could not load %s\n',s.stem1);
      return;
    end
    s.vol1 = s.vol1(:,:,:,1); % first frame
  end
  
  if(~isempty(s.stem2))
    fprintf('Loading %s\n',s.stem2);
    s.vol2 = fast_ldbslice(s.stem2,-1);
    if(isempty(s.vol2))
      fprintf('ERROR: could not load %s\n',s.stem2);
      return;
    end
    s.vol2 = s.vol2(:,:,:,1); % first frame
  end
  
  % Rescale Vol1 to be between 1 and (s.ncmap-1) %
  minbase = min(reshape1d(s.vol1));
  maxbase = max(reshape1d(s.vol1));
  if(s.rescale)
    s.vol1 = floor((s.ncmap-2)*(s.vol1-minbase)/(maxbase-minbase)) + 1;
  end
  
  % Rescale Vol2 to be between 1 and (s.ncmap-1) %
  minbase = min(reshape1d(s.vol2));
  maxbase = max(reshape1d(s.vol2));
  if(s.rescale)
    s.vol2 = floor((s.ncmap-2)*(s.vol2-minbase)/(maxbase-minbase)) + 1;
  end
  
  % Create the color map
  s.cmap = gray(s.ncmap); % Start with gray

  s.curvox(1) = round(size(s.vol1,1)/2);
  s.curvox(2) = round(size(s.vol1,2)/2);
  s.curvox(3) = round(size(s.vol1,3)/2);
  s.prevvox = s.curvox;

  s.displayimg1 = s.vol1(:,:,s.curvox(3));
  s.displayimg2 = s.vol2(:,:,s.curvox(3));
  s.displayimg = s.displayimg1;
  s.curpoint(1) = s.curvox(1);
  s.curpoint(2) = s.curvox(2);
  
  % Image the display image %
  figure;
  s.himage = image(s.displayimg1);
  set(s.himage,'EraseMode','none');
  axis image;
  s.hfig = gcf;
  s.haxis = gca;
  colormap(s.cmap);
  set(gcf,'pointer','crosshair');

  % Set up the call-back functions %
  set(gcf,'KeyPressFcn',          'swapview(''kbd'');');
  set(gcf,'WindowButtonDownFcn',  'swapview(''wbd'');');
  %set(gcf,'WindowButtonUpFcn',    'swapview(''wbu'');');
  set(gcf,'WindowButtonMotionFcn','swapview(''wbm'');');

  if(isempty(s.title)) s.title = 'SwapView'; end
  set(gcf,'Name',s.title);

  set(gcf,'Interruptible','off'); % smooth animations %
  set(gcf,'DoubleBuffer','on');   % smooth animations %
  set(gcf,'BusyAction','cancel'); % dont build up a lot of events %
  set(gcf,'renderer','painters'); % seems to be the best
  
  % Set up a menus %
  s.viewmenu = uimenu('Label','View');
  uimenu(s.viewmenu,'Label','Mosaic','Callback','swapview(''mosaic'');');
  uimenu(s.viewmenu,'Label','Slice', 'Callback','swapview(''slice'');');
  hc = get(s.viewmenu,'children');
  set(hc(1),'checked','on');
   
  % Set up a menus %
  s.curvolmenu = uimenu('Label','CurVol');
  uimenu(s.curvolmenu,'Label','Vol1','Callback','swapview(''vol1'');','Tag','CurVolMenuVol1');
  uimenu(s.curvolmenu,'Label','Vol2','Callback','swapview(''vol2'');','Tag','CurVolMenuVol2');
  %s.curvolmenu.vol2 = uimenu(s.curvolmenu,'Label','Vol2','Callback','swapview(''vol2'');');
  %s.curvolmenu.vol1 = uimenu(s.curvolmenu,'Label','Vol1','Callback','swapview(''vol1'');');
  %set(s.curvolmenu.vol1,'checked','on');
  h = findobj(s.curvolmenu,'Tag','CurVolMenuVol1');
  set(h,'checked','on');

  s.mosbut = uicontrol('Style', 'pushbutton', 'String', 'Mosaic','Position', ...
	    [1   1 50 50], 'Callback', 'swapview(''mostoggle'');');
  s.swapbut = uicontrol('Style', 'pushbutton', 'String', 'Swap','Position', ...
	    [1  50 50 50], 'Callback', 'swapview(''voltoggle'');');
  uicontrol('Style', 'pushbutton', 'String', 'State','Position', ...
	    [1 100 50 50], 'Callback', 'swapview(''state'');');
  s.curpostxt = uicontrol('Style', 'text','Position',  [1 150 60 20]);
  s.mousepostxt = uicontrol('Style', 'text','Position',[1 170 60 20]);
  s.vol2valtxt = uicontrol('Style', 'text','Position',  [1 190 60 20]);
  s.vol1valtxt = uicontrol('Style', 'text','Position',  [1 210 60 20]);
  s.curvoltxt = uicontrol('Style', 'text','Position',  [1 230 60 20]);

  nslices = size(s.vol1,3)
  if(nslices > 1) 
    d = 1/(nslices-1);
    s.sliceslider = uicontrol('Style','slider','Min',1,'Max',nslices,...
			      'SliderStep',[d 3*d],...
			      'value',s.curvox(3),...
			      'position', [1 260 20 120],...
			      'callback','swapview(''sliceslider'');');
  else 
    s.sliceslider = [];
  end
  
  s.uicontour = ...
      uicontrol('Style', 'checkbox', ...
		'String', 'Contour',...
		'Position', [60 1 75 25], ...
		'Callback', 'swapview(''contour_cb'');',...
		'tooltipstring','Toggle display of contour.',...
		'value',s.contour);
  
  fprintf('\n');
  fprintf(' ----------------------------------\n');
  fprintf(' For help press "h"\n');
  if(s.verbose)
    fprintf(' ------------ Help ----------------\n');
    printhelp;
    fprintf('\n');
    fprintf(' ------ Current State -------------\n');
    printstate(s);
    fprintf('\n');
  end

  if(~isempty(s.twfstemlist))
    s.htwf = rawplot;
    rawplot('AddVolumes',s.twfstemlist,s.htwf);
  end

  set(gcf,'UserData',s);

  % force an event to get things rolling %
  swapview('r');
  
  return;
end
%---------------------------------------------------------%
%---------------------------------------------------------%
%---------------------------------------------------------%

%----------- Parse the call-back function ----------%
s = get(gcf,'UserData');
if(~isfield(s,'swapview')) return; end
figure(s.hfig);
redraw = 0;
switch(flag)
  
 case {'mosaic'}
  if(s.mosview == 1) return; end
  if(size(s.vol1,3)==1)  return; end
  hc = get(s.viewmenu,'children');
  set(hc(1),'checked','off');
  set(hc(2),'checked','on');
  set(s.mosbut,'string','Slice');
  s.mosview = 1;
   
 case {'slice'}
  if(s.mosview ~= 1) return; end
  hc = get(s.viewmenu,'children');
  set(hc(1),'checked','on');
  set(hc(2),'checked','off');
  set(s.mosbut,'string','Mosaic');
  s.mosview = 0;
   
 case {'mostoggle'}
  s.mosview = ~s.mosview;
  if(s.mosview) swapview('mosaic');
  else          swapview('slice');
  end
  %fprintf('mosview = %d %d\n',s.mosview,s.mosviewprev); 
  return;
 
 case {'vol1'}
  if(s.curvol == 1) return; end
  h = findobj(s.curvolmenu,'Tag','CurVolMenuVol1');
  set(h,'checked','on');
  h = findobj(s.curvolmenu,'Tag','CurVolMenuVol2');
  set(h,'checked','off');
  s.curvol = 1;
   
 case {'vol2'}
  if(s.curvol == 2) return; end
  h = findobj(s.curvolmenu,'Tag','CurVolMenuVol1');
  set(h,'checked','off');
  h = findobj(s.curvolmenu,'Tag','CurVolMenuVol2');
  set(h,'checked','on');
  s.curvol = 2;
  
 case {'voltoggle'}
  if(s.curvol == 1) swapview('vol2');
  else              swapview('vol1');
  end
  return;
  
 case {'setslice'}
  if(~s.mosview) 
    s.displayimg1 = s.vol1(:,:,s.curvox(3),s.curvox(4));
    s.displayimg2 = s.vol2(:,:,s.curvox(3),s.curvox(4));
  else
    [r c] = volsub2mossub(s.curvox(1), s.curvox(2), s.curvox(3), size(s.vol1));
    s.curpoint = [r c];
  end
  redraw = 1;
   
 case{'upslice'}
  if(s.mosview) return; end
  s.curvox(3) = s.curvox(3) + 1;
  if(s.curvox(3) > size(s.vol1,3)) s.curvox(3) = 1; end
  set(gcf,'UserData',s);
  swapview('setslice');
  return;
  
 case{'downslice'}
  if(s.mosview) return; end
  s.curvox(3) = s.curvox(3) - 1;
  if(s.curvox(3) < 1) s.curvox(3) = size(s.vol1,3); end
  set(gcf,'UserData',s);
  swapview('setslice');
  return;
 
 case {'sliceslider'}
  v = round(get(s.sliceslider,'value'));
  %fprintf('v = %g\n',v);
  s.curvox(3) = v;
  set(gcf,'UserData',s);
  swapview('setslice');
  return;
   
 case {'contour_cb'}
  s.contour = ~s.contour;
  if(~s.contour & ishandle(s.hcontour)) 
    delete(s.hcontour); 
    s.hcontour = -1;
    set(gcf,'UserData',s);
    return;
  elseif(~s.mosview)
    hold on;
    if(ishandle(s.hcontour)) delete(s.hcontour); end
    [cm s.hcontour] = contour(s.vol1(:,:,s.curvox(3)),'r');
    hold off;
  end
  set(gcf,'UserData',s);
  drawnow;
  return;
   
 case{'state'}
  fprintf('\n\n');
  printstate(s);
  fprintf('\n\n');
  return;

 case {'r'}, % refresh
  redraw = 1;
   
 %case {'wbm'} % -------Window Button Move ------------ %
 
 case {'wbd','wbm'} % -------Window Button Down ------------ %
  
  xyz = get(gca,'CurrentPoint');
  c = round(xyz(1,1));
  r = round(xyz(1,2));
  if(r < 1 | r > size(s.displayimg,1) | c < 1 | c > size(s.displayimg,2))
    set(gcf,'pointer','arrow');
    return;
  end
  set(gcf,'pointer','crosshair');
  if(~s.mosview)
    curvox = [round([r c]) s.curvox(3) s.curvox(4) ];
  else
    [rv cv sv] = mossub2volsub(r, c, size(s.vol1));
    if(isempty(rv)) return; end % clicked zero padding
    curvox = [rv cv sv s.curvox(4)];
  end
  if(strcmp(flag,'wbm')) 
    curvoxstr = sprintf('%d %d %d',curvox(1),curvox(2),curvox(3));
    set(s.mousepostxt,'string',curvoxstr);
    return;
  end
  s.curpoint = [r c]; 
  s.curvox = curvox;
  if(s.verbose)
    fprintf('Current Point: ');
    fprintf('%2d ',s.curpoint);
    fprintf('\n');
  end
  if(~isempty(s.htwf)) 
    rawplot('SetCurPoint',s.curvox(1),s.curvox(2),s.curvox(3),s.htwf);
  end
  redraw = 1;
  
 case {'kbd'} %----------- Keyboard -------------%
  c = get(s.hfig,'CurrentCharacter'); 
  %fprintf('c = %s (%d)\n',c,c);
  switch(c)
   
   case {'t'},
    swapview('voltoggle');
    return;
   
   case {'m'},
    swapview('mostoggle');
    return;
    
   case {'u',30}, % up arrow
    swapview('upslice');
    return;
    
   case {'d',31}, % down arrow
    swapview('downslice');
    return;
    
   case {'h'},
    fprintf('\n');
    printhelp;
    fprintf('\n');
    printstate(s);
    fprintf('\n');
   
   case {'o'},
    s.MarkerOn = ~s.MarkerOn;
    if(s.MarkerOn == 0)
      if(~isempty(s.hMarker)) delete(s.hMarker); end
      s.hMarker = [];
    else
      hold on;
      s.hMarker = plot(s.curpoint(2),s.curpoint(1),'g+');
      %ud.hMarkerRow = plot(ud.CurPixel(2),1,'gv');
      %ud.hMarkerCol = plot(1,ud.CurPixel(1),'g>');
      hold off;
    end
    drawnow;
   
   case {'q'},
    close(s.hfig);
    r = 0;
    return;
   
   case {'v'},
    swapview('state');
    return;
   
   case {'z'},
    zoom;
    s.zoomstate = ~ s.zoomstate;
    set(gcf,'UserData',s);
    if(s.zoomstate)
      fprintf('Zoom is on\n');
    else
      fprintf('Zoom is off\n');
    end
    return;
  end
  
end % --- switch(flag) ----- %
%-------------------------------------------------%
%-------------->>>>>>><<<<<<<<<<<<<<--------------%

% Change from mos to slice or the other 
if(s.mosviewprev ~= s.mosview)
  if(s.mosview) 
    set(gcf,'pointer','watch');
    s.displayimg1 = vol2mos(s.vol1(:,:,:,s.curvox(4)));
    s.displayimg2 = vol2mos(s.vol2(:,:,:,s.curvox(4)));
    set(gcf,'pointer','crosshair');
    [r c] = volsub2mossub(s.curvox(1),s.curvox(2),s.curvox(3), ...
			  size(s.vol1));
    s.curpoint = [r c];
  else
    s.displayimg1 = vol2mos(s.vol1(:,:,s.curvox(3),s.curvox(4)));
    s.displayimg2 = vol2mos(s.vol2(:,:,s.curvox(3),s.curvox(4)));
    s.curpoint = [s.curvox(1) s.curvox(2)];
  end
  s.mosviewprev = s.mosview;
  redraw = 2;
end

% Change the current volume %
if(s.prevvol ~= s.curvol | redraw > 0)
  if(s.curvol == 1) s.displayimg = s.displayimg1;
  else              s.displayimg = s.displayimg2;
  end
  s.prevvol = s.curvol;
  redraw = max(redraw,1);
end

% redraw %
if(redraw > 0)
  if(s.verbose) fprintf('redraw = %d\n',redraw); end
  figure(s.hfig);
  axes(s.haxis);
  if(redraw == 1)
    set(s.himage,'CData',s.displayimg);
  else
    if(ishandle(s.hMarker)) delete(s.hMarker);  end
    s.himage = image(s.displayimg);    
    axis image;
    s.hMarker = [];
  end
  if(s.MarkerOn)
    hold on;
    %if(~isempty(s.hMarker)) delete(s.hMarker);  end
    if(ishandle(s.hMarker)) delete(s.hMarker);  end
    s.hMarker = plot(s.curpoint(2),s.curpoint(1),'g+');
    %ud.hMarkerRow = plot(ud.CurPixel(2),1,'gv');
    %ud.hMarkerCol = plot(1,ud.CurPixel(1),'g>');
    hold off;
  end
  if(~s.mosview & s.contour)
    hold on;
    if(ishandle(s.hcontour)) delete(s.hcontour); end
    [cm s.hcontour] = contour(s.vol1(:,:,s.curvox(3)),'r');
    hold off;
  end

  drawnow;
end

curvoxstr = sprintf('%d %d %d',s.curvox(1),s.curvox(2),s.curvox(3));
set(s.curpostxt,'string',curvoxstr);
if(s.curvol == 1)  set(s.curvoltxt,'string','Vol1');
else               set(s.curvoltxt,'string','Vol2');
end
vol1valstr = sprintf('Vol1: %g',s.vol1(s.curvox(1),s.curvox(2),s.curvox(3)));
set(s.vol1valtxt,'string',vol1valstr);
vol2valstr = sprintf('Vol2: %g',s.vol2(s.curvox(1),s.curvox(2),s.curvox(3)));
set(s.vol2valtxt,'string',vol2valstr);

set(gcf,'UserData',s);

return;
%--------------->>>>>>><<<<<<<<<<<<<<--------------%
%--------------------------------------------------%
%-- end main --------------------------------------%
%--------------------------------------------------%

%--------------------------------------------------%
%% Print Usage 
function print_usage(dummy)
  fprintf(1,'USAGE:\n');
  fprintf(1,'  swapview\n');
  fprintf(1,'     -init \n');
  fprintf(1,'     -v1  volume1  \n');
  fprintf(1,'     -v2  volume2\n');
  fprintf(1,'     -stem1  volume1  \n');
  fprintf(1,'     -stem2  volume2\n');
  fprintf(1,'     -title  title\n');
  fprintf(1,'     -norescale\n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = main_struct
  s.swapview       = 1;
  s.MarkerOn       = 1;
  s.hMarker        = [];
  s.vol1           = [];
  s.vol2           = [];
  s.stem1           = [];
  s.stem2           = [];
  s.hfig           = [];
  s.haxis          = [];
  s.himage         = [];
  s.displayimg     = [];
  s.displayimg1    = [];
  s.displayimg2    = [];
  s.contour        = 0;
  s.hcontour       = [];
  s.curpoint       = [1 1]; % current point in the display img
  s.curvox         = [1 1 1 1]; % current vox index [r c s f]
  s.prevvox        = [1 1 1 1]; % prev vox index [r c s f]
  s.curvol         = 1; % current view volume
  s.prevvol        = 1; % previous current view volume
  s.curview        = 1; % current plane view 1 = rc
  s.mosview        = 0; % 1 = view mos, 0 = view vol  
  s.mosviewprev    = 0; % previous mosview
  s.zoomstate      = 0;
  s.ncmap          = 64;
  s.twfstemlist    = ''; % list of twf volumes to plot
  s.htwf           = [];
  s.verbose        = 0;
  s.rescale        = 1;
  s.title          = '';
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

      case {'-v1'},
        arg1check(flag,narg,ninputargs);
        s.vol1 = inputargs{narg};
        narg = narg + 1;

      case {'-v2'},
        arg1check(flag,narg,ninputargs);
        s.vol2 = inputargs{narg};
        narg = narg + 1;

      case {'-stem1'},
        arg1check(flag,narg,ninputargs);
        s.stem1 = inputargs{narg};
        narg = narg + 1;

      case {'-stem2'},
        arg1check(flag,narg,ninputargs);
        s.stem2 = inputargs{narg};
        narg = narg + 1;

      case {'-title'},
        arg1check(flag,narg,ninputargs);
        s.title = inputargs{narg};
        narg = narg + 1;

      case {'-twf'},
        arg1check(flag,narg,ninputargs);
        s.twfstemlist = strvcat(s.twfstemlist,inputargs{narg});
        narg = narg + 1;

      case '-verbose',
        s.verbose = 1;
     
     case '-rescale',
        s.rescale = 1;

     case '-norescale',
        s.rescale = 0;

      case {'-debug','-echo','-init'}, % ignore

      otherwise
        fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
        s = [];
        return;

    end % --- switch(flag) ----- %

  end % while(narg <= ninputargs)

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check that there is at least one more argument %%
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)
%fprintf(1,'Checking Parameters\n');

if( isempty(s.vol1) & isempty(s.stem1) )
  fprintf(2,'ERROR: No volume 1 specified\n');
  s=[]; return;
end

if( isempty(s.vol2) & isempty(s.stem2) )  
  fprintf(2,'ERROR: No volume 2 specified\n');
  s=[]; return;
end

if( isempty(s.vol1) & isempty(s.vol2)) return; end

for n = 1:length(size(s.vol1))
  if( size(s.vol1,n) ~= size(s.vol2,n) )
    fprintf(2,'ERROR: No dimension mismatch\n');
    size(s.vol1)
    size(s.vol2)
    s=[]; return;
  end
end
  
return;

%--------------------------------------------------%
function printstate(s)
  fprintf('Current Vol: %d\n',s.curvol);
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
  fprintf('Contour: %d\n',s.contour);
return;

%--------------------------------------------------%
function printhelp
  fprintf('m - toggle mosaic\n');
  fprintf('t - toggle volume\n');
  fprintf('o - toggle marker\n');
  fprintf('q - quit/exit\n');
  fprintf('v - print the current state\n');
  fprintf('z - toggle zoom state. When zoom is on, use left\n');
  fprintf('    button to zoom in and right to zoom out\n');
return;
