function r = swapview(varargin)
% r = swapview(varargin)

version = '$Id: swapview.m,v 1.1 2003/07/31 02:34:00 greve Exp $';
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

  s.volsize = size(s.vol1);

  % Rescale Vol1 to be between 1 and (s.ncmap-1) %
  minbase = min(reshape1d(s.vol1));
  maxbase = max(reshape1d(s.vol1));
  s.vol1 = floor((s.ncmap-2)*(s.vol1-minbase)/(maxbase-minbase)) + 1;

  % Rescale Vol2 to be between 1 and (s.ncmap-1) %
  minbase = min(reshape1d(s.vol2));
  maxbase = max(reshape1d(s.vol2));
  s.vol2 = floor((s.ncmap-2)*(s.vol2-minbase)/(maxbase-minbase)) + 1;
  
  % Create the color map
  s.cmap = gray(s.ncmap); % Start with gray

  s.curvox(1) = round(size(s.vol1,1)/2);
  s.curvox(2) = round(size(s.vol1,2)/2);
  s.curvox(3) = round(size(s.vol1,3)/2);

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
  %set(gcf,'WindowButtonMotionFcn','swapview(''wbm'');');
  set(gcf,'UserData',s);

  if(isempty(s.title)) s.title = 'SwapView'; end
  set(gcf,'Name',s.title);

  set(gcf,'Interruptible','off'); % smooth animations %
  set(gcf,'DoubleBuffer','on');   % smooth animations %
  set(gcf,'BusyAction','cancel'); % dont build up a lot of events %
  set(gcf,'renderer','painters'); % seems to be the best
  
  fprintf('\n');
  fprintf(' ----------------------------------\n');
  fprintf(' For help press "h"\n');
  fprintf(' ------------ Help ----------------\n');
  printhelp;
  fprintf('\n');
  fprintf(' ------ Current State -------------\n');
  printstate(s);
  fprintf('\n');

  % force an event to get things rolling %
  swapview('r');
  
  return;
end
%---------------------------------------------------------%
%---------------------------------------------------------%
%---------------------------------------------------------%

%----------- Parse the call-back function ----------%
s = get(gcf,'UserData');
figure(s.hfig);
redraw = 0;
switch(flag)
  
 case {'r'}, % refresh
  redraw = 1;
   
 case {'wbm'} % -------Window Button Move ------------ %
  return;
  
 case {'wbd'} % -------Window Button Down ------------ %
  
  xyz = get(gca,'CurrentPoint');
  c = round(xyz(1,1));
  r = round(xyz(1,2));
  if(r < 1 | r > size(s.displayimg,1) | c < 1 | c > size(s.displayimg,2))
    return;
  end
  s.curpoint = [r c];
  if(~s.mosview)
    s.curvox(1:2) = round([r c]);
  else
    [rv cv sv] = mossub2volsub(r, c, size(s.vol1));
    s.curvox(1:3) = [rv cv sv];
  end
  if(s.verbose)
    fprintf('Current Point: ');
    fprintf('%2d ',s.curpoint);
    fprintf('\n');
  end
  redraw = 1;
  
 case {'kbd'} %----------- Keyboard -------------%
  c = get(s.hfig,'CurrentCharacter'); 
  %fprintf('c = %s (%d)\n',c,c);
  switch(c)
   
   case {'t'},
    if(s.curvol == 1) s.curvol = 2;
    else              s.curvol = 1;
    end
    redraw = 1;
   
   case {'m'},
    s.mosview = ~s.mosview;
    if(s.mosview) 
      s.displayimg1 = vol2mos(s.vol1(:,:,:,s.curvox(4)));
      s.displayimg2 = vol2mos(s.vol2(:,:,:,s.curvox(4)));
      [r c] = volsub2mossub(s.curvox(1),s.curvox(2),s.curvox(3), ...
				 size(s.vol1));
      s.curpoint = [r c];
    else
      s.displayimg1 = vol2mos(s.vol1(:,:,s.curvox(3),s.curvox(4)));
      s.displayimg2 = vol2mos(s.vol2(:,:,s.curvox(3),s.curvox(4)));
      s.curpoint = [s.curvox(1) s.curvox(2)];
    end
    redraw = 2;
    
   case {30}, % up arrow
    if(s.mosview) return; end
    s.curvox(3) = s.curvox(3) + 1;
    if(s.curvox(3) > size(s.vol1,3)) s.curvox(3) = 1; end
    s.displayimg1 = s.vol1(:,:,s.curvox(3),s.curvox(4));
    s.displayimg2 = s.vol2(:,:,s.curvox(3),s.curvox(4));
    redraw = 1;
    
   case {31}, % down arrow
    if(s.mosview) return; end
    s.curvox(3) = s.curvox(3) - 1;
    if(s.curvox(3) < 1) s.curvox(3) = size(s.vol1,3); end
    s.displayimg1 = s.vol1(:,:,s.curvox(3),s.curvox(4));
    s.displayimg2 = s.vol2(:,:,s.curvox(3),s.curvox(4));
    redraw = 1;
    
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
    fprintf('\n');
    fprintf('\n');
    printstate(s);
    fprintf('\n');
    fprintf('\n');
   
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

if(redraw > 0)
  figure(s.hfig);
  axes(s.haxis);
  if(s.curvol == 1) s.displayimg = s.displayimg1;
  else              s.displayimg = s.displayimg2;
  end
  if(redraw == 1)
    set(s.himage,'CData',s.displayimg);
  else
    s.himage = image(s.displayimg);    
    axis image;
    s.hMarker = [];
  end
  if(s.MarkerOn)
    hold on;
    if(~isempty(s.hMarker)) delete(s.hMarker); end
    s.hMarker = plot(s.curpoint(2),s.curpoint(1),'g+');
    %ud.hMarkerRow = plot(ud.CurPixel(2),1,'gv');
    %ud.hMarkerCol = plot(1,ud.CurPixel(1),'g>');
    hold off;
  end
  drawnow;
end

set(gcf,'UserData',s);

return;
%--------------------------------------------------%
%--------------------------------------------------%
%--------------------------------------------------%

%--------------------------------------------------%
%% Print Usage 
function print_usage(dummy)
  fprintf(1,'USAGE:\n');
  fprintf(1,'  swapview\n');
  fprintf(1,'     -init \n');
  fprintf(1,'     -v1  volume1  \n');
  fprintf(1,'     -v2  volume2\n');
  fprintf(1,'     -title  title\n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = main_struct
  s.MarkerOn       = 1;
  s.hMarker        = [];
  s.volsize        = [];
  s.mossize        = [];
  s.vol1           = [];
  s.vol2           = [];
  s.hfig           = [];
  s.haxis          = [];
  s.himage         = [];
  s.displayimg     = [];
  s.displayimg1    = [];
  s.displayimg2    = [];
  s.curpoint       = [1 1]; % current point in the display img
  s.curvox         = [1 1 1 1]; % current point [r c s f]
  s.curvol         = 1; % current volume
  s.curview        = 1; % current plane view 1 = rc
  s.mosview        = 0; % 1 = view mos, 0 = view vol  
  s.zoomstate      = 0;
  s.ncmap          = 64;
  s.verbose        = 0;
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

      case {'-title'},
        arg1check(flag,narg,ninputargs);
        s.title = inputargs{narg};
        narg = narg + 1;

      case '-verbose',
        s.verbose = 1;

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

if( isempty(s.vol1) )
  fprintf(2,'ERROR: No volume 1 specified\n');
  s=[]; return;
end

if( isempty(s.vol2) )
  fprintf(2,'ERROR: No volume 2 specified\n');
  s=[]; return;
end

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
return;

%--------------------------------------------------%
function printhelp
  fprintf('b - change brush size (1x1,3x3,5x5,7x7,9x9,11x11,13x13,...)\n');
  fprintf('h - print help (this message)\n');
  fprintf('m - change edit mode\n');
  fprintf('    No Edit - button down does nothing\n');
  fprintf('    Set     - button down turns on mask at voxel\n');
  fprintf('    Unset   - button down turns off mask at voxel\n');
  fprintf('    Toggle  - button down toggels mask at voxel\n');
  fprintf('q - quit/exit\n');
  fprintf('s - save mask\n');
  fprintf('t - toggle mask on and off\n');
  fprintf('v - print the current state\n');
  fprintf('z - toggle zoom state. When zoom is on, use left\n');
  fprintf('    button to zoom in and right to zoom out\n');
return;
