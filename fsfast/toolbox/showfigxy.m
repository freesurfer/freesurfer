function showfigxy(varargin)
% Displays the last clicked point and current mouse point
% in the current figure. If there is an image in the window,
% the value at the point is also printed out.
%
% To use, create the figure (eg, using plot or image), then, with that
% figure as the current figure, run showfigxy;
%
% Run showfigxy('help') for key press commands or hit 'h'
% when the cursor is in the figure to print help to terminal.
%
% To change the crosshair color, run 
%    showfigxy('crosshaircolor','newcolor') 
% where newcolor is a single letter color spec accepted by plot.
%
% Can also run:
%    showfigxy('init',vol) 
% To get values from vol instead of figure image
%
% Note: this will take over keyboard and mousing callbacks!
%
%
%


%
% showfigxy.m
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

if(nargin == 0)  event = 'init'; 
else             event = varargin{1};
end

if(strcmp(event,'help'))
  printhelp;
  return;
end

if(~strcmp(event,'init'))  hfig = gcf; end
if(strcmp(event,'init') & nargin ~= 2)
  hfig = gcf; 
end
hfig = gcf; 

%--------------------------------------------------------------%
switch(event)
 
 case 'init'
  set(hfig,'KeyPressFcn','showfigxy(''kbd'')');
  set(hfig,'WindowButtonDownFcn','showfigxy(''wbd'')');
  set(hfig,'WindowButtonMotionFcn','showfigxy(''wbm'')');
  ud = get(hfig,'UserData');
  ud.showfigxy = 1;
  ud.curxytxt = uicontrol('Style', 'text','Position',  [1 1 250 30]);
  ud.mvxytxt = uicontrol('Style', 'text','Position',  [260 1 250 30]);
  ax = axis;
  ud.curxy = [(ax(1)+ax(2))/2 (ax(3)+ax(4))/2];
  ud.mvxy  = [(ax(1)+ax(2))/2 (ax(3)+ax(4))/2];
  ud.hcursor1 = [];
  ud.hcursor2 = [];
  ud.zoom  = 0;
  ud.xzoom = 0;
  ud.yzoom = 0;
  ud.usecrosshair = 1;
  ud.crosshaircolor = 'g';
  ud.v = [];
  if(nargin == 2) ud.v = varargin{2}; end
  ud = drawcrosshair(ud);
  set(hfig,'UserData',ud);
  
 case 'crosshaircolor'
  if(nargin == 2) 
    ud = get(hfig,'UserData');
    if(~isshowfigxy(ud)) return; end
    ud.crosshaircolor = varargin{2};
    ud = drawcrosshair(ud); 
    set(hfig,'UserData',ud);
  end
  
 case 'wbd';
  ud = get(hfig,'UserData');
  if(~isshowfigxy(ud)) return; end
  xyz = get(gca,'CurrentPoint');
  x = xyz(1,1);
  y = xyz(1,2);  
  ax = axis;
  if(x < ax(1) | x > ax(2) | y < ax(3) | y > ax(4)) return; end
  ud.curxy = [x y];
  set(hfig,'UserData',ud);
  setxystring(hfig,'cur');
  %fprintf('x = %g, y = %g\n',x,y);
  ud = drawcrosshair(ud);
  set(hfig,'UserData',ud);
  
 case 'wbm';
  ud = get(hfig,'UserData');
  if(~isshowfigxy(ud)) return; end
  xyz = get(gca,'CurrentPoint');
  x = xyz(1,1);
  y = xyz(1,2);  
  ax = axis;
  if(x < ax(1) | x > ax(2) | y < ax(3) | y > ax(4)) return; end
  ud.mvxy = [x y];
  set(hfig,'UserData',ud);
  setxystring(hfig,'mv');
  %fprintf('x = %g, y = %g\n',x,y);
  
 case 'kbd';

  ud = get(hfig,'UserData');
  if(~isshowfigxy(ud)) return; end
  c = get(hfig,'CurrentCharacter'); 
  %fprintf('c = %s (%d)\n',c,c);

  switch(c)
   
   case 'c', 
    ud.usecrosshair = ~ud.usecrosshair;
    ud = drawcrosshair(ud);
     
   case 'd', 
    db = get(hfig,'DoubleBuffer');
    if(strcmp(db,'on'))
      set(hfig,'DoubleBuffer','off');
    else
      set(hfig,'DoubleBuffer','on');
    end     
   
   case 'h', 
    printhelp;
   
   case 'x', 
    ud.xzoom = ~ud.xzoom ;
    if(ud.xzoom) zoom xon;
    else         zoom off;
    end
   
   case 'y', zoom;
    ud.yzoom = ~ud.yzoom ;
    if(ud.yzoom) zoom yon;
    else         zoom off;
    end
   
   case 'z', zoom;
    ud.zoom = ~ud.zoom ;
    if(ud.zoom)   
      ud.xzoom = 1; 
      ud.yzoom = 1; 
      zoom on;
    else
      ud.xzoom = 0; 
      ud.yzoom = 0; 
      zoom off;
    end
  
  end

  set(hfig,'UserData',ud);
end

return;  
%---------------------------------------------------------%

%---------------------------------------------------------%
function setxystring(hfig,type)
  ud = get(hfig,'UserData');
  if(~isempty(ud.v)) img = ud.v;
  else               img = isimage;
  end
  
  if(strcmp(type,'cur'))
    x = ud.curxy(1); 
    y = ud.curxy(2);
    htxt = ud.curxytxt;
  else
    x = ud.mvxy(1); 
    y = ud.mvxy(2);
    htxt = ud.mvxytxt;
  end
  r = round(y);
  c = round(x);

  xystring = sprintf('x = %g, y = %g ',x,y);
  if(~isempty(img) & r > 0 & r <= size(img,1) & c > 0 & c <= size(img,2) )
    v = img(round(y),round(x));
    xystring = sprintf('%s, v = %g',xystring,v);
  end

  set(htxt,'string',xystring);
return;  

%---------------------------------------------------------%
function img = isimage
  % Image if the current axis has an image in it.
  % Bug: not correct for bar graphs
  img = [];
  
  chlist = get(gca,'children');
  nc = length(chlist);
  for n = 1:nc
    c = chlist(n);
    s = get(c);
    if(isfield(s,'CData')) 
      img = s.CData; 
      return; 
    end
  end

return;

%---------------------------------------------------------%
function ud = drawcrosshair(ud)
  if(ud.usecrosshair)
    if(~isempty(ud.hcursor1) & ishandle(ud.hcursor1))
      delete(ud.hcursor1);
      delete(ud.hcursor2);
    end
    a = axis;
    x = ud.curxy(1);
    y = ud.curxy(2);
    hold on;
    ud.hcursor1 = plot([a(1) a(2)],[y y],ud.crosshaircolor);
    ud.hcursor2 = plot([x x],[a(3) a(4)],ud.crosshaircolor);
    hold off;
  else
    if(~isempty(ud.hcursor1) & ishandle(ud.hcursor1))
      delete(ud.hcursor1);
      delete(ud.hcursor2);
    end
  end
   
return;

%---------------------------------------------------------%
function r = isshowfigxy(ud)
  if(~isfield(ud,'showfigxy')) r = 0;
  else r = 1;
  end
return;

%---------------------------------------------------------%
function printhelp
  fprintf('\n');
  help showfigxy;
  fprintf('c - toggle crosshair\n');
  fprintf('d - toggle double buffering\n');
  fprintf('h - print help\n');
  fprintf('x - toggle x zoom\n');
  fprintf('y - toggle y zoom\n');
  fprintf('z - toggle zoom\n');
  fprintf('\n');

return;


