function showfigxy(event,hfig)
% Displays the last clicked point and current mouse point
% in the current figure. To use, create the figure 
% (eg, using plot or image), then, with that figure as 
% the current figure, run showfigxy;
%
% $Id: showfigxy.m,v 1.1 2004/01/14 22:01:08 greve Exp $
%

if(nargin == 0)  event = 'init'; end
if(nargin < 2)   hfig = gcf; end

switch(event)
 
 case 'init'
  set(hfig,'WindowButtonDownFcn','showfigxy(''wbd'')');
  set(hfig,'WindowButtonMotionFcn','showfigxy(''wbm'')');
  ud = get(hfig,'UserData');
  ud.curxytxt = uicontrol('Style', 'text','Position',  [1 1 200 20]);
  ud.curxy = [0 0];
  ud.mvxytxt = uicontrol('Style', 'text','Position',  [210 1 200 20]);
  ud.mvxy = [0 0];
  ud.hcursor1 = [];
  ud.hcursor2 = [];
  set(hfig,'UserData',ud);
  
  case 'wbd';
   ud = get(hfig,'UserData');
   xyz = get(gca,'CurrentPoint');
   x = xyz(1,1);
   y = xyz(1,2);  
   ax = axis;
   if(x < ax(1) | x > ax(2) | y < ax(3) | y > ax(4)) return; end
   ud.curxy = [x y];
   set(hfig,'UserData',ud);
   setxystring(hfig,'cur');
   %fprintf('x = %g, y = %g\n',x,y);
   a = axis;
   if(~isempty(ud.hcursor1) & ishandle(ud.hcursor1))
     delete(ud.hcursor1);
     delete(ud.hcursor2);
   end
   hold on;
   ud.hcursor1 = plot([a(1) a(2)],[y y],'g-');
   ud.hcursor2 = plot([x x],[a(3) a(4)],'g-');
   hold off;
   set(hfig,'UserData',ud);
   
   
  case 'wbm';
   ud = get(hfig,'UserData');
   xyz = get(gca,'CurrentPoint');
   x = xyz(1,1);
   y = xyz(1,2);  
   %ax = axis;
   %if(x < ax(1) | x > ax(2) | y < ax(3) | y > ax(4)) return; end
   ud.mvxy = [x y];
   set(hfig,'UserData',ud);
   setxystring(hfig,'mv');
   %fprintf('x = %g, y = %g\n',x,y);

end

return;  
%---------------------------------------------------------%

%---------------------------------------------------------%
function setxystring(hfig,type)
  ud = get(hfig,'UserData');
  if(strcmp(type,'cur'))
    xystring = sprintf('x = %g, y = %g ',ud.curxy(1),ud.curxy(2));
    set(ud.curxytxt,'string',xystring);
  else
    xystring = sprintf('x = %g, y = %g ',ud.mvxy(1),ud.mvxy(2));
    set(ud.mvxytxt,'string',xystring);
  end
return;  
