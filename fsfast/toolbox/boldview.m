function boldview(cbflag,data)


%
% boldview.m
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

global imgvol;

if(strcmp(cbflag,'subject'))
  ud = new_boldview_userdata;
  ud.subject = data

  %ud.cor = tstimg;
  %ud.cor = fmri_ldcor('dougg');
  %imgvol = fmri_ldcor('dougg');
  %ud.subject = 'nouchine';
  imgvol = fmri_ldcor(ud.subject);
  volmin = min(reshape1d(imgvol));
  volmax = max(reshape1d(imgvol));
  imgvol = uint8((imgvol - volmin) * 63/(volmax-volmin) + 1);

  screensize = get(0,'ScreenSize');
  %w = floor(.8*min(screensize(3),screensize(4)));
  w = 518;

  h = figure;
  set(gcf,'Interruptible','Off');
  set(gcf,'BusyAction','Cancel');
  cp = get(h,'Position');
  set(h,'Position',[cp(1) cp(2) w w]);
  ud.bvfig = h;

  ud.hsag = subplot(2,2,2);
  set(gca,'Units','Pixels');
  set(gca,'Position',[258 258 256 256]);
  axis image; set(gca,'XTick',[]); set(gca,'YTick',[]);

  ud.hcor = subplot(2,2,1);
  set(gca,'Position',[0 258 256 256]);
  set(gca,'Units','Pixels');
  axis image; set(gca,'XTick',[]); set(gca,'YTick',[]);

  ud.haxial = subplot(2,2,3);
  set(gca,'Position',[0 0 256 256]);
  set(gca,'Units','Pixels');
  axis image; set(gca,'XTick',[]); set(gca,'YTick',[]);

  ud.hlast = ud.hcor;

  set(gcf,'WindowButtonDownFcn','boldview(''wbd'')');
  set(gcf,'KeyPressFcn',        'boldview(''kbd'')');
  set(gcf,'Pointer','crosshair');

  ud = redraw(ud);
  set(gcf,'UserData',ud);
  return;
end

%---------------------------------------------------------------%
ud = get(gcf,'UserData');
switch (cbflag)

  case {'wbd'}
     xyz = get(gca,'CurrentPoint');
     x = xyz(1,1);
     y = xyz(1,2);
     if(x < 1 | x > 256 | y < 1 | y > 256) return; end
     fprintf('WBD: x=%g, y=%g  ',x,y);
     switch(gca)
       case ud.hcor, 
         fprintf('Corronal Slice\n');
         ud.curss(1) = floor(x); 
         ud.curss(3) = floor(y); 

       case ud.hsag, 
         fprintf('Sagital Slice\n');
         ud.curss(2) = floor(x); 
         ud.curss(3) = floor(y); 

       case ud.haxial, 
         fprintf('Axial Slice\n');
         ud.curss(1) = floor(x); 
         ud.curss(2) = 257-floor(y); 
     end
     ud.hlast = gca;
     ud = redraw(ud);

  case 'kbd', 
    c = get(gcf,'CurrentCharacter'); 
    fprintf(1,'Key %s\n',c);
    switch(c)
      case 'q', 
        close(gcf); 
        return;
      case 'z', 
        ud.zoom = ~ud.zoom;
        zoom;
      case {'+','='}
        if(ud.hlast == ud.hcor)       n = 2;
        elseif(ud.hlast == ud.hsag)   n = 1;
        elseif(ud.hlast == ud.haxial) n = 3;
        else return; 
        end
        if(c == '+') d = 10;
        else         d =  1;
        end
        ud.curss(n) = min(ud.curss(n)+d,256);
        ud = redraw(ud);
        fprintf('n = %d, d= %2d, (%3d,%3d,%3d)\n',n,d,...
                 ud.curss(1),ud.curss(2),ud.curss(3));
      case {'-','_'}
        if(ud.hlast == ud.hcor)       n = 2;
        elseif(ud.hlast == ud.hsag)   n = 1;
        elseif(ud.hlast == ud.haxial) n = 3;
        else return;
        end
        if(c == '_') d = 10;
        else         d =  1;
        end
        ud.curss(n) = max(ud.curss(n)-d,1);
        ud = redraw(ud);
        fprintf('n = %d, d= %2d, (%3d,%3d,%3d)\n',n,d,...
                 ud.curss(1),ud.curss(2),ud.curss(3));
      end

end

set(gcf,'UserData',ud);
return;

%----------------------------------------------------%
function ud = redraw(ud)
  global imgvol;

  sagimg = squeeze(imgvol(ud.curss(1),:,:))'; %'
  corimg = squeeze(imgvol(:,ud.curss(2),:))'; %'
  axiimg = flipud(squeeze(imgvol(:,:,ud.curss(3)))'); %'
     %sagimg = squeeze(ud.cor(ud.curss(1),:,:))'; %'
     %corimg = squeeze(ud.cor(:,ud.curss(2),:))'; %'
     %axiimg = flipud(squeeze(ud.cor(:,:,ud.curss(3)))'); %'

     axes(ud.hsag);
     if(~ud.zoom) zoom(1); end;
     set(gca,'Interruptible','Off');
     set(gca,'BusyAction','Cancel');
     %set(gca,'DrawMode','Fast');
     image(sagimg); colormap(gray);
     set(gca,'Units','Pixels');
     set(gca,'Position',[258 258 256 256]);
     set(gca,'NextPlot','Add');
     plot(ud.curss(2),ud.curss(3),'r+');
     set(gca,'NextPlot','Replace');
     axis image; set(gca,'XTick',[]); set(gca,'YTick',[]);

     axes(ud.hcor);
     if(~ud.zoom) zoom(1); end;
     set(gca,'Interruptible','Off');
     set(gca,'BusyAction','Cancel');
     set(gca,'DrawMode','Fast');
     image(corimg); colormap(gray);
     set(gca,'Units','Pixels');
     set(gca,'Position',[0 258 256 256]);
     set(gca,'NextPlot','Add');
     plot(ud.curss(1),ud.curss(3),'r+');
     set(gca,'NextPlot','Replace');
     axis image; set(gca,'XTick',[]); set(gca,'YTick',[]);

     axes(ud.haxial);
     if(~ud.zoom) zoom(1); end;
     set(gca,'Interruptible','Off');
     set(gca,'BusyAction','Cancel');
     set(gca,'DrawMode','Fast');
     image(axiimg); colormap(gray);
     set(gca,'Units','Pixels');
     set(gca,'Position',[0 0 256 256]);
     set(gca,'NextPlot','Add');
     plot(ud.curss(1),257-ud.curss(2),'r+'); 
     set(gca,'NextPlot','Replace');
     axis image; set(gca,'XTick',[]); set(gca,'YTick',[]);
return;

%----------------------------------------------------%
function ud = new_boldview_userdata
  ud.bvfig = [];
  ud.curss = [128 128 128];
  ud.hsag  = [];
  ud.hcor  = [];
  ud.haxial  = [];
  ud.hlast   = [];
  ud.subject = '';
  ud.cor    = []; % Corronal   volume %
  ud.fun    = []; % Functional volume %
  ud.Qcor   = []; % Subscript in cor given mm in cor %
  ud.Qfun   = []; % Subscript in fun given mm in fun %
  ud.Sfc    = []; % Subscript in fun given subscript in cor %
  ud.M12    = []; % mm in fun given mm in cor (register.dat) %
  ud.zoom   = 0;

return;
%----------------------------------------------------%
function img = tstimg

img = zeros(256,256,256);

img( 96:160, 64:192,   64:192) = 64;
img( 96:160, 192:200,  70:80)  = 64;
img(112:144, 64:192, 112:144) = 0;

return;
