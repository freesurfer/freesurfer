function yak(varargin)
% a1,underlay,pmin)
%
% yak displays images.  It can display a significance
% map over a structural.  Multiple images can be viewed
% by hitting the '+' and '-' keys.  The sign of the tail
% can be toggled by hitting 't'.  Text in the figure gives
% the current status of the image, including number of 
% rows, columns, and depth as well as the location and 
% pvalue of the mouse.  It is assumed that the pvalues
% are ln(p).
% 
% yak('init',img)
% yak('init',img,'-A',ActImg)
% yak('init',img,'-h',HsaImg) 
%
% yak(cbstring) % for callback functions
%
%


%
% yak.m
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

if(nargin == 0)
  msg = 'USAGE: hfig = yak(flag,options)';
  qoe(msg);error(msg);
end

Init = 0;
narg = 1;
while(narg <= nargin)
  cbflag = varargin{narg};

  if(~isstr(cbflag))
    msg = sprintf('Argument %d must be a string',narg);
    qoe(msg);error(msg);
  end
  
  narg = narg + 1;

  switch(cbflag)
    case 'init',
      arg1check(cbflag,narg,nargin);
      Init = 1;
      ud = new_user_data;
      ud.UnderlayImg = varargin{narg};
      basemean = mean(reshape1d(ud.UnderlayImg));
      imask = find(ud.UnderlayImg > .5*basemean);
      ud.MaskImg = zeros(size(ud.UnderlayImg));
      ud.MaskImg(imask) = 1;
      ud.MaskImg = repmat(ud.MaskImg, [1 1 3]);
    case '-A'
      arg1check(cbflag,narg,nargin);
      ud.ActImg = varargin{narg};
    case '-Nnnc'
      arg1check(cbflag,narg,nargin);
      ud.Nnnc = varargin{narg};
    case '-nskip'
      arg1check(cbflag,narg,nargin);
      ud.Nskip = varargin{narg};
    case '-Ch'
      arg1check(cbflag,narg,nargin);
      ud.Ch = varargin{narg};
    case '-hdat'
      arg1check(cbflag,narg,nargin);
      ud.hdat = varargin{narg};
    case '-dof'
      arg1check(cbflag,narg,nargin);
      ud.hdof = varargin{narg};
    case '-h'
      arg1check(cbflag,narg,nargin);
      if(~isempty(varargin{narg}));
          [ud.havg ud.hstd ud.Nh ] = ...
              hsa_convert(varargin{narg},ud.Nnnc);
      end
    case '-raw'
      arg1check(cbflag,narg,nargin);
      if(~isempty(varargin{narg}))
        ud.raw = varargin{narg};
        ud.yraw = squeeze(ud.raw(1,1,:));
        ud.traw = [1:length(ud.yraw)];
      end
    case '-rawfunc'
      arg1check(cbflag,narg,nargin);
      if(~isempty(varargin{narg}))
        ud.rawfunc = varargin{narg};
      end
    case '-off'
      arg1check(cbflag,narg,nargin);
      if(~isempty(varargin{narg}))
        ud.hoffset = varargin{narg};
      end
    case '-thresh'
    case '-pmin'
      arg1check(cbflag,narg,nargin);
      ud.PMin = varargin{narg};
    case '-maxthresh'
    case '-pmax'
      arg1check(cbflag,narg,nargin);
      ud.PMax = varargin{narg};
    case '-TR'
      arg1check(cbflag,narg,nargin);
      ud.TR = varargin{narg};
      if(ud.TR == 0) ud.TR = 1; end
      fprintf('Set TR = %g\n',ud.TR);
    case '-tPreStim'
      arg1check(cbflag,narg,nargin);
      ud.tPreStim = varargin{narg};
    otherwise
      ud = get(gcf,'UserData'); 
      if( ~isfield(ud,'hThisYakFig')) return; end
      handle_cb(cbflag);

  end % --- switch(cbflag) ----- %

  narg = narg + 1;

end %--------- while(narg <= nargin) ---------- %

if(Init)
  ncmap = 64;
  cmgray = gray(ncmap);
  fprintf('Rescaling Intensity Image\n');
  t1min = min(reshape1d(ud.UnderlayImg));
  t1max = max(reshape1d(ud.UnderlayImg));
  if(t1max ~= t1min)
    ud.UnderlayImg = round((ncmap-1)*(ud.UnderlayImg-t1min)/(t1max-t1min))+1;
    ud.UnderlayImgTC = cmgray(ud.UnderlayImg,:);
    ud.UnderlayImgTC = reshape(ud.UnderlayImgTC , [size(ud.UnderlayImg) 3]);
  else
    ud.UnderlayImgTC = ones([size(ud.UnderlayImg) 3]);
  end


  if(~isempty(ud.ActImg))
    % ud.ovmin = -log10(ud.PMin);
    % ud.ovmax = -log10(ud.PMax);
    ud.ovmin = ud.PMin;
    ud.ovmax = ud.PMax;
    fprintf('Computing Overlay ');
    %[ud.DisplayImg ud.CMap ud.CScale ud.ovmax] = ...
    %    imgoverlay(ud.UnderlayImg,ud.ActImg,ud.ovmin,ud.ovmax);
    [ud.DisplayImg ud.CMap ud.CScale ] = ...
        imgoverlaytc2(ud.UnderlayImgTC,...
          ud.ActImg(:,:,ud.CurPlaneNo),ud.ovmin,ud.ovmax,ud.tail,ud.interp);
    fprintf('... Done\n');
    % ud.PMax = 10^(-abs(ud.ovmax));
    fprintf('PMin = %4.2f, PMax = %4.2f\n',ud.PMin,ud.PMax);
    %image(ud.DisplayImg(:,:,ud.CurPlaneNo));
    image(ud.DisplayImg);
    colormap(ud.CMap);
    ud.hcbar = colorbar;
    set_colorbar_scale(ud.hcbar,ud.CMap,ud.CScale);
  else %% There is no overlay %%
    ud.DisplayImg = ud.UnderlayImg;
    colormap('gray');
    ud.CMap = colormap;
    imagesc(ud.DisplayImg(:,:,ud.CurPlaneNo));
    colormap(gray);
    ud.hcbar = colorbar;
  end
  set(gca,'FontUnits','points');
  set(gca,'FontSize',10);
  %set(gcf,'Name',ud.Title);
  ud.hThisYakFig = gcf;

  set(gcf,'UserData',ud);
  ud.hcbar = colorbar;

  set(gcf,'WindowButtonDownFcn','yak(''wbd'')');
  %set(gcf,'WindowButtonUpFcn','yak(''wbu'')');
  %set(gcf,'WindowButtonMotionFcn','yak(''wbm'')');
  set(gcf,'KeyPressFcn',        'yak(''kbd'')');
  set(gcf,'DeleteFcn',          'yak(''delete'')');

  set(gcf,'Interruptible','off'); % smooth animations %
  set(gcf,'DoubleBuffer','on');   % smooth animations %
  set(gcf,'BusyAction','cancel'); % dont build up a lot of events %
  set(gcf,'renderer','painters'); % seems to be the best
  
  s = sprintf('%5d, %3s, %4.2f, %4.2f',...
              ud.CurPlaneNo,ud.OverlaySign,ud.PMin,ud.PMax);

  xr =  get(gca,'XLim');
  yr =  get(gca,'YLim');
  h = text(xr(2),yr(1),s);
  set(h,'Tag','txStatus');
  set(h,'VerticalAlignment','Bottom');
  set(h,'HorizontalAlignment','Right');
  set(h,'FontUnits','points');
  set(h,'FontSize',10);
  set(h,'EraseMode','background');

  h = text(xr(1),yr(1),'-*-*-*-');
  set(h,'Tag','txMouse');
  set(h,'VerticalAlignment','Bottom');
  set(h,'FontUnits','points');
  set(h,'FontSize',10);
  set(h,'EraseMode','background');

  unix('date');

  return;
end % Initialize %

return;



%%%%%%%%%%%%%% Callback Routines %%%%%%%%%%%%%%%%%%
function handle_cb(cbflag)
ud = get(gcf,'UserData'); 
switch (cbflag)
  case 'delete', 
    if(ishandle2(ud.hHDR)) close(ud.hHDR); end
    close(gcf);
    if(exist('QuitOnError')) quit; end

  case {'wbm','wbu'}, return;

  case {'wbd'}
    %switch (cbflag)
    %  case 'wbm', if( ~ ud.DragOn ) return; end
    %  case 'wbd', ud.DragOn = 1;     %fprintf('DragOn\n');
    %  case 'wbu', ud.DragOn = 0;     %fprintf('DragOff\n');
    %end

    xyz = get(gca,'CurrentPoint');
    x = xyz(1,1);
    y = xyz(1,2);
    c = round(x);
    r = round(y);

    [nrows ncols nplanes] = size(ud.DisplayImg);
    if (r<1 | r > nrows | c < 1 | c > ncols) return; end;

    if(~isempty(ud.ActImg))
      [rf cf] = sind2find(r,c,ud);
      v =  ud.ActImg(rf,cf,ud.CurPlaneNo);
      % v =  sign(v)*10^(-abs(v));
      s = sprintf('r = %3d, c = %3d, (%2d,%2d), ovl = %6.4f',r,c,rf,cf,v);
    elseif(ishandle2(ud.hHDR))
      [rf cf] = sind2find(r,c,ud);
      s = sprintf('r = %3d, c = %3d, (%2d,%2d)',r,c,rf,cf);
    else
      v =  ud.UnderlayImg(r,c,ud.CurPlaneNo);
      s = sprintf('r = %3d, c = %3d, v = %g',r,c,v);
    end
    ud.CurPixel = [r c];

    if(ud.MarkerOn)
      hold on;
      if(isempty(ud.hMarker))
        ud.hMarker = plot(ud.CurPixel(2),ud.CurPixel(1),'g+');
        ud.hMarkerRow = plot(ud.CurPixel(2),1,'gv');
        ud.hMarkerCol = plot(1,ud.CurPixel(1),'g>');
      else
        set(ud.hMarker,'xdata',ud.CurPixel(2));
        set(ud.hMarker,'ydata',ud.CurPixel(1));
        set(ud.hMarkerRow,'xdata',ud.CurPixel(2));
        set(ud.hMarkerRow,'ydata',1);
        set(ud.hMarkerCol,'xdata',1);
        set(ud.hMarkerCol,'ydata',ud.CurPixel(1));
      end
      hold off;
    end

    fprintf('CurrentPlane: %2d, %s\n',ud.CurPlaneNo,s);

    % mosaic talaiarch stuff %
    ud.MosSlice = -1;
    figtitle = get(gcf,'Name');
    mtch0 = findstr(figtitle,'00-15');
    if(~isempty(mtch0)) ud.MosSlice = 0;
    else
      mtch1 = findstr(figtitle,'16-31');
      if(~isempty(mtch1)) ud.MosSlice = 1;end
    end

    if(ud.MosSlice > -1) 
      nmosimg = size(ud.UnderlayImg,1)/4;
      mosrow  = floor(r/nmosimg);
      moscol  = floor(c/nmosimg);
      imgnum  = 4*mosrow + moscol;
      slicenum = 16 * ud.MosSlice + imgnum;
      slicerow = rem(r-1,nmosimg) + 1;
      slicecol = rem(c-1,nmosimg) + 1 ;
      res = 256/nmosimg;
      yres = 8;
      xoff = 256/2  + res/2;
      yoff = -256/2 + yres/2;
      zoff = 256/2  + res/2;
      xtal =  -(xoff -  res * (slicecol));
      ytal =  yoff + yres * (slicenum);
      ztal =  zoff -  res * (slicerow);
    
      %fprintf('Mos Slice Num: %d\n',ud.MosSlice);
      %fprintf('Indices: num = %d, row = %d, col = %d\n',...
      %        slicenum,slicerow,slicecol);
      fprintf('Tal Coords (mm): x = %g, y = %g, z = %g \n',xtal,ytal,ztal);
    end


    h = findobj(gcf,'Tag','txMouse');
    set(h,'String',s);

    ud = redraw_hdr(ud);

    if(~isempty(ud.raw))
      szraw = size(ud.raw);
      szraw = szraw(1:2);
      szdisp = size(ud.DisplayImg);
      szdisp = szdisp(1:2);
      [rraw craw] = rcmap_eqfov(szdisp,r,c,szraw);
      ud.yraw = squeeze(ud.raw(rraw,craw,:));
      ny = length(ud.yraw);
      ind = [(ud.Nskip+1):ny];
      ud.yraw = ud.yraw(ind);
      ud.traw = reshape1d(ud.TR*(ind-1));
    end

    if(~ud.DragOn & ishandle2(ud.hRaw))
      hyak = gcf;
      figure(ud.hRaw);
  
    if(isempty(ud.rawfunc))
        plot(ud.traw,ud.yraw,'b+-');
        tmp = sprintf('Raw Time Courses (%d,%d)',rraw,craw);
        title(tmp);
        set(gca,'xtick',[0:5*ud.TR:(ny-1)*ud.TR]);
        set(gca,'xlim',[0 (ny-1)*ud.TR]);
      else
        cmd = sprintf('%s(''plot'',ud.traw,ud.yraw);',ud.rawfunc);
        %cmd = sprintf('%s(ud.traw,ud.yraw)',ud.rawfunc);
        eval(cmd);
      end
      figure(hyak);
    end

  case 'kbd', 
    c = get(gcf,'CurrentCharacter'); 
    %fprintf(1,'Key %s\n',c);

    switch(c)

     case {'+','='},
       ud.CurPlaneNo = ud.CurPlaneNo + 1;
       if(~isempty(ud.ActImg))
         if(ud.CurPlaneNo > size(ud.ActImg,3)) 
            ud.CurPlaneNo = 1;
         end
       else
         if(ud.CurPlaneNo > size(ud.UnderlayImg,3)) 
            ud.CurPlaneNo = 1;
         end
       end
       ud = redraw(ud);
       setstatus(ud);
       if(ishandle2(ud.hHDR))
         hyak = gcf;
         hdrviewlite('tmarker',ud.hHDR,ud.CurPlaneNo);
         figure(hyak);
       end

     case {'-','_'},
       ud.CurPlaneNo = ud.CurPlaneNo - 1;
       if(ud.CurPlaneNo < 1) 
         ud.CurPlaneNo = size(ud.ActImg,3);
       end
       ud = redraw(ud);
       setstatus(ud);
       if(ishandle2(ud.hHDR))
         hyak = gcf;
         hdrviewlite('tmarker',ud.hHDR,ud.CurPlaneNo);
         figure(hyak);
       end

    case 'h' % Help
      s = sprintf('Keypress Commands: \n'); 
      s = sprintf(' %s+ - step to the next image\n',s);
      s = sprintf(' %s- - step to the previous image\n',s);
      s = sprintf(' %sc - toggle on crosshair\n',s);
      s = sprintf(' %sd - apply FDR thresholding\n',s);
      s = sprintf(' %sf - goto a row/col\n',s);
      s = sprintf(' %sh - help (this message)\n',s);
      s = sprintf(' %si - toggle interpolation\n',s);
      s = sprintf(' %sg - view hemodynamic waveforms\n',s);
      s = sprintf(' %sk - mask\n',s);
      s = sprintf(' %sm - toggle marker\n',s);
      s = sprintf(' %so - adjust overlay thresholds\n',s);
      s = sprintf(' %sq - close and exit yakview \n',s);
      s = sprintf(' %sr - show raw time course \n',s);
      s = sprintf(' %st - change sign to view (+,-,+/-)\n',s);
      s = sprintf(' %sv - save raw timecourse to ascii file\n',s);
      s = sprintf(' %sw - toggle display of overlay\n',s);
      s = sprintf(' %sz - toggle zoom state\n',s);
      msgbox(s,'Yakview Help','help');

      case 'i', 
      ud.interp = ~ud.interp;
      if(~exist('imresize') & ud.interp)
        s = sprintf('Cannot Interpolate:\n no image processing toolbox');
        msgbox(s,'Interpolation');
        ud.interp = 0;
      else
        ud = redraw(ud);
      end
      setstatus(ud);

    case 'k' 
      ud.Mask = ~ud.Mask;
      ud = redraw(ud);
      setstatus(ud);

    case 'a' 
      ud.MultiRaw = ~ud.MultiRaw;
      setstatus(ud);

    case 'm' 
      ud.MarkerOn = ~ud.MarkerOn;
      if(~ud.MarkerOn) 
        ud.hMarker = []; 
        ud.hMarkerRow = []; 
        ud.hMarkerCol = []; 
      end
      ud = redraw(ud);

    case 'o' 
      if(~isempty(ud.ActImg))
        tit  = 'Adjust Overlay Threshold';
        prompt = {'Max Threshold:','Min Threshold'};
        lines  = 1;
        def = {sprintf('%7.4f',ud.PMax),sprintf('%7.4f',ud.PMin)};
        answer   = inputdlg(prompt,tit,lines,def);
        amax = sscanf(answer{1},'%f');
        amin = sscanf(answer{2},'%f');
        if(amax <= amin)
          msg = sprintf('Max threshold (%f) cannot be less than min (%f)',...
                         amax,amin);
          errordlg(msg);
        elseif(amax ~= ud.PMax | amin ~= ud.PMin)
          ud.PMin = amin;
          ud.PMax = amax;
          ud.ovmax = ud.PMax;
          ud.ovmin = ud.PMin;
          [ud.DisplayImg ud.CMap ud.CScale ] = ...
                 imgoverlaytc2(ud.UnderlayImgTC,ud.ActImg(:,:,ud.CurPlaneNo),...
                 ud.ovmin,ud.ovmax,ud.tail,ud.interp);
          ud = redraw(ud);
          setstatus(ud);
        end 
      else
        msg = sprintf('There is no overlay image for which to adjust thresholds.');
        errordlg(msg);
      end

         case 'd' 
      if(~isempty(ud.ActImg))
        tit  = 'False Discovery Rate';
        prompt = {'FDR:'};
        lines  = 1;
        def = {sprintf('%7.4f',ud.FDR)};
        answer   = inputdlg(prompt,tit,lines,def);
        ud.FDR = sscanf(answer{1},'%f');
	p = 10.^(-abs(ud.ActImg(:,:,ud.CurPlaneNo)));
	pthresh = fast_fdrthresh(p,ud.FDR);
	ud.PMin = -log10(pthresh);
	ud.PMax = ud.PMin + 3;
	
	ud.ovmax = ud.PMax;
	ud.ovmin = ud.PMin;
	[ud.DisplayImg ud.CMap ud.CScale ] = ...
	    imgoverlaytc2(ud.UnderlayImgTC,ud.ActImg(:,:,ud.CurPlaneNo),...
			  ud.ovmin,ud.ovmax,ud.tail,ud.interp);
	ud = redraw(ud);
	setstatus(ud);
      else
        msg = sprintf('There is no overlay image for which to adjust thresholds.');
        errordlg(msg);
      end

     case 'f' 
      tit  = 'Goto Row/Col';
      prompt = {'Row:','Col'};
      lines  = 1;
      def = {sprintf('%d',ud.CurPixel(1)),sprintf('%d',ud.CurPixel(2))};
      answer   = inputdlg(prompt,tit,lines,def);
      arow = sscanf(answer{1},'%d');
      acol = sscanf(answer{2},'%d');
      fprintf('new pos: %d %d\n',arow,acol);
      if(arow < 1 | arow > size(ud.DisplayImg,1) | ...
         acol < 1 | acol > size(ud.DisplayImg,2))
        msg = sprintf('Specified point (%d,%d) is out of bounds',arow,acol);
        errordlg(msg);
      else
        ud.GoToPixel = [arow acol];
        ud.CurPixel = [arow acol];
        ud = redraw(ud);
        setstatus(ud);
      end

    case 'c' 
      ud.CrossHairOn = ~ud.CrossHairOn;
      ud = redraw(ud);

    case 'p' % pink
      ncmap = 64;
      cm = pink(ncmap);
      ud.UnderlayImgTC = cm(ud.UnderlayImg,:);
      ud.UnderlayImgTC = reshape(ud.UnderlayImgTC , [size(ud.UnderlayImg) 3]);
      ud = redraw(ud);
      setstatus(ud);

    case 'b' % bone
      ncmap = 64;
      cm = bone(ncmap);
      ud.UnderlayImgTC = cm(ud.UnderlayImg,:);
      ud.UnderlayImgTC = reshape(ud.UnderlayImgTC , [size(ud.UnderlayImg) 3]);
      ud = redraw(ud);
      setstatus(ud);

     case {'t'},
       if(~isempty(ud.ActImg))
        switch(ud.OverlaySign)
          case '+/-', 
             ud.OverlaySign = '+';
             img   = ud.ActImg;
             ud.tail = 'pos';
          case '+', 
             ud.OverlaySign = '-';
             img   = -ud.ActImg;
             ud.tail = 'neg';
          case '-', 
             ud.OverlaySign = 'abs';
             img = abs(ud.ActImg);
             ud.tail = 'abs';
          case 'abs', 
             ud.OverlaySign = '+/-';
             img = abs(ud.ActImg);
             ud.tail = 'both';
          end % ------ switch(ud.OverlaySign) ----- %

          % ud.ovmax = -log10(ud.PMax);
          % ud.ovmin = -log10(ud.PMin);
          ud.ovmax = ud.PMax;
          ud.ovmin = ud.PMin;
          %[ud.DisplayImg ud.CMap ud.CScale] = ...
          %          imgoverlay(ud.UnderlayImg,img,ud.ovmin,ud.ovmax);

          [ud.DisplayImg ud.CMap ud.CScale ] = ...
               imgoverlaytc2(ud.UnderlayImgTC,ud.ActImg(:,:,ud.CurPlaneNo),...
               ud.ovmin,ud.ovmax,ud.tail,ud.interp);

          fprintf('Sign Changed to %s\n',ud.OverlaySign);
          ud = redraw(ud);
          setstatus(ud);
        end % if %

     case {'g'},
      if(~isempty(ud.havg) & ~ishandle2(ud.hHDR))
        hyak = gcf;
        ud.t = ud.TR*[0:ud.Nh-1] - ud.tPreStim;
        ud.hHDR = figure;
        hdrviewlite('init',ud.hHDR,ud.t,ud.hdof);
        hdrviewlite('tmarker',ud.hHDR,ud.CurPlaneNo);
        havg = squeeze(ud.havg(1,1,:,:));
        hstd = squeeze(ud.hstd(1,1,:,:));

        if(~isempty(ud.hoffset))
          hoffset = ud.hoffset(1,1);
          hdrviewlite('plot',ud.hHDR,havg,hstd,[1 1],hoffset);
        else
          hdrviewlite('plot',ud.hHDR,havg,hstd,[1 1]);
        end
        
        set(ud.hHDR,'Name','HDRView');
        figure(hyak);
      end

     case {'z'},
       zoom;
       ud.ZoomState = ~ud.ZoomState;
       fprintf('Toggling Zoom State to %d\n',ud.ZoomState);

     case {'w'},
       ud.ovshow = ~ud.ovshow;
       fprintf('Toggling Show Overlay to %d\n',ud.ovshow);
       ud = redraw(ud);
       setstatus(ud);

     case {'r'},
       if(~isempty(ud.raw) & ~ishandle2(ud.hRaw))
         hyak = gcf;
	 %ud.yraw = squeeze(ud.raw(1,1,:));
         %ind = [(ud.Nskip+1):length(ud.yraw)];
	 %ud.yraw = ud.yraw(ind);
         %ud.traw = ud.TR*(ind-1);
         ud.hRaw = figure;
         if(isempty(ud.rawfunc))
           plot(ud.traw,ud.yraw,'b+-');
           % tit('Raw Time Courses');
         else
  	   cmd = sprintf('%s(''init'',ud.traw,ud.yraw);',ud.rawfunc);
           eval(cmd);
         end
         figure(hyak);
      end

    case 'v' 
      if(isempty(ud.raw))
        msgbox('Error: No raw time course loaded','','help');
      else
        curdir = pwd;
        cd(ud.SvDir);
        [fname pname] = uiputfile(ud.SvFile,'Save Raw Time Course');
        if(fname ~= 0)
          ud.SvFile = fname;
          if(isempty(pname)) ud.SvDir  = '.';
          else               ud.SvDir  = pname;
          end
          svname = strcat(ud.SvDir,ud.SvFile);
          fid = fopen(svname,'w');
          fprintf(fid,'%g\n',ud.yraw);
          fclose(fid);
        end
      end

     case {'q'},
      hyak = gcf;      
      if(ishandle2(ud.hHDR)) close(ud.hHDR); end
      if(ishandle2(ud.hRaw)) close(ud.hRaw); end
      close(hyak);
      if(exist('QuitOnError')) quit; end
      return;
    end %------- switch(c) --------- %

    otherwise,
      msg = sprintf('flag %s unrecognized',cbflag);
      qoe(msg);error(msg);

 end %--------- switch (cbflag) ------------%

 set_text(gca,ud);
 set(gcf,'UserData',ud);
return;

%%%------------------------------------------%%%%
function ud = redraw(ud)
  %set(gco,'CData',ud.DisplayImg(:,:,ud.CurPlaneNo));
  % colorbar;

  if(~isempty(ud.ActImg))
    if(ud.ovshow) 
      ovminuse = ud.ovmin;
      ovmaxuse = ud.ovmax;
    else
      ovminuse = 10^10;
      ovmaxuse = 10^11;
    end
    [ud.DisplayImg ud.CMap ud.CScale ] = ...
	imgoverlaytc2(ud.UnderlayImgTC,ud.ActImg(:,:,ud.CurPlaneNo),...
		      ovminuse,ovmaxuse,ud.tail,ud.interp);
      
    if(ud.Mask) ud.DisplayImg = ud.DisplayImg.*ud.MaskImg; end
    image(ud.DisplayImg);
    colormap(ud.CMap);
    ud.hcbar = colorbar;
    set_colorbar_scale(ud.hcbar,ud.CMap,ud.CScale);
    %image(ud.DisplayImg(:,:,ud.CurPlaneNo));
  else
    tmpimg = ud.DisplayImg(:,:,ud.CurPlaneNo);
    if(ud.Mask) tmpimg = tmpimg.*ud.MaskImg; end
    imagesc(tmpimg);
    colormap(gray);
    colorbar;
  end

  if(ud.MarkerOn)
    hold on;
    ud.hMarker = plot(ud.CurPixel(2),ud.CurPixel(1),'g+');
    ud.hMarkerRow = plot(ud.CurPixel(2),1,'gv');
    ud.hMarkerCol = plot(1,ud.CurPixel(1),'g>');
    hold off;
  end

  if(ud.CrossHairOn)
    set(gcf,'pointer','fullcrosshair');
  else
    set(gcf,'pointer','crosshair');
  end

  ud = redraw_hdr(ud);

  set_text(gca,ud);
return

ud.hHDR 
ishandle2(ud.hHDR)
%%%------------------------------------------%%%%
function ud = redraw_hdr(ud)
 if(ishandle2(ud.hHDR))
   r = ud.CurPixel(1);
   c = ud.CurPixel(2);
   [rh ch] = sind2hind(r,c,ud);
   fprintf('HDR: r=%d, c=%d  rh=%d, ch=%d\n',r,c,rh,ch);
   hyak = gcf;
   havg = squeeze(ud.havg(rh,ch,:,:));
   hstd = squeeze(ud.hstd(rh,ch,:,:));
      if(~isempty(ud.hoffset))
        hoffset = ud.hoffset(rh,ch);
        hdrviewlite('plot',ud.hHDR,havg,hstd,[rh ch],hoffset);
      else
        hdrviewlite('plot',ud.hHDR,havg,hstd,[rh ch]);
      end
      yaktitle = get(hyak,'Name');
      hdrtitle = sprintf('HDRView - %s',yaktitle);
      set(gcf,'Name',hdrtitle);
      figure(hyak);
  end
return;

%%%------------------------------------------%%%%
function set_colorbar_scale(hcb,cmap,cscale)
  ncmap = size(cmap,1);
  nyticks = 10;
  dytick = (ncmap-1)/(nyticks-1);
  yticks = round([1:dytick:ncmap]);
  set(hcb,'YTick',yticks);
  yticknum   = cscale(yticks);
  yticklabel = splitstring(sprintf('%7.2f ',yticknum)) ;
  set(hcb,'YTickLabel',yticklabel);
return

%%%------------------------------------------%%%%
function ud = new_user_data
   ud = struct('hThisYakFig',     [],...
               'DisplayImg',    [], ...
               'CurPlaneNo',     1, ...
               'CurPixel',     [1 1], ...
               'GoToPixel',     [], ...
               'OverlaySign',   '+/-', ...
               'ActImg',      [], ...
               'MaskImg',     [], ...
               'Mask',         0, ...
               'UnderlayImg',   [], ...
               'UnderlayImgTC', [], ...
               'interp',        0, ...
               'PMin',         .01, ...
               'PMax',         .00001, ...
               'FDR',          .01, ...
               'ovmin',        log(.01),...
               'ovmax',        log(.00001),...
               'tail',         'both',...
               'ovshow',       1,...
               'hcbar',        [],... % colorbar
               'CMap',          [], ...
               'CScale',        [], ...
               'hHDR',          [], ...
               'hdat',        [], ...
               'havg',        [], ...
               'hstd',        [], ...
               'hoffset',     [],...
               'Ch',          [], ...
               'evar',        [], ...
               'hdof',        [], ...
               'Nnnc',        [], ...
               'Nh',           0, ...
               'TR',           1,...
               'raw',         [],...
               'rawfunc',     [],...
               'yraw',         [],...
               'traw',         [],...
               'Nskip',       0,...
               'ShowRaw',     0,...
               'MultiRaw',    0,...
               'hRaw',        [],...
               'SvFile', '', ...
               'SvDir', '.', ...
               'ZoomState',    0,...
               'MosSlice',     -1,...
               'tPreStim',     0,...
               'Title',     'Yakview',...
               'hMarker',    [], ...
               'hMarkerRow',    [], ...
               'hMarkerCol',    [], ...
               'MarkerOn',   1, ...
               'CrossHairOn',   0, ...
               'DragOn',        0);

return

%%%------------------------------------------%%%%
function [rf, cf] = sind2find(rs,cs,ud)
  [nrs ncs tmp] = size(ud.UnderlayImg);
  [nrf ncf tmp] = size(ud.ActImg);

  rf = ceil(nrf * rs / nrs);
  cf = ceil(ncf * cs / ncs);
return;
%%%------------------------------------------%%%%
function [rh, ch] = sind2hind(rs,cs,ud)

  [nrh nch tmp] = size(ud.havg);
  if(~isempty(ud.UnderlayImg))
    [nrs ncs tmp] = size(ud.UnderlayImg);
  else
    if(~isempty(ud.ActImg))
      [nrs ncs tmp] = size(ud.ActImg);
    else
      [nrs ncs tmp] = size(ud.havg);
    end
  end

  rh = ceil(nrh * rs / nrs);
  ch = ceil(nch * cs / ncs);
return;
%%%------------------------------------------%%%%
function setstatus(ud)
  s = sprintf('%2d, %3s, %4.2f, %4.2f',...
              ud.CurPlaneNo,ud.OverlaySign,ud.PMin,ud.PMax);
  h = findobj(gcf,'Tag','txStatus');
  set(h,'String',s);

  r = ud.CurPixel(1);
  c = ud.CurPixel(1);

    if(~isempty(ud.ActImg))
      [rf cf] = sind2find(r,c,ud);
      v =  ud.ActImg(rf,cf,ud.CurPlaneNo);
      % v =  sign(v)*10^(-abs(v));
      s = sprintf('r = %3d, c = %3d, (%2d,%2d), ovl = %6.4f',r,c,rf,cf,v);
    elseif(ishandle2(ud.hHDR))
      v =  ud.UnderlayImg(r,c,ud.CurPlaneNo);
      [rf cf] = sind2find(r,c,ud);
      s = sprintf('r = %3d, c = %3d, (%2d,%2d), ovl = %6.4f',r,c,rf,cf,v);
    else
      v =  ud.UnderlayImg(r,c,ud.CurPlaneNo);
      s = sprintf('r = %3d, c = %3d, v = %g',r,c,v);
    end

return
%%%------------------------------------------%%%%
function arg1check(cbflag,nflag,nmax)
  if(nflag>nmax) 
    msg = sprintf('Flag %s needs one argument',cbflag);
    qoe(msg);error(msg);
  end
return

%%%------------------------------------------%%%%
function [havg, hstd, Nh] = hsa_convert(hsa,Nnnc)
  Nc = Nnnc+1;

  if(size(hsa,2) == 1 & length(size(hsa)) == 2)
    ud.nRows = 1;
    ud.nCols = 1;
    Nch2 = size(hsa,1);
  else
    [ud.nRows ud.nCols Nch2] = size(hsa);
  end
  Nh = Nch2/(2*Nc);

  hsa2 = permute(hsa, [3 2 1]); % t,c,r
  clear hsa;
  hsa3 = reshape(hsa2, [Nh 2 Nc ud.nCols ud.nRows ]); % h stat cond col row
  clear hsa2;

  hsa4 = permute(hsa3, [5 4 3  1 2 ]); % row col cond h stat
  clear hsa3;

  havg = (hsa4(:,:,:,:,1)); % col row cond havg 
  hstd = (hsa4(:,:,:,:,2)); % col row cond hstd

  %havg = squeeze(hsa4(:,:,:,:,1)); % col row cond havg 
  %hstd = squeeze(hsa4(:,:,:,:,2)); % col row cond hstd
  clear hsa4;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function set_text(hYak,ud)

  s = sprintf('%5d,%3s, %4.2f, %4.2f',...
              ud.CurPlaneNo,ud.OverlaySign,ud.PMin,ud.PMax);
  xr =  get(hYak,'XLim');
  yr =  get(hYak,'YLim');

  h = gettaghandle(hYak,'txStatus');
  if(~isempty(h))
    set(h,'string',s);
    set(h,'position',[xr(2) yr(1) 0]);
  else
    h = text(xr(2),yr(1),s);
    set(h,'Tag','txStatus');
    set(h,'VerticalAlignment','Bottom');
    set(h,'HorizontalAlignment','Right');
    set(h,'FontUnits','points');
    set(h,'FontSize',10);
    set(h,'EraseMode','background');
  end

    r = ud.CurPixel(1);
    c = ud.CurPixel(2);
    if(~isempty(ud.ActImg))
      [rf cf] = sind2find(r,c,ud);
      v =  ud.ActImg(rf,cf,ud.CurPlaneNo);
      % v =  sign(v)*10^(-abs(v));
      s = sprintf('r = %3d, c = %3d, (%2d,%2d), ovl = %6.4f',r,c,rf,cf,v);
    elseif(~ishandle2(ud.hHDR))
      v =  ud.UnderlayImg(r,c,ud.CurPlaneNo);
      [rf cf] = sind2find(r,c,ud);
      s = sprintf('r = %3d, c = %3d, (%2d,%2d), ovl = %6.4f',r,c,rf,cf,v);
    else
      v =  ud.UnderlayImg(r,c,ud.CurPlaneNo);
      s = sprintf('r = %3d, c = %3d, v = %g',r,c,v);
    end

  h = text(xr(1),yr(1),s);
  set(h,'Tag','txMouse');
  set(h,'VerticalAlignment','Bottom');
  set(h,'FontUnits','points');
  set(h,'FontSize',10);
  set(h,'EraseMode','background');

return

% This is for when h is null (ishandle
% crashes in this case).
function ish = ishandle2(h)
  if(isempty(h)) 
    ish = 0;
    return;
  end
  ish = zeros(size(h));
  for n = 1:length(h),
    if(ishandle(h(n))) ish(n) = 1; end
  end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [r2, c2] = rcmap_eqfov(szimg1,r1,c1,szimg2)
% Maps the row and column in one image to the closest
% in another image assume that the fields-of-view are
% the same in both.  szimg1 = [nrows1 ncols1], ...

tmp = round( (szimg2 .* [r1 c1] ./ szimg1));
r2 = tmp(1);
c2 = tmp(2);

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function h = gettaghandle(hparent,tag)

  hchildren = get(hparent,'Children');

  for n = 1:length(hchildren)
    tmptag = get(hchildren(n),'tag');
    if(strcmp(tmptag,tag)) 
      h = hchildren(n);
      return;
    end
  end

  h = [];

return
