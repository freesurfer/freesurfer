function hdrviewlite(varargin)
%
% hdrview(string,options) 
% hdrview('init',hFig,t,dof);
% hdrview('plot',hFig,havg,hstd, CurrPixel);
% hdrview('tmarker',tmarker)
% hdrview('kbd')
% hdrview('wbd')
%
%


%
% hdrviewlite.m
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
  msg = 'USAGE: hdrviewlite(string,options)';
  qoe(msg);error(msg);
end

cbflag = varargin{1};

if(~isstr(cbflag))
  msg = 'First argument must be a string';
  qoe(msg);error(msg);
end

switch(cbflag)
  case 'init',
    if(nargin ~= 4)
      msg = 'USAGE: hdrview(init,hFig,t,dof)';
      qoe(msg);error(msg);
    end

    ud = new_hdr_data;
    ud.ThisHDRView = varargin{2};
    ud.t    = varargin{3};
    ud.hdof = varargin{4};
    ud.Nc   = length(ud.hdof);
    ud.DisplayCond  = [1:ud.Nc-1];
    ud.Nh   = length(ud.t);
    ud.havg = zeros(ud.Nc,ud.Nh);
    ud.hstddev = zeros(ud.Nc,ud.Nh);
    ud.hstderr = zeros(ud.Nc,ud.Nh);
    ind = find(ud.hdof==0);
    ud.hdof(ind) = 1;

    set(ud.ThisHDRView,'UserData',ud);
    set(ud.ThisHDRView,'WindowButtonDownFcn','hdrviewlite(''wbd'')');
    set(ud.ThisHDRView,'KeyPressFcn',        'hdrviewlite(''kbd'')');

    hdr_redraw(ud);

  case 'tmarker',
    if(nargin ~= 3)
      msg = 'USAGE: hdrview(tmarker,hFig,tmarker)';
      qoe(msg);error(msg);
    end
    ud = get(varargin{2},'UserData');
    ud.ThisHDRView = varargin{2};
    ud.tmarker = varargin{3};
    set(ud.ThisHDRView,'UserData',ud);
    hdr_redraw(ud);

  case 'wbd',
    xyz = get(gca,'CurrentPoint');
    x = xyz(1,1);
    y = xyz(1,2);
    %fprintf('t = %g, v = %g\n',x,y);
    ud = get(gcf,'UserData');
    ud.CurPoint = [x y];
    set(gcf,'UserData',ud);
    hdr_redraw(ud);

  case 'kbd',
    %fprintf(1,'Keypress bd\n');
    hdr_keypress;

  % hdrview('plot',hFig,havg,hstd, CurrPixel);
  case 'plot',
    if(nargin ~= 5 & nargin ~= 6)
      msg = 'USAGE: hdrview(plot,hFig,havg,hstd, CurrPixel, <hoffset>)';
      qoe(msg);error(msg);
    end
    ud = get(varargin{2},'UserData');
    ud.havg = varargin{3};
    n = length(find(ud.havg(1,:)==0));
    if(n==ud.Nh) ud.SrcType = 'selxavg';
    else         ud.SrcType = 'selavg';
    end
    ud.hstddev = varargin{4};
    if(0)
     fprintf('hdof = ');
     fprintf('%d ',ud.hdof);
     fprintf('\n');
     fprintf('std = ');
     fprintf('%g ',ud.hstddev);
     fprintf('\n');
    end
    ud.hstderr = ud.hstddev./repmat(reshape1d(sqrt(ud.hdof)), [1 ud.Nh]);
    ud.CurPixel   = varargin{5};
    if(nargin == 6)
      ud.hoffset   = varargin{6};
      if(ud.ShowPercent & ud.hoffset == 0)
        fprintf('WARNING: current offset = 0, percent not implemented\n',1,1);
      end
    else
      ud.hoffset   = [];
    end

    set(ud.ThisHDRView,'UserData',ud);
    hdr_redraw(ud);

  otherwise,
    msg = sprintf('Option %s unknown',cbflag);
    qoe(msg);error(msg);

end % ---- switch(cbflag) ---- %

return;
%%%------------------------------------------%%%%
%%%------- End of hdrview -------------------%%%%
%%%------------------------------------------%%%%

function hdr_keypress
  ud =   get(gcf,'UserData');
  CurAxis = gca;
  c = get(ud.ThisHDRView,'CurrentCharacter'); 
  n = str2num(c);
  % fprintf('c = %s, n = %d\n',c,n);

   % if(~isempty(n))
  switch(c)
    case {'0','1','2','3','4','5','6','7','8','9'}
     if(n < 0 | n > ud.Nc)
       msg = sprintf('Condition %d is invalid (0,%d)',n,ud.Nc);
       yr =  get(CurAxis,'YLim');
       xr =  get(CurAxis,'XLim');
       xm =  (xr(2)-xr(1))/2 + xr(1);
       h = text(xm , yr(2), msg);
       set(h,'VerticalAlignment','Bottom');
       set(h,'FontUnits','points');
       set(h,'FontSize',10);
       return;
     end
     ind = find(ud.DisplayCond==n);
     if(isempty(ind)) 
        ud.DisplayCond = [ud.DisplayCond n];
        ud.DisplayCond = sort(ud.DisplayCond);
        fprintf(1,'Activating Display of Condition %d\n',n);
     else
       ind = find(ud.DisplayCond ~= n);
       if(isempty(ind))
         msg = 'Cannot clear the only trace';
         fprintf(2,'%s\n',msg);
         %qoe(msg);error(msg);
         return;
       end
       ud.DisplayCond = ud.DisplayCond(ind);
       % fprintf(1,'Removing Display of Condition %d\n',n);
     end
     set(ud.ThisHDRView,'UserData',ud);
     hdr_redraw(ud);
     return;

    case 'h' % Help
      s = sprintf('Keypress Commands: \n'); 
      s = sprintf(' %sd - toggle subtraction of condition 0\n',s);
      s = sprintf(' %se - toggle display of error bars\n',s);
      s = sprintf(' %sh - help (this message)\n',s);
      s = sprintf(' %sn - toggle display of condition n\n',s);
      s = sprintf(' %sp - display percent signal change \n',s);
      s = sprintf(' %sq - close HDRView \n',s);
      s = sprintf(' %ss - toggle display of std bars\n',s);
      s = sprintf(' %sv - save graph as ascii\n',s);
      s = sprintf(' %sz - subtract prestim avg for all conditions\n',s);
      msgbox(s,'HDR Help','help');
      % fprintf(s);

    case 'q' % quit
      fprintf(1,'Quiting HDRView\n');
      close(ud.ThisHDRView);

    case 'i' % toggle show of ideal
      ud.ShowIdeal = ~ud.ShowIdeal;
      fprintf(1,'Toggling Ideal %d\n',ud.ShowIdeal);
      set(ud.ThisHDRView,'UserData',ud);
      hdr_redraw(ud);

    case 'p' % toggle use of percent signal change
      if(isempty(ud.hoffset))
        s = sprintf('Percent cannot be used because\nno offset has been given.');
        msgbox(s,'Error','error');
      else
        ud.ShowPercent = ~ud.ShowPercent;
        fprintf(1,'Toggling Percent %d\n',ud.ShowPercent);
        set(ud.ThisHDRView,'UserData',ud);
        hdr_redraw(ud);
      end

    case 'z' % Set baseline to zero
      ud.BaselineZero = ~ud.BaselineZero;
      fprintf(1,'Toggling baseline to %d\n',ud.BaselineZero);
      set(ud.ThisHDRView,'UserData',ud);
      if(ud.BaselineZero)
        ud.BLZRange = find(ud.t <= 0);
        fprintf('Baseline computed from the first %d elements\n',length(ud.BLZRange));
      end
      hdr_redraw(ud);

    case 'd' % Subtract Condition 0
      ud.SubtrCond0 = ~ud.SubtrCond0;
      fprintf(1,'Toggling SubtrCond0 to %d\n',ud.SubtrCond0);
      set(ud.ThisHDRView,'UserData',ud);
      hdr_redraw(ud);

    case 'e' % Toggle StdError Bars
      if(ud.ShowBars == 2) ud.ShowBars = 0;
      else                 ud.ShowBars = 2;
      end
      % fprintf(1,'Toggling Display of Standard Error Bars (%d)\n',ud.ShowBars);
      set(ud.ThisHDRView,'UserData',ud);
      hdr_redraw(ud);

    case 'o' 
      ud.FixYRange = 1;
      ud.YRange = get(CurAxis,'YLim');
      title  = 'Adjust YMin and YMax';
      prompt = {'YMax:','YMin'};
      lines  = 1;
      def = {sprintf('%7.4f',ud.YRange(2)),sprintf('%7.4f',ud.YRange(1))};
      answer   = inputdlg(prompt,title,lines,def);
      ud.YRange(2) = sscanf(answer{1},'%f');
      ud.YRange(1) = sscanf(answer{2},'%f');
      fprintf(1,'Setting YRange to %g %g\n',ud.YRange(1),ud.YRange(2));
      set(ud.ThisHDRView,'UserData',ud);
      hdr_redraw(ud);


    case 's' % Toggle StdDev  Bars
      if(ud.ShowBars == 1) ud.ShowBars = 0;
      else                 ud.ShowBars = 1;
      end
      % fprintf(1,'Toggling Display of Standard Dev Bars (%d)\n',ud.ShowBars);
      set(ud.ThisHDRView,'UserData',ud);
      hdr_redraw(ud);

    case 'y' % Toggle Hold on Axis Range
      ud.FixYRange = ~ud.FixYRange;
      fprintf(1,'Toggling Hold on YRange (%d)\n',ud.FixYRange);
      ud.YRange = get(CurAxis,'YLim');
      set(ud.ThisHDRView,'UserData',ud);
      hdr_redraw(ud);

    case 'v' % Save HDRs in a column file
      curdir = pwd;
      cd(ud.SvDir);
      [fname pname] = uiputfile(ud.SvFile,'Save Hemodynamic Responses');
      cd(curdir);
      if(fname ~= 0)
        ud.SvFile = fname;
        if(isempty(pname)) ud.SvDir  = '.';
        else               ud.SvDir  = pname;
        end

        set(ud.ThisHDRView,'UserData',ud);
        svname = strcat(ud.SvDir,ud.SvFile);

if(0)
        r = ud.CurPixel(1);
        c = ud.CurPixel(2);
        havg = ud.havg;
        if(ud.SubtrCond0)
          havg0 = repmat(havg(1,:), [size(havg,1) 1]);
          havg = havg - havg0;
        end
        if(ud.BaselineZero)
          ud.BLZRange = find(ud.t <= 0);
	  b = mean(havg(:,[ud.BLZRange]),2);
	  blc = repmat(b, [1 size(havg,2)]);
          havg = havg - blc;
        end
        m = ud.Nc*3+1;
        hstddev = ud.hstddev;
        hstderr = ud.hstderr;
        tmp = zeros(m,ud.Nh);
        tmp(1,:) = ud.t;
        tmp(2:3:m,:) = havg;
        tmp(3:3:m,:) = hstddev;
        tmp(4:3:m,:) = hstderr;
        n = size(tmp,1);
        fmt = repmat('%g ',[1 n]);
        fmt = strcat(fmt,'\n');
end
        fid = fopen(svname,'w');
        if(fid == -1)
          fprintf(2,'ERROR: cannot open %s\n',svname);
        else
          % fprintf(fid,'# Hemodynamic: Row = %3d, Col = %3d\n',r,c); 
          % if(ud.BaselineZero)
          %   fprintf(fid,'# Baseline Zero: ');
          %   fprintf(fid,'%g ',blc(:,1));
          %   fprintf(fid,'\n');
          % end
          % fprintf(fid,'# Time Avg1 Std1 Avg2 Std2 ...\n');
          % fprintf(fid,fmt,tmp); %
          for n = 1:length(ud.t)
            fprintf(fid,'%4.1f ',ud.t(n));
            for c = 1:size(ud.tmp1,1)
              %fprintf(fid,'%7.3f ',ud.tmp1(c,n));
              fprintf(fid,'%g ',ud.tmp1(c,n));
            end
            for c = 1:size(ud.tmp1,1)
              %fprintf(fid,'%7.3f ',ud.tmp2(c,n));
              fprintf(fid,'%g ',ud.tmp2(c,n));
            end
  	    fprintf(fid,'\n');
          end
          fclose(fid);
        end

      end

    otherwise,
      msg = sprintf('Unkown keystroke %c',c);

  end %--- switch(c) --- %
return

%%%------------------------------------------%%%%
function hdr_redraw(ud)
  r = ud.CurPixel(1);
  c = ud.CurPixel(2);
  Nh = size(ud.havg,2);

  if(ud.ShowPercent & ~isempty(ud.hoffset) & ud.hoffset ~= 0)
    rescale = 100/abs(ud.hoffset);
  else
    rescale = 1;
  end
  % fprintf('rescale = %g\n',rescale);

  nDC = length(ud.DisplayCond); % number of display conditions
  havg = rescale*squeeze(ud.havg(ud.DisplayCond+1,:));

  if(ud.SubtrCond0)
    havg0 = repmat(ud.havg(1,:), [nDC 1]);
    havg = havg - havg0;
  end

  % subtract the first component, if desired %
  if(ud.BaselineZero)
    ud.BLZRange = find(ud.t <= 0);
    b = mean(havg(:,ud.BLZRange),2);
    blc = repmat(b, [1 size(havg,2)]);
    havg = havg - blc;
  end

  if(Nh ~= 1) trange = [min(ud.t) max(ud.t)];
  else        trange = [-1 1];   ud.t = 0;
  end

  hstddev = rescale*squeeze(ud.hstddev(ud.DisplayCond+1,:));
  hstderr = rescale*squeeze(ud.hstderr(ud.DisplayCond+1,:));

  figure(ud.ThisHDRView);
  if(ud.ShowBars==1)
    t = repmat(ud.t, [nDC 1]);
    hc = errorbar(t',havg',hstddev','o-');  %'
  elseif(ud.ShowBars==2)
    t = repmat(ud.t, [nDC 1]);
    hc = errorbar(t',havg',hstderr','o-');  %'
  else  
    if( size(havg,1) == size(havg,2))
      hc = plot(ud.t,havg','o-'); % 'This is a hack %
    else
      hc = plot(ud.t,havg,'o-');
    end
  end

  hold on;
  plot(ud.t,zeros(size(ud.t)),'k-.');
  yr =  get(gca,'YLim');
  plot([ud.t(ud.tmarker) ud.t(ud.tmarker)], yr, 'k-.');

  if(ud.ShowIdeal)
    hideal = 50*fmri_hemodyn(ud.t,ud.Delta,ud.Tau);
    plot(ud.t,hideal,'r-.');
  end
  hold off;

  for n=nDC:-1:1 % reverse order
    ndc = rem(ud.DisplayCond(n),7)+1;
    %fprintf('n=%d, ndc = %d\n',n,ndc);
    set(hc(n),'Color',ud.CondColor(ndc,:));
    if(ud.ShowBars ~= 0) 
      m = nDC+n;
      set(hc(m),'Color',ud.CondColor(ndc,:));
    end
  end
  set(gca,'FontUnits','points');
  set(gca,'FontSize',10);
  set(gca,'XLim', trange);
  if(ud.FixYRange)
    set(gca,'YLim', ud.YRange);
  end

  s = printstatus(ud);
  yr =  get(gca,'YLim');
  h = text(0, yr(2), s);
  set(h,'FontUnits','points');
  set(h,'FontSize',10);
  set(h,'VerticalAlignment','Bottom');
  set(h,'EraseMode','background');
  xlabel('Post-stimulus Delay (sec)')

  ud.tmp1 = havg;
  ud.tmp2 = hstderr;
  set(ud.ThisHDRView,'UserData',ud);

return
%%%------------------------------------------%%%%
function s = printstatus(ud)
  s = repmat('*',[1 ud.Nc]);

  for n = 1:length(ud.DisplayCond),
    m = rem(ud.DisplayCond(n),10);
    s(ud.DisplayCond(n)+1) = num2str(m);
  end

  if(ud.FixYRange) s = sprintf('%s YHold',s);
  else             s = sprintf('%s YFree',s);
  end

  if(ud.BaselineZero) s = sprintf('%s BLZero',s);
  else                s = sprintf('%s BLFree',s);
  end

  s = sprintf('%s (%2d,%2d)',s,ud.CurPixel);

  if(ud.ShowPercent)
    s = sprintf('%s %%',s);
  else
    s = sprintf('%s ',s);
  end

  if(~isempty(ud.hoffset))
    s = sprintf('%s offset=%7.2f',s,ud.hoffset);
  end

  hstderrmn = mean(reshape1d(ud.hstderr(2:ud.Nc,:)));

  s = sprintf('%s stderr=%5.3f ',s,hstderrmn);
  s = sprintf('%s t=%5.3f, v=%5.3f',s,ud.CurPoint);

  %fprintf('Current Status: %s\n',s);

return

%%%------------------------------------------%%%%
function ud = new_hdr_data
   ud = struct('hThisHDRView',    [], ...
               'havg',    [], ...
               'tmp1',    [], ...
               'tmp2',    [], ...
               'hoffset', [], ...
               'hstderr', [], ...
               'hstddev',    [], ...
               'hdof',    [], ...
               't',       [],...
               'tmarker',  1,...
               'CurPoint',  [0 0],...
               'CondColor', [],...
               'Nc',        1, ...
               'Nh',          0, ...
               'CurPixel',   [1 1], ...
               'DisplayCond', [1], ...
               'ShowPercent',    0,...
               'FixYRange',      0,...
               'BaselineZero',   0,...
               'BLZRange',   1,...
               'SubtrCond0',   1,...
               'ShowIdeal',   0,...
               'Delta',   2.25,...
               'Tau',   1.25,...
               'YRange',   [],...
               'SvFile',   '*',...
               'SvDir',    '.',...
               'ShowBars',        0,...
               'SrcType',   'selxavg');

    ud.CondColor = ...
    [    0         0    1.0000;...
         0    0.5000         0;...
    1.0000         0         0;...
         0    0.7500    0.7500;...
    0.7500         0    0.7500;...
    0.7500    0.7500         0;...
    0.2500    0.2500    0.2500];

return

