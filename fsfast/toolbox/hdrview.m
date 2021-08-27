function hfig = hdrview(varargin)
%
% hdrview(string,options) 
% hfig = hdrview('init',hsa,Nnnc,<TR>,<tPreStim>)
% hdrview('CurPixel',hFig,[r c]);
%
%


%
% hdrview.m
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
  msg = 'USAGE: hdrview(string,options)';
  qoe(msg);error(msg);
end

cbflag = varargin{1};

if(~isstr(cbflag))
  msg = 'First argument must be a string';
  qoe(msg);error(msg);
end

switch(cbflag)
  case 'init',
    if(nargin < 3)
      msg = sprintf('%s requires at least two arguments',cbflag);
      qoe(msg);error(msg);
    end
    ud = init_hdrview(varargin);
    hfig = figure;
    ud.hThisHDRView = hfig;
    set(hfig,'UserData',ud);
    hdr_redraw(hfig);
    set(hfig,'WindowButtonDownFcn','hdrview(''wbd'')');
    set(hfig,'KeyPressFcn',        'hdrview(''kbd'')');
    set(hfig,'DeleteFcn',          'hdrview(''delete'')');

  case 'reinit',
      ud = hdr_reinit(varargin);
      set(ud.hThisHDRView,'UserData',ud);
      hdr_redraw(ud.hThisHDRView);

  case 'wbd',
    %fprintf(1,'Window bd\n');

  case 'kbd',
    %fprintf(1,'Keypress bd\n');
    hdr_keypress;

  case 'delete',
    fprintf(1,'Closing HDR View\n');

  case 'CurPixel',
    if(nargin ~= 3)
      msg = 'CurPixel requires 3 arguments';
      qoe(msg);error(msg);
    end
    hfig = varargin{1,2};
    if(~ishandle(hfig))
      msg = sprintf('Invalid figure handle (%d)',hfig);
      qoe(msg);error(msg);
    end
    ud = get(hfig,'UserData');
    CurPixel = varargin{1,3};
    if(CurPixel(1) > ud.nRows | CurPixel(2) > ud.nCols)
      msg = sprintf('CurPixel (%d,%d) exceeds dimensions (%d,%d)',...
            CurPixel,ud.nRows,ud.nCols);
      qoe(msg);error(msg);
    end
    if(CurPixel(1) < 0 | CurPixel(2) < 0)
      msg = sprintf('CurPixel (%d,%d) cannot be less than zero',CurPixel);
      qoe(msg);error(msg);
    end
    ud.CurPixel = CurPixel;
    set(hfig,'UserData',ud);
    hdr_redraw(hfig);

  otherwise,
    msg = sprintf('Option %s unknown',cbflag);
    qoe(msg);error(msg);

end % ---- switch(cbflag) ---- %

return;
%%%------------------------------------------%%%%
%%%------- End of hdrview -------------------%%%%
%%%------------------------------------------%%%%

function hdr_keypress
  c = get(gcf,'CurrentCharacter'); 
  ud =   get(gcf,'UserData');
  n = str2num(c);

  if(~isempty(n))
     if(n < 1 | n > ud.Nnnc)
       msg = sprintf('Condition %d is invalid (1,%d)',n,ud.Nnnc);
       yr =  get(gca,'YLim');
       xr =  get(gca,'XLim');
       xm =  (xr(2)-xr(1))/2 + xr(1);
       h = text(xm , yr(2), msg);
       set(h,'VerticalAlignment','Bottom');
       return;
       %qoe(msg);error(msg);
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
       fprintf(1,'Removing Display of Condition %d\n',n);
     end
     set(gcf,'UserData',ud);
     hdr_redraw(gcf);
     return;
  end

  switch(c)
    case 'b' % Set baseline to zero
      ud.BaselineZero = ~ud.BaselineZero;
      fprintf(1,'Toggling baseline to %d\n',ud.BaselineZero);
      set(gcf,'UserData',ud);
      hdr_redraw(gcf);

    case 'e' % Toggle Error Bars
      ud.ShowStd = ~ud.ShowStd;
      fprintf(1,'Toggling Display of Standard Deviation (%d)\n',ud.ShowStd);
      set(gcf,'UserData',ud);
      hdr_redraw(gcf);

    case 'h' % Toggle Hold on Axis Range
      ud.FixYRange = ~ud.FixYRange;
      fprintf(1,'Toggling Hold on YRange (%d)\n',ud.FixYRange);
      ud.YRange = get(gca,'YLim');
      set(gcf,'UserData',ud);
      hdr_redraw(gcf);

    case 's' % Save HDRs in a column file
      curdir = pwd;
      cd(ud.SvDir);
      [fname pname] = uiputfile(ud.SvFile,'Save Hemodynamic Responses');
      cd(curdir);
      if(fname ~= 0)
        ud.SvFile = fname;
        if(isempty(pname)) ud.SvDir  = '.';
        else               ud.SvDir  = pname;
        end

        set(gcf,'UserData',ud);
        svname = strcat(ud.SvDir,ud.SvFile);

        r = ud.CurPixel(1);
        c = ud.CurPixel(2);
        havg = squeeze(ud.havg(r,c,:,:));
        if(ud.BaselineZero)
          blc = repmat(havg(:,1), [1 size(havg,2)]);
          havg = havg - blc;
        end

        hstd = squeeze(ud.hstd(r,c,:,:));
        n = ud.Nnnc*2+1;
        tmp = zeros(n,ud.Nh);
        tmp(1,:) = ud.t;
        tmp(2:2:n,:) = havg;
        tmp(3:2:n,:) = hstd;
        fmt = repmat('%g ',[1 n]);
        fmt = strcat(fmt,'\n');
        fid = fopen(svname,'w');
        if(fid == -1)
          fprintf(2,'ERROR: cannot open %s\n',svname);
        else
          fprintf(fid,'# Hemodynamic: Row = %3d, Col = %3d\n',r,c); 
          if(ud.BaselineZero)
            fprintf(fid,'# Baseline Zero: ');
            fprintf(fid,'%g ',blc(:,1));
            fprintf(fid,'\n');
          end
          fprintf(fid,'# Time Avg1 Std1 Avg2 Std2 ...\n');
          fprintf(fid,fmt,tmp); %'
          fclose(fid);
        end

      end

    otherwise,
      msg = sprintf('Unkown keystroke %c',c);

  end %--- switch(c) --- %
return

%%%------------------------------------------%%%%
function hdr_redraw(hfig)
  if(~ishandle(hfig))
    msg = sprintf('Invalid figure handle (%d)',hfig);
    qoe(msg);error(msg);
  end
  ud = get(hfig,'UserData');
  r = ud.CurPixel(1);
  c = ud.CurPixel(2);
  figure(hfig);

  havg = squeeze(ud.havg(r,c,ud.DisplayCond,:));

  % subtract the first component, if desired %
  if(ud.BaselineZero)
    % blc = repmat(havg(1,:), [size(havg,1) 1]);
    blc = repmat(havg(:,1), [1 size(havg,2)]);
    havg = havg - blc;
  end

  nDC = length(ud.DisplayCond); % number of display conditions

  trange = [(min(ud.t)-ud.TR/2) (max(ud.t)+ud.TR/2)];
  if(ud.ShowStd) 
    hstd = squeeze(ud.hstd(r,c,ud.DisplayCond,:));
    t = repmat(ud.t, [nDC 1]);
    hc = errorbar(t',havg',hstd','o-');  %'
  else  
    if( size(havg,1) == size(havg,2))
      hc = plot(ud.t,havg','o-'); % This is a hack %
    else
      hc = plot(ud.t,havg,'o-');
    end
  end
  hold on;
  plot(ud.t,zeros(size(ud.t)),'k-.');
  hold off;

  for n=nDC:-1:1 % reverse order
    set(hc(n),'Color',ud.CondColor(ud.DisplayCond(n),:));
    if(ud.ShowStd) 
      m = nDC+n;
      set(hc(m),'Color',ud.CondColor(ud.DisplayCond(n),:));
    end
  end

  set(gca,'XLim', trange);
  if(ud.FixYRange)
    set(gca,'YLim', ud.YRange);
  end

  s = printstatus(ud);
  yr =  get(gca,'YLim');
  h = text(0, yr(2), s);
  set(h,'VerticalAlignment','Bottom');
  set(h,'FontSize',10);
  set(h,'EraseMode','background');

return
%%%------------------------------------------%%%%
function ud = hdr_reinit(varargin)
  hfig = varargin{1}{2};
  ud = get(hfig,'UserData');

  hsa  = varargin{1}{3};
  Nch = size(hsa,3);
  [ud.havg ud.hstd ud.Nh] = hsa_convert(hsa,ud.Nnnc);

  if(length(varargin{1})>4) 
    ud.TR = varargin{1}{5}; 
  end
  if(length(varargin{1})>5) 
    ud.tPreStim = varargin{1}{6}; 
  end
  % ud.DisplayCond = [1:ud.Nnnc];
  ud.t = ud.TR*[0:ud.Nh-1] - ud.tPreStim;
  [ud.nRows ud.nCols dummy] = size(hsa);

return

%%%------------------------------------------%%%%
function ud = init_hdrview(varargin)
  ud = new_hdr_data;
  hsa = varargin{1}{2};
  ud.Nnnc = varargin{1}{3};
  Nch = size(hsa,3);
  [ud.havg ud.hstd ud.Nh] = hsa_convert(hsa,ud.Nnnc);

  if(length(varargin{1})>3) 
    ud.TR = varargin{1}{4}; 
  end
  if(length(varargin{1})>4) 
    ud.tPreStim = varargin{1}{5}; 
  end
  ud.DisplayCond = [1:ud.Nnnc];
  ud.t = ud.TR*[0:ud.Nh-1] - ud.tPreStim;
  [ud.nRows ud.nCols dummy] = size(hsa);
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

  hsa2 = permute(hsa, [3 2 1]);
  hsa3 = reshape(hsa2, [Nh 2 Nc ud.nRows ud.nCols]);

  havg = squeeze(hsa3(:,1,:,:,:));
  havg = havg(:,[2:Nc],:,:) - repmat(havg(:,1,:,:),[1 Nnnc 1 1]);
  havg = permute(havg, [4 3  2 1]);

  hstd = squeeze(hsa3(:,2,:,:,:));
  hstd = hstd(:,[2:Nc],:,:);
  hstd = permute(hstd, [4 3 2 1]);
return
%%%------------------------------------------%%%%
function s = printstatus(ud)
  s = repmat('*',[1 ud.Nnnc]);
  for n = 1:length(ud.DisplayCond),
    s(ud.DisplayCond(n)) = num2str(ud.DisplayCond(n));
  end

  if(ud.FixYRange) s = sprintf('%s YHold',s);
  else             s = sprintf('%s YFree',s);
  end

  if(ud.BaselineZero) s = sprintf('%s BLZero',s);
  else                s = sprintf('%s BLFree',s);
  end

  s = sprintf('%s (%2d,%2d)',s,ud.CurPixel);

  % fprintf('Current Status: %s\n',s);

return

%%%------------------------------------------%%%%
function ud = new_hdr_data
   ud = struct('hThisHDRView',    [], ...
               'havg',    [], ...
               'hstd',    [], ...
               't',       [],...
               'CondColor', [],...
               'Nnnc',        1, ...
               'Nh',          0, ...
               'nRows',       0, ...
               'nCols',       0, ...
               'tPreStim',    0, ...
               'TR',         2, ...
               'CurPixel',   [1 1], ...
               'DisplayCond', [1], ...
               'FixYRange',      0,...
               'BaselineZero',   0,...
               'YRange',   [],...
               'SvFile',   '*',...
               'SvDir',    '.',...
               'ShowStd',        0);

    ud.CondColor = ...
    [    0         0    1.0000;...
         0    0.5000         0;...
    1.0000         0         0;...
         0    0.7500    0.7500;...
    0.7500         0    0.7500;...
    0.7500    0.7500         0;...
    0.2500    0.2500    0.2500];

return

