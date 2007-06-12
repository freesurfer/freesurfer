function msg = mkcontrast_gui(varargin)
% mkcontrast_gui('init',flac);
% mkcontrast_gui(cbstring);

msg = [];
if(nargin == 0)
  msg = 'mkcontrast_gui(''init'',flac);';
  fprintf('%s',msg);
  return;
end

Init = 0;
narg = 1;
while(narg <= nargin)
  cbflag = varargin{narg};
  narg = narg + 1;

  switch(cbflag)
    case 'init',
      argNcheck(cbflag,narg,nargin,1);
      Init = 1;
      ud = new_user_data;
      ud.flac = varargin{narg};
    otherwise
      ud = get(gcf,'UserData'); 
      if( ~isfield(ud,'hMkConGUI')) return; end
      handle_cb(cbflag);
  end % --- switch(cbflag) ----- %

  narg = narg + 1;
end %--------- while(narg <= nargin) ---------- %

if(Init)

  % Compute nregressors per condition
  if(ud.flac.ana.gammafit) 
    ud.nregressors = 1; 
  end
  if(ud.flac.ana.spmhrffit) 
    ud.nregressors = ud.flac.ana.nspmhrfderiv + 1;
  end
  if(ud.flac.ana.firfit)
    ud.nregressors = round(ud.flac.ana.timewindow/ud.flac.ana.TER);
  end

  ud.hMkConGUI = figure;

  figpos = get(gcf,'position');
  dx = figpos(3);
  dy = figpos(4);
  
  h = uicontrol('style', 'text','position',  [1 dy-40 100 20]);
  set(h,'string','Contrast Name','tag','txConName');
  set(h,'tooltip','Contrast Name');
  ud.txConName = h;
  
  h = uicontrol('style', 'edit','position', [55 dy-40 100 20]);
  set(h,'string','conname','tag','ebConName');
  set(h,'tooltip','Contrast Name');
  ud.ebConName = h;
  align([ud.txConName ud.ebConName],'Fixed',5,'Bottom');
  
  h = uicontrol('style', 'checkbox','position',  [1 dy-80 140 20]);
  set(h,'string','Sum Conditions','tag','cbSumConditions','value',1);
  ud.cbSumConditions = h;
  
  h = uicontrol('style', 'checkbox','position',  [1 dy-80 160 20]);
  set(h,'string','Sum Delays/Regressors','tag','cbSumRegressors','value',0);
  if(0 & ud.nregressors < 2) set(h,'visible','off'); end
  ud.cbSumRegressors = h;
  
  h = uicontrol('style', 'checkbox','position',  [1 dy-80 200 20]);
  set(h,'string','Remove Prestim','tag','cbRmPreStim','value',0);
  if(0 & ~ud.flac.ana.firfit) set(h,'visible','off'); end
  ud.cbRmPreStim = h;
  
  h = uicontrol('style', 'checkbox','position',  [1 dy-100 300 20]);
  set(h,'string','Set Condition Weights Manually',...
	'tag','cbManConWeights','value',0);
  ud.cbManConWeights = h;
  
  align([ud.txConName ud.cbSumConditions],'Left',5,'None');
  align([ud.cbSumConditions ud.cbSumRegressors ud.cbRmPreStim],'Fixed',5,'Bottom');

  for nth = 1:ud.flac.ana.nconditions
    ndy = dy-130-(nth-1)*35;

    hpan = uipanel('units','pixels','position',[1 ndy 250 30]);
    ud.panCondition(nth) = hpan;

    h = uicontrol('parent',hpan,'style','text',...
		  'position',[1 1 100 20]);
    hstring = sprintf('Condition %d',nth);
    htag = sprintf('txCondition%02d',nth);
    set(h,'string',hstring,'tag',htag);
    ud.txCondition(nth) = h;
  
    bgh = uibuttongroup('parent',hpan,'units','pixels','position',[140 1 68 25]);
    align([ud.txCondition(nth) bgh],'Fixed',5,'None');
    ud.bgCondition(nth) = bgh;

    h = uicontrol(bgh,'style','radiobutton','position',[1 1 20 20]);
    ud.rbConditionAct(nth) = h;

    h = uicontrol(bgh,'style','radiobutton','position',[20 1 20 20]);
    ud.rbConditionCtl(nth) = h;
    
    h = uicontrol(bgh,'style','radiobutton','position',[40 1 20 20]);
    ud.rbConditionIgnore(nth) = h;
    
    h = uicontrol(bgh,'style','edit','position',[75 1 30 20]);
    ud.ebConditionW(nth) = h;
    
  end
  ndy = ndy - 80;
  
  h = uicontrol('style', 'checkbox','position',  [1 ndy 160 20]);
  set(h,'string','Manually Set Delays/Regressor Weight',...
	'tag','cbManRegWeights','value',1);
  if(0 & ud.nregressors < 2) set(h,'visible','off'); end
  ud.cbManRegWeights = h;
  ndy = ndy - 60;
  
  hpan = uipanel('units','pixels','position',[1 ndy 400 60]);
  ud.panRegWeights = hpan;
  for nth = 1:ud.nregressors
    x = 25*(nth-1);
    h = uicontrol('parent',hpan,'style','text','position',[x 30 20 20]);
    hstring = sprintf('%d',nth);
    htag = sprintf('txRegWeight%02d',nth);
    set(h,'string',hstring,'tag',htag);
    ud.txRegWeight(nth) = h;
    
    h = uicontrol('parent',hpan,'style','edit','position',[x 1 20 20]);
    htag = sprintf('ebRegWeight%02d',nth);
    set(h,'string','1','tag',htag);
    ud.ebRegWeight(nth) = h;
  end
  
  set(gcf,'UserData',ud);

  return;
end % Initialize %

return;
%%%------------------------------------------%%%%
%%%------------------------------------------%%%%
%%%------------------------------------------%%%%






%%%------------------------------------------%%%%
function argNcheck(cbflag,nflag,nmax,nargs)
  if(nflag + nargs - 1 > nmax) 
    msg = sprintf('Flag %s needs %d arguments',cbflag,nargs);
    qoe(msg);error(msg);
  end
return

%-----------------------------------------%
function ud = new_user_data
ud.hMkConGUI = [];
return
