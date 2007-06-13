function hfig = mkcontrast_gui(varargin)
% h = mkcontrast_gui('init',hparent,flac,<connum>);
% mkcontrast_gui(cbstring);
% Creates a cspec field in the hparent UserData struct
%  If Cancel is hit, this field is there but empty

msg = [];
ud = [];
if(nargin == 0)
  msg = 'mkcontrast_gui(''init'',flac);';
  fprintf('%s',msg);
  return;
end

% ------ Parse the input args ----------
Init = 0;
narg = 1;
while(narg <= nargin)
  cbflag = varargin{narg};
  narg = narg + 1;

  switch(cbflag)
    case 'init',
      argNcheck(cbflag,narg,nargin,2);
      Init = 1;
      ud = new_user_data;
      ud.hparent = varargin{narg};
      ud.flac = varargin{narg+1};
      narg = narg + 2;
      if(nargin == 4)
	ud.newcon = 0;
	ud.connumber = varargin{narg};
	narg = narg + 1;
      else
	ud.newcon = 1;
	ud.connumber = length(ud.flac.ana.con) + 1;
	ud.flac.ana.con(ud.connumber).cspec = cspecinit(ud);
      end
      
    otherwise
      ud = get(gcf,'UserData'); 
      if(~isfield(ud,'hMkConGUI')) return; end
      handle_cb(varargin);
  end % --- switch(cbflag) ----- %

end %--------- while(narg <= nargin) ---------- %

% --------------- Initialize -------------------------------
if(Init)
  % Compute nregressors per condition
  hfig = figure;
  ud.hMkConGUI = hfig;
  set(ud.hMkConGUI,'DeleteFcn',  'mkcontrast_gui(''delete'');');
  set(ud.hMkConGUI,'MenuBar','none');

  nregressors = ud.flac.ana.nregressors;

  figpos = get(gcf,'position');
  dx = figpos(3);
  dy = figpos(4);
  ndy = dy - 40;
  
  % Contrast Name ----------------------
  h = uicontrol('style', 'text','position',  [1 ndy 100 20]);
  set(h,'string','Contrast Name','tag','txConName');
  set(h,'tooltip','Contrast Name');
  ud.txConName = h;
  h = uicontrol('style', 'edit','position', [55 ndy 100 20]);
  set(h,'string','conname','tag','ebConName');
  set(h,'tooltip','Contrast Name');
  set(h,'callback','mkcontrast_gui(''ebConName'');');
  ud.ebConName = h;

  % Done --------------------------
  h = uicontrol('style', 'pushbutton','position', [200 ndy 100 20]);
  set(h,'string','Done','tag','pbDone');
  set(h,'callback','mkcontrast_gui(''pbDone'');');
  ud.pbDone = h;
  
  % Cancel --------------------------
  h = uicontrol('style', 'pushbutton','position', [300 ndy 100 20]);
  set(h,'string','Cancel','tag','pbDone');
  set(h,'callback','mkcontrast_gui(''pbCancel'');');
  ud.pbCancel = h;
  
  align([ud.txConName ud.ebConName ud.pbDone ud.pbCancel],'Fixed',5,'Bottom');
  ndy = ndy - 25;

  % Conditions -------------------------
  h = uicontrol('style', 'checkbox','position',  [1 ndy 140 20]);
  set(h,'string','Normalize Rows','tag','cbCNorm');
  set(h,'callback','mkcontrast_gui(''cbCNorm'');');
  ud.cbCNorm = h;
  ndy = ndy - 20;
  
  h = uicontrol('style', 'checkbox','position',  [1 ndy 140 20]);
  set(h,'string','Sum Conditions','tag','cbSumConditions');
  set(h,'callback','mkcontrast_gui(''cbSumConditions'');');
  ud.cbSumConditions = h;
  ndy = ndy - 25;
  
  h = uicontrol('style', 'checkbox','position',  [1 ndy 300 20]);
  set(h,'string','Set Condition Weights Manually',...
	'tag','cbManConWeights','value',0);
  set(h,'callback','mkcontrast_gui(''cbManConWeights'');');
  ud.cbManConWeights = h;
  ndy = ndy - 30;
  
  for nth = 1:ud.flac.ana.nconditions
    ndy = ndy-(nth-1)*30;

    hpan = uipanel('units','pixels','position',[1 ndy 300 30]);
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
    cback = sprintf('mkcontrast_gui(''rbConditionAct'',%d);',nth);
    set(h,'callback',cback);
    ud.rbConditionAct(nth) = h;

    h = uicontrol(bgh,'style','radiobutton','position',[20 1 20 20]);
    cback = sprintf('mkcontrast_gui(''rbConditionCtl'',%d);',nth);
    set(h,'callback',cback);
    ud.rbConditionCtl(nth) = h;
    
    h = uicontrol(bgh,'style','radiobutton','position',[40 1 20 20]);
    cback = sprintf('mkcontrast_gui(''rbConditionIgnore'',%d);',nth);
    set(h,'callback',cback);
    ud.rbConditionIgnore(nth) = h;
    
    h = uicontrol(bgh,'style','edit','position',[75 1 60 20]);
    cback = sprintf('mkcontrast_gui(''ebConditionW'',%d);',nth);
    set(h,'callback',cback,'horizontalalignment','right');
    ud.ebConditionW(nth) = h;
    
  end
  ndy = ndy - 80;
  
  % ---------- Handle Regressor Weights -----------------------
  h = uicontrol('style', 'checkbox','position',  [1 ndy 260 20]);
  set(h,'string','Set Delay/Regressor Weights Manually',...
	'tag','cbManRegWeights','value',1);
  set(h,'callback','mkcontrast_gui(''cbManRegWeights'');');
  if(0 & nregressors < 2) set(h,'visible','off'); end
  ud.cbManRegWeights = h;
  ndy = ndy - 20;
  
  h = uicontrol('style', 'checkbox','position',  [1 ndy 160 20]);
  set(h,'string','Sum Delays/Regressors','tag','cbSumRegressors','value',0);
  set(h,'callback','mkcontrast_gui(''cbSumRegressors'');');
  if(0 & nregressors < 2) set(h,'visible','off'); end
  ud.cbSumRegressors = h;
  ndy = ndy - 20;
  
  h = uicontrol('style', 'checkbox','position',  [1 ndy 200 20]);
  set(h,'string','Remove Prestim','tag','cbRmPreStim','value',0);
  set(h,'callback','mkcontrast_gui(''cbRmPreStim'');');
  if(~ud.flac.ana.firfit) 
    set(h,'visible','off'); 
    set(h,'enable','off'); 
  else
    ndy = ndy - 20;
  end
  ud.cbRmPreStim = h;
  
  ndy = ndy - 40;
  hpan = uipanel('units','pixels','position',[1 ndy 400 60]);
  ud.panRegWeights = hpan;
  wdx = 40;
  for nth = 1:nregressors
    x = wdx*(nth-1);
    h = uicontrol('parent',hpan,'style','text','position',[x 30 wdx 20]);
    hstring = sprintf('%d',nth);
    htag = sprintf('txRegWeight%02d',nth);
    set(h,'string',hstring,'tag',htag);
    ud.txRegWeight(nth) = h;
    
    h = uicontrol('parent',hpan,'style','edit','position',[x 1 wdx 20]);
    htag = sprintf('ebRegWeight%02d',nth);
    set(h,'string','1','tag',htag,'horizontalalignment','right');
    cback = sprintf('mkcontrast_gui(''ebRegWeight'',%d);',nth);
    set(h,'callback',cback,'horizontalalignment','right');
    ud.ebRegWeight(nth) = h;
  end

  ud = setstate(ud);
  set(gcf,'UserData',ud);

  return;
end % Initialize %

return;
%%%------------------------------------------%%%%
%%%------------------------------------------%%%%
%%%------------------------------------------%%%%

function ud = handle_cb(varargin)
cbflag = varargin{1}{1};
ud = get(gcf,'UserData'); 
flac = ud.flac;
ana = flac.ana;
cspec = flac.ana.con(ud.connumber).cspec;
switch (cbflag)
 case 'ebConName', 
  tmp = get(ud.ebConName,'string');
  ud.flac.ana.con(ud.connumber).cspec.name = tmp;
  ud = setstate(ud);
 case 'cbSumConditions', 
  ud.flac.ana.con(ud.connumber).cspec.sumconds = ...
      get(ud.cbSumConditions,'value');
  ud = setstate(ud);
 case 'cbSumRegressors', 
  ud.flac.ana.con(ud.connumber).cspec.sumdelays = ...
      get(ud.cbSumRegressors,'value');
  ud = setstate(ud);
 case 'cbRmPreStim', 
  ud.flac.ana.con(ud.connumber).cspec.RmPreStim = ...
      get(ud.cbRmPreStim,'value');
  ud = setstate(ud);
 case 'cbManConWeights', 
  ud.flac.ana.con(ud.connumber).cspec.setwcond = ...
      get(ud.cbManConWeights,'value');
  ud = setstate(ud);
 case 'cbManRegWeights', 
  ud.flac.ana.con(ud.connumber).cspec.setwdelay = ...
      get(ud.cbManRegWeights,'value');
  ud = setstate(ud);
 case 'cbCNorm', 
  ud.flac.ana.con(ud.connumber).cspec.CNorm = ...
      get(ud.cbCNorm,'value');
  ud = setstate(ud);
 case 'rbConditionAct',
  c = varargin{1}{2};
  ud.flac.ana.con(ud.connumber).cspec.CondState(c) = 1;
  ud.flac.ana.con(ud.connumber).cspec.WCond(c) = 1;
  ud = setstate(ud);
 case 'rbConditionCtl',
  c = varargin{1}{2};
  ud.flac.ana.con(ud.connumber).cspec.CondState(c) = 2;
  ud.flac.ana.con(ud.connumber).cspec.WCond(c) = -1;
  ud = setstate(ud);
 case 'rbConditionIgnore',
  c = varargin{1}{2};
  ud.flac.ana.con(ud.connumber).cspec.CondState(c) = 0;
  ud.flac.ana.con(ud.connumber).cspec.WCond(c) = 0;
  ud = setstate(ud);
 case 'ebConditionW',
  c = varargin{1}{2};
  val = sscanf(get(ud.ebConditionW(c),'string'),'%f');
  ud.flac.ana.con(ud.connumber).cspec.WCond(c) = val;
  ud = setstate(ud);
 case 'ebRegWeight',
  w = varargin{1}{2};
  val = sscanf(get(ud.ebRegWeight(w),'string'),'%f')
  ud.flac.ana.con(ud.connumber).cspec.WDelay(w) = val;
  ud = setstate(ud);
 case 'pbDone',
  fprintf('Done\n');
  delete(ud.hMkConGUI);
  return;
 case 'pbCancel',
  fprintf('Cancel\n');
  pud = get(ud.hparent,'userdata');
  pud.cspec = [];
  set(ud.hparent,'userdata',pud);
  delete(ud.hMkConGUI);
  return;
 case 'delete',
  fprintf('Deleting\n');
  return;
 
 otherwise,
  msg = sprintf('flag %s unrecognized',cbflag);
  fprintf('%s\n',msg);

end %--------- switch (cbflag) ------------%

set(ud.hMkConGUI,'UserData',ud);
return;
 
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

%%%------------------------------------------%%%%
function cspec = cspecinit(ud)

flac = ud.flac;
ana = flac.ana;
Nc = ana.nconditions;

cspec.name = sprintf('contrast%02d',ud.connumber);
%cspec.name = '';
cspec.CNorm = 1;
cspec.sumconds = 1;
cspec.setwcond = 0;
cspec.CondState = zeros(1,Nc);
cspec.CondState(1) = 1;

cspec.WCond = zeros(1,ana.nconditions);
cspec.WCond(1) = 1;

cspec.setwdelay = 0;
cspec.sumdelays = 0;
cspec.RmPreStim = 0;
cspec.WDelay = ones(1,ana.nregressors);

cspec.NCond    = Nc;
cspec.TER      = ana.TER;
% This part is pretty hideous
if(ana.firfit)
  cspec.TPreStim = ana.prestim;
  cspec.TimeWindow = ana.timewindow;
else
  cspec.TPreStim = 0;
  cspec.TimeWindow = ana.nregressors * ana.TER;
end

cspec.nircorr = 0;
cspec.rdelta = [0 0];
cspec.rtau   = [0 0];
cspec.ContrastMtx_0 = fast_contrastmtx(cspec);

return;

%%%------------------------------------------%%%%
function ud = setstate(ud)

ana = ud.flac.ana;
cspec = ana.con(ud.connumber).cspec;
set(ud.ebConName,'string',cspec.name);
if(isempty(cspec.name))
  set(ud.ebConName,'backgroundcolor','red');
else
  set(ud.ebConName,'backgroundcolor','white');
end

set(ud.cbManConWeights,'value',cspec.setwcond);
set(ud.cbManRegWeights,'value',cspec.setwdelay);
set(ud.cbSumConditions,'value',cspec.sumconds);
set(ud.cbSumRegressors,'value',cspec.sumdelays);
set(ud.cbRmPreStim,'value',cspec.RmPreStim);
set(ud.cbCNorm,'value',cspec.CNorm);

for c = 1:ana.nconditions
  set(ud.ebConditionW(c),'string',sprintf('%6.4f',cspec.WCond(c)));
  if(cspec.setwcond)
    set(ud.rbConditionAct(c),'enable','off');
    set(ud.rbConditionCtl(c),'enable','off');
    set(ud.rbConditionIgnore(c),'enable','off');
    set(ud.ebConditionW(c),'enable','on');
  else
    set(ud.rbConditionAct(c),'enable','on');
    set(ud.rbConditionCtl(c),'enable','on');
    set(ud.rbConditionIgnore(c),'enable','on');    
    set(ud.ebConditionW(c),'enable','off');
    if(cspec.CondState(c) == 0)
      set(ud.rbConditionIgnore(c),'value',1);
    elseif(cspec.CondState(c) == 1)
      set(ud.rbConditionAct(c),'value',1);
    else
      set(ud.rbConditionCtl(c),'value',1);
    end
  end
end

for w = 1:ana.nregressors
  tmp = sprintf('%4.2f',cspec.WDelay(w));
  set(ud.ebRegWeight(w),'string',tmp);
  if(cspec.setwdelay)
    set(ud.ebRegWeight(w),'enable','on');
  else
    set(ud.ebRegWeight(w),'enable','off');
  end
end
cspec.ContrastMtx_0 = fast_contrastmtx(cspec);
ud.flac.ana.con(ud.connumber).cspec = cspec;

pud = get(ud.hparent,'userdata');
pud.cspec = cspec;
set(ud.hparent,'userdata',pud);

ud.cspec = cspec; % For easier debuggin
cspec.ContrastMtx_0
  
return

