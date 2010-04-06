function varargout = mkanalysis_gui(varargin)
% MKANALYSIS_GUI M-file for mkanalysis_gui.fig
%      MKANALYSIS_GUI, by itself, creates a new MKANALYSIS_GUI or raises the existing
%      singleton*.
%
%      H = MKANALYSIS_GUI returns the handle to a new MKANALYSIS_GUI or the handle to
%      the existing singleton*.
%
%      MKANALYSIS_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MKANALYSIS_GUI.M with the given input arguments.
%
%      MKANALYSIS_GUI('Property','Value',...) creates a new MKANALYSIS_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mkanalysis_gui_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mkanalysis_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mkanalysis_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @mkanalysis_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before mkanalysis_gui is made visible.
function mkanalysis_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mkanalysis_gui (see VARARGIN)
global MkAnalysisName;
global MkAnalysisClone;

% Choose default command line output for mkanalysis_gui
handles.output = hObject;
handles.version = '$Id: mkanalysis_gui.m,v 1.22 2010/04/06 23:08:34 greve Exp $';
handles.saveneeded = 1;
handles.flac = [];
handles.clone = '';
handles.hrfplot = [];

% These two should agree with values in mkanalysis-sess
handles.DefaultTER     = .050; % 50 ms
handles.DefaultTWindow = 40; % 40 sec
handles.DefaultACFBins = 10; 

%handles = parse_args(handles,varargin);
if(isempty(MkAnalysisName))
  errordlg('Global Variable MkAnalysisName not set.');  
  fprintf('Global Variable MkAnalysisName not set.\n');  
  return;
end

set(handles.ebAnalysisName,'string',MkAnalysisName);
AnaLoaded = '';
if(isempty(MkAnalysisClone))
  if(exist(MkAnalysisName,'dir'))
    fprintf('Loading %s\n',MkAnalysisName);
    handles.flac = fast_ldanaflac(MkAnalysisName);
    if(isempty(handles.flac))
      msg = sprintf('Could not load %s\n',MkAnalysisName);
      errordlg(msg);
      fprintf('ERROR: %s\n',msg);
      return;
    end
    handles.originalflac = handles.flac;
    handles = setstate(handles);
    AnaLoaded = MkAnalysisName;
  end
else
  if(exist(MkAnalysisClone,'dir'))
    fprintf('Cloning %s\n',MkAnalysisClone);
    handles.flac = fast_ldanaflac(MkAnalysisClone);
    if(isempty(handles.flac))
      msg = sprintf('Could not load %s\n',MkAnalysisClone);
      errordlg(msg);
      fprintf('ERROR: %s\n',msg);
      return;
    end
    handles.originalflac = handles.flac;
    handles.flac.name = MkAnalysisName;
    handles.flac.ana.analysis = MkAnalysisName;
    handles.clone = MkAnalysisClone;
    handles = setstate(handles);
  end
  AnaLoaded = MkAnalysisClone;
end

% Determine whether this is from an old GUI so we know whether 
% to change rescaling target.
if(~isempty(AnaLoaded))
  tmpfile = sprintf('%s/fsfast.flac',AnaLoaded);
  if(~exist(tmpfile,'file'))
    % fsfast.flac does not exist, could be from an old gui or from a
    % cmd-line call to mkanalysis-sess
    if(~isempty(handles.flac.con))
      if(isfield(handles.flac.ana.con(1).cspec,'CondState'))
	handles.flac.inorm = 100;
      end
    end
  end
end

if(isempty(handles.flac))
  fprintf('Using default parameters\n');
  handles.flac = flacinit;
  handles.flac.name = MkAnalysisName;
  handles.flac.ana.analysis = MkAnalysisName;
  handles.originalflac = [];
  handles = setstate(handles);
end
handles.flac.creator = 'GUI-01';

% Fill contrast list box
ncontrasts = length(handles.flac.ana.con);
clear tmpstr;
for nthcon = 1:ncontrasts
  tmpstr{nthcon} = handles.flac.ana.con(nthcon).cspec.name;
end
tmpstr{nthcon+1} = 'Add Contrast';
set(handles.lbContrast,'string',tmpstr);
handles.editing_contrast = 0;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes mkanalysis_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mkanalysis_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%=====================================================================
function ebFuncStem_Callback(hObject, eventdata, handles)
handles.flac.funcstem = get(hObject,'String');
handles = setstate(handles);
guidata(hObject, handles);
return;

%=======================================================================
function ebTR_Callback(hObject, eventdata, handles)
tmp = get(handles.ebTR,'string');
if(isempty(tmp))
  errordlg('You cannot set TR to be empty');  
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end
handles.flac.TR = sscanf(tmp,'%f');
handles.flac.ana.TER = handles.flac.TR;
%if(handles.flac.ana.TER == -1) 
%  %handles.flac.ana.TER = handles.flac.TR/round(handles.flac.TR/handles.DefaultTER);  
%end
handles.flac.ana.timewindow = handles.flac.TR * round(handles.DefaultTWindow/handles.flac.TR);
handles.flac.ana.prestim = handles.flac.TR * round(4/handles.flac.TR);
handles = setstate(handles);
guidata(hObject, handles);
return;

% ========================================================================
function cbERBlock_Callback(hObject, eventdata, handles)
if(get(handles.cbERBlock,'value') == 1)
  handles.flac.ana.designtype = 'event-related';
else
  handles.flac.ana.designtype = 'abblocked';  
end
handles = setstate(handles);
guidata(hObject, handles);
return;

% ========================================================================
function ebParFile_Callback(hObject, eventdata, handles)
handles.flac.parfile = get(handles.ebParFile,'string');
handles = setstate(handles);
guidata(hObject, handles);
return;

% =========================================================
function cbHRFGamma_Callback(hObject, eventdata, handles)
newval = get(handles.cbHRFGamma,'value');
if(handles.flac.ana.gammafit & newval == 0)
  % Tried to uncheck; recheck and return.
  set(handles.cbHRFGamma,'value',1);
  return;
end
if(0)
  cont = DeleteContrastsDialog(handles);
  if(~cont)
    handles = setstate(handles);
    guidata(hObject, handles);
  return;
  end
  handles.flac.ana.con = [];
end
handles.flac.ana.firfit = 0;
handles.flac.ana.gammafit = 1;
handles.flac.ana.spmhrffit = 0;
handles = setstate(handles);
guidata(hObject, handles);
return;

% =========================================================
function cbHRFFIR_Callback(hObject, eventdata, handles)
newval = get(handles.cbHRFFIR,'value');
if(handles.flac.ana.firfit & newval == 0)
  % Tried to uncheck; recheck and return.
  set(handles.cbHRFFIR,'value',1);
  return;
end
if(0)
  cont = DeleteContrastsDialog(handles);
  if(~cont)
    handles = setstate(handles);
    guidata(hObject, handles);
    return;
  end
  handles.flac.ana.con = [];
end

handles.flac.ana.firfit = 1;
handles.flac.ana.gammafit = 0;
handles.flac.ana.spmhrffit = 0;
handles = setstate(handles);
guidata(hObject, handles);
return;

% =========================================================
function cbHRFSPMHRF_Callback(hObject, eventdata, handles)
newval = get(handles.cbHRFSPMHRF,'value');
if(handles.flac.ana.spmhrffit & newval == 0)
  % Tried to uncheck; recheck and return.
  set(handles.cbHRFSPMHRF,'value',1);
  return;
end
if(0)
  cont = DeleteContrastsDialog(handles);
  if(~cont)
    handles = setstate(handles);
    guidata(hObject, handles);
    return;
  end
  handles.flac.ana.con = [];
end
handles.flac.ana.firfit = 0;
handles.flac.ana.gammafit = 0;
handles.flac.ana.spmhrffit = 1;
handles = setstate(handles);
guidata(hObject, handles);
return;

% ===========================================================
function ebFIRTotTimeWin_Callback(hObject, eventdata, handles)
if(checkTR(handles))
  tmp = sprintf('%f',handles.flac.ana.timewindow);
  set(handles.ebFIRTotTimeWin,'string',tmp);
  return;
end

win = sscanf(get(handles.ebFIRTotTimeWin,'string'),'%f');
if(win == handles.flac.ana.timewindow) return; end

ter = handles.flac.ana.TER;
if(rem(win,ter) ~= 0)
  errordlg('Your Time Window must be an integer multiple of the TER');
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end

if(0)
  cont = DeleteContrastsDialog(handles);
  if(~cont)
    handles = setstate(handles);
    guidata(hObject, handles);
    return;
  end
  handles.flac.ana.con = [];
end

handles.flac.ana.timewindow = win;
handles = setstate(handles);
guidata(hObject, handles);

% ===========================================================
function ebFIRPreStim_Callback(hObject, eventdata, handles)
if(checkTR(handles))
  tmp = sprintf('%f',handles.flac.ana.timewindow);
  set(handles.ebFIRTotTimeWin,'string',tmp);
  return;
end

prestim = sscanf(get(handles.ebFIRPreStim,'string'),'%f');
if(prestim == handles.flac.ana.prestim) return; end

ter = handles.flac.ana.TER;
if(rem(prestim,ter) ~= 0)
  errordlg('Your PreStim must be an integer multiple of the TER');
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end

% Note: this operation does not change the number of regressors

handles.flac.ana.prestim = prestim;
handles = setstate(handles);
guidata(hObject, handles);
return;

% ===========================================================
function ebFIRTER_Callback(hObject, eventdata, handles)
if(checkTR(handles))
  tmp = sprintf('%f',handles.flac.ana.timewindow);
  set(handles.ebFIRTotTimeWin,'string',tmp);
  return;
end

ter = sscanf(get(handles.ebFIRTER,'string'),'%f');
if(ter == handles.flac.ana.TER) return; end

cont = DeleteContrastsDialog(handles);
if(~cont)
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end

win = sscanf(get(handles.ebFIRTotTimeWin,'string'),'%f');
prestim = sscanf(get(handles.ebFIRPreStim,'string'),'%f');
if(rem(win,ter) ~= 0)
  % This is not quite the right thing to do, should offer to
  % Change the window to make it consitent
  errordlg('Your Time Window must be an integer multiple of the TER');
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end

handles.flac.ana.con = [];
handles.flac.ana.TER = ter;
handles = setstate(handles);
guidata(hObject, handles);

return;

%=====================================================================
function ebFSD_Callback(hObject, eventdata, handles)
tmp = get(handles.ebFSD,'string');
if(isempty(tmp))
  errordlg('You cannot set FSD to be empty');  
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end
if(~isempty(find(tmp == ' ')))
  errordlg('FSD cannot have blanks');  
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end
handles.flac.fsd = tmp;
handles = setstate(handles);
guidata(hObject, handles);
return;

%=====================================================================
function ebRLF_Callback(hObject, eventdata, handles)
tmp = get(handles.ebRLF,'string');
if(~isempty(find(tmp == ' ')))
  errordlg('Run List File cannot have blanks');  
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end
handles.flac.runlistfile = tmp;
handles = setstate(handles);
guidata(hObject, handles);
return;

%=====================================================================
function ebTPExcludeFile_Callback(hObject, eventdata, handles)
tmp = get(handles.ebTPExcludeFile,'string');
if(~isempty(find(tmp == ' ')))
  errordlg('TP Exclude File cannot have blanks');  
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end
handles.flac.tpexcfile = tmp;
handles = setstate(handles);
guidata(hObject, handles);
return;


% ==============================================================
function ebGammaDelta_Callback(hObject, eventdata, handles)
tmp = get(handles.ebGammaDelta,'string');
val = sscanf(tmp,'%f');
handles.flac.ana.gamdelay = val;
handles = setstate(handles);
guidata(hObject, handles);
return;

% ==============================================================
function ebGammaTau_Callback(hObject, eventdata, handles)
tmp = get(handles.ebGammaTau,'string');
val = sscanf(tmp,'%f');
handles.flac.ana.gamtau = val;
handles = setstate(handles);
guidata(hObject, handles);
return;

% ==============================================================
function ebGammaAlpha_Callback(hObject, eventdata, handles)
tmp = get(handles.ebGammaAlpha,'string');
val = sscanf(tmp,'%f');
handles.flac.ana.gamexp = val;
handles = setstate(handles);
guidata(hObject, handles);
return;

% ==============================================================
function ebNCycles_Callback(hObject, eventdata, handles)
tmp = get(handles.ebNCycles,'string');
val = sscanf(tmp,'%d');
handles.flac.ana.ncycles = val;
handles = setstate(handles);
guidata(hObject, handles);
return;

% ==============================================================
function ebAnalysisName_Callback(hObject, eventdata, handles)
tmp = get(handles.ebAnalysisName,'string');
if(isempty(tmp))
  errordlg('You cannot set Analysis Name to be empty');  
  set(handles.ebAnalysisName,'string',handles.flac.name);
  return;
end
if(~isempty(find(tmp == ' ')))
  errordlg('Analysis Name cannot have blanks');  
  set(handles.ebAnalysisName,'string',handles.flac.name);
  return;
end
d = dir(tmp);
if(~isempty(d))
  qstring = sprintf('Analysis %s already exists. ',tmp);
  qstring = [qstring 'It will be overwritten if you continue.'];
  qstring = [qstring 'Do you want to continue with this operation?'];
  button = questdlg(qstring,'WARNING','yes','no','no');
  if(strcmp(button,'no'))
    set(hObject,'Value',handles.flac.ana.nconditions);
    set(handles.ebAnalysisName,'string',handles.flac.name);
    return;
  end
end
handles.flac.name = tmp;
handles.flac.ana.analysis = tmp;
handles = setstate(handles);
guidata(hObject, handles);
return;

% ==============================================================
function cbPeriodicDesign_Callback(hObject, eventdata, handles)
if(get(handles.cbPeriodicDesign,'value') == 1)
  handles.flac.ana.designtype = 'abblocked';  
else
  handles.flac.ana.designtype = 'event-related';
end
handles = setstate(handles);
guidata(hObject, handles);
return;

% ==============================================================
function pbSave_Callback(hObject, eventdata, handles)
ok = OKToSave(handles);
if(~ok) return; end
if(handles.flac.ana.gammafit | handles.flac.ana.spmhrffit)
  % These match the default values for mkanalysis-sess cmd line
  % Added on 4/9/09. Modified on 3/23/10
  handles.flac.ana.prestim = 0;
  TER = handles.flac.TR/round(handles.flac.TR/handles.DefaultTER);  
  handles.flac.ana.TER = TER;
  handles.flac.ana.timewindow = TER*floor(handles.DefaultTWindow/TER);
end
fprintf('Saving %s ... ',handles.flac.name);
fast_svana(handles.flac.name,handles.flac);
fprintf(' done\n');
handles.originalflac = handles.flac;
handles.clone = '';
handles = setstate(handles);
guidata(hObject, handles);
return;

% ==============================================================
function pbQuit_Callback(hObject, eventdata, handles)
if(handles.saveneeded)
  qs = 'You have edited this analysis without saving .';
  qs = sprintf('%sContinuing this operation will lose those edits.',qs);
  qs = sprintf('%sDo you want to continue?.',qs);
  button = questdlg(qs,'WARNING','yes','no','no');
  if(strcmp(button,'no')) return; end
end
delete(handles.figure1);
return;

% =======================================================
function pbPlot_Callback(hObject, eventdata, handles)
fprintf('here\n');
handles.hrfplot  
if(isempty(handles.hrfplot) | ~ishandle(handles.hrfplot))  
  fprintf('here2\n');
  handles.hrfplot = figure; 
  figure(handles.figure1);
  handles = setstate(handles);
  guidata(hObject, handles);
end
return;


% =======================================================
function cbMCExtReg_Callback(hObject, eventdata, handles)
val = get(handles.cbMCExtReg,'value');
if(val) handles.flac.ana.extreg = 'mcextreg';
else    handles.flac.ana.extreg = '';
end
handles = setstate(handles);
guidata(hObject, handles);
return;


% =======================================================
function cbWhiten_Callback(hObject, eventdata, handles)
val = get(handles.cbWhiten,'value');
if(val) handles.flac.acfbins = handles.DefaultACFBins;
else    handles.flac.acfbins = 0;
end
handles = setstate(handles);
guidata(hObject, handles);
return;

% =======================================================
function cbINorm_Callback(hObject, eventdata, handles)
val = get(handles.cbINorm,'value');
if(val) handles.flac.inorm = 100;
else    handles.flac.inorm = 0;
end
handles = setstate(handles);
guidata(hObject, handles);
return;

% --- Executes on button press in pbRestore.
function pbRestore_Callback(hObject, eventdata, handles)
if(~isempty(handles.clone)) name = handles.flac.name; end
handles.flac = handles.originalflac;
if(~isempty(handles.clone)) handles.flac.name = name; end
handles = setstate(handles);
guidata(hObject, handles);
return;

%-----------------------------------------------------
function handles = parse_args(handles,varargin)

% Convert to a strvcat'ed string array
inputargs = char(varargin{1});
if(isempty(inputargs))
  handles.flac = fast_ldflac;
  return;
end
ninputargs = size(inputargs,1);

narg = 1;
while(narg <= nargin)

  flag = deblank(inputargs(narg,:));
  narg = narg + 1;
  %fprintf(1,'Argument: %s\n',flag);
  if(~isstr(flag))
    flag
    fprintf(1,'ERROR: All Arguments must be a string\n');
    error;
  end

  switch(flag)

   case '-a',
    arg1check(flag,narg,ninputargs);
    analysisdir = deblank(inputargs(narg,:));
    handles.flac = fast_ldanaflac(analysisdir);
    set(handles.ebAnalysisName,'string',analysisdir);
    narg = narg + 1;

   case '-clone',
    arg1check(flag,narg,ninputargs);
    analysisdir = deblank(inputargs(narg,:));
    handles.flac = fast_ldanaflac(analysisdir);
    narg = narg + 1;

   otherwise
    fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
    s = [];
    return;

  end % --- switch(flag) ----- %

end % while(narg <= ninputargs)

return;
%--------------------------------------------------%
%% Check that there is at least one more argument %%
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;

% SETSTATE --------------------------------------------------%
function handles = setstate(handles)

flac = handles.flac;
ana = flac.ana;
set(handles.ebAnalysisName,'string',ana.analysis);
set(handles.ebFuncStem,'string',flac.funcstem);
if(~isempty(flac.funcstem))  
  set(handles.ebFuncStem,'backgroundcolor','white')
else
  set(handles.ebFuncStem,'backgroundcolor','red')
end
if(~isempty(flac.TR))  
  set(handles.ebTR,'backgroundcolor','white')
else
  set(handles.ebTR,'backgroundcolor','red')
end

set(handles.slNConditions,'value',flac.ana.nconditions);
tmp = sprintf('NConditions %d',flac.ana.nconditions);
set(handles.txNConditions,'string',tmp);

set(handles.ebFSD,'string',flac.fsd);
set(handles.ebRLF,'string',flac.runlistfile);
set(handles.ebTR,'string',flac.TR);
set(handles.ebRED,'string',flac.RefEventDur);
set(handles.ebTPExcludeFile,'string',flac.tpexcfile);
if(flac.inorm ~= 0) set(handles.cbINorm,'value',1);
else                set(handles.cbINorm,'value',0);
end
if(strcmp(ana.extreg,'mcextreg')) set(handles.cbMCExtReg,'value',1);
else                              set(handles.cbMCExtReg,'value',0);
end

if(flac.acfbins ~= 0) set(handles.cbWhiten,'value',1);
else                  set(handles.cbWhiten,'value',0);
end

tmp = sprintf('Polynomial Order %d',ana.PolyOrder);
set(handles.txPolyFit,'string',tmp);
set(handles.slPolyFit,'value',ana.PolyOrder);

set(handles.ebFIRTotTimeWin,'string',flac.ana.timewindow);
set(handles.ebFIRPreStim,'string',flac.ana.prestim);
if(flac.TR > 0 & flac.ana.TER < 0) flac.ana.TER = flac.TR; end
set(handles.ebFIRTER,'string',flac.ana.TER);

set(handles.cbAutoStimDur,'value',flac.autostimdur);
set(handles.txParFile,'enable','on');
set(handles.ebParFile,'enable','on');

if(strcmp(ana.designtype,'event-related') | ...
   strcmp(ana.designtype,'blocked'))
  EnableERBlock = 1;
  EnableABBlock = 0;
  % Enable ERBlock ----------------------
  set(handles.cbERBlock,'value',1);
  set(handles.txNConditions,'enable','on');
  set(handles.slNConditions,'enable','on');
  set(handles.cbHRFFIR,'enable','on');
  set(handles.cbHRFGamma,'enable','on');
  set(handles.cbHRFSPMHRF,'enable','on');
  set(handles.lbContrast,'enable','on');  
  set(handles.cbAutoStimDur,'enable','on');
  % Disable ABBlocked ------------------
  set(handles.cbPeriodicDesign,'value',0);
  set(handles.ebNCycles,'enable','off');
  set(handles.txNCycles,'enable','off');
else
  EnableERBlock = 0;
  EnableABBlock = 1;
  % Disable ERBlock ----------------------
  set(handles.cbERBlock,'value',0);
  set(handles.txNConditions,'enable','off');
  set(handles.slNConditions,'enable','off');
  set(handles.cbHRFFIR,'enable','off');
  set(handles.cbHRFGamma,'enable','off');
  set(handles.cbHRFSPMHRF,'enable','off');
  set(handles.lbContrast,'enable','off');  
  set(handles.cbAutoStimDur,'enable','off');
  % Enable ABBlocked ------------------
  set(handles.cbPeriodicDesign,'value',1);
  set(handles.ebNCycles,'enable','on');
  set(handles.txNCycles,'enable','on');
end
set(handles.ebNCycles,'string',handles.flac.ana.ncycles);

set(handles.cbHRFFIR,'value',ana.firfit);
set(handles.cbHRFGamma,'value',ana.gammafit);
set(handles.cbHRFSPMHRF,'value',ana.spmhrffit);

if(EnableERBlock && ana.firfit)   FIREnable   = 'on';
else                              FIREnable   = 'off';
end
if(EnableERBlock && ana.gammafit)  GammaEnable   = 'on';
else                               GammaEnable   = 'off';
end
if(EnableERBlock && ana.spmhrffit) SPMHRFEnable   = 'on';
else                               SPMHRFEnable   = 'off';
end

if(EnableERBlock && isempty(handles.flac.parfile))
  set(handles.ebParFile,'backgroundcolor','red')
else
  set(handles.ebParFile,'backgroundcolor','white')
end
set(handles.ebParFile,'string',flac.parfile);

set(handles.txFIRTotTimeWin,'enable',FIREnable);
set(handles.ebFIRTotTimeWin,'enable',FIREnable);
set(handles.txFIRPreStim,'enable',FIREnable);
set(handles.ebFIRPreStim,'enable',FIREnable);
set(handles.txFIRTER,'enable',FIREnable);
set(handles.ebFIRTER,'enable',FIREnable);

set(handles.txGammaDelta,'enable',GammaEnable);
set(handles.ebGammaDelta,'enable',GammaEnable);
set(handles.txGammaTau,'enable',GammaEnable);
set(handles.ebGammaTau,'enable',GammaEnable);
set(handles.txGammaAlpha,'enable',GammaEnable);
set(handles.ebGammaAlpha,'enable',GammaEnable);

set(handles.txSPMHRFNDeriv,'enable',SPMHRFEnable);
set(handles.txSPMHRFNDeriv,'string',...
  sprintf('NDerivatives = %d',handles.flac.ana.nspmhrfderiv));
set(handles.slSPMHRFNDeriv,'enable',SPMHRFEnable);

if(~isfield(handles.flac.ana,'con'))
  handles.flac.ana.con = [];
end


ncontrasts = length(handles.flac.ana.con);
clear tmpstr;
for nthcon = 1:ncontrasts
  tmpstr{nthcon} = handles.flac.ana.con(nthcon).cspec.name;
end
tmpstr{nthcon+1} = 'Add Contrast';
set(handles.lbContrast,'string',tmpstr);

if(~isempty(handles.clone))
  tmp = sprintf('Cloned from: %s\n',handles.clone);
  set(handles.txClone,'string',tmp);
end

% nregressors per factor
if(ana.gammafit)  nregressors = 1; end
if(ana.spmhrffit) nregressors = ana.nspmhrfderiv + 1; end
if(ana.firfit)    nregressors = round(ana.timewindow/ana.TER); end

% Check whether the number of regessors per factor has changed.
if(handles.flac.ana.nregressors ~= nregressors)
  for nthcon = 1:ncontrasts
    cspec = handles.flac.ana.con(nthcon).cspec;

    if(handles.flac.ana.nregressors < nregressors)
      % Number of regessors per factor has increased
      % If manually set reg weights, fill new ones in with 0
      % Otherwise fill with 1s
      if(cspec.setwdelay) cspec.WDelay(end:nregressors) = 0;
      else cspec.WDelay(1:nregressors) = 1;
      end
    else
      % Number of regessors per factor has decreased
      % Just throw away to now unused ones
      cspec.WDelay = cspec.WDelay(1:nregressors);
    end

    % This part is pretty hideous
    if(ana.firfit)
      cspec.TPreStim   = ana.prestim;
      cspec.TimeWindow = ana.timewindow;
    else
      cspec.TPreStim = 0;
      cspec.TimeWindow = nregressors * ana.TER;
    end
    
    cspec.ContrastMtx_0 = fast_contrastmtx(cspec);
    handles.flac.ana.con(nthcon).cspec = cspec;
  end
  handles.flac.ana.nregressors = nregressors;
end


if(ishandle(handles.hrfplot))
  t = [0:.1:32];
  hgamma = fmri_hemodyn(t, ana.gamdelay, ana.gamtau, ana.gamexp);
  hgamma = hgamma/max(hgamma);
  hspmhrf = fast_spmhrf(t);
  hspmhrf = hspmhrf/max(hspmhrf);
  figure(handles.hrfplot);
  plot(t,hgamma,t,hspmhrf);
  legend('Gamma','SPM HRF');
  xlabel('Time (sec)');
  title('Hemodynamic Response Function');
  figure(handles.figure1);  
end

handles.saveneeded = SaveNeeded(handles);
%fprintf('Returned from SaveNeeded = %d\n',handles.saveneeded);
if(handles.saveneeded)
  set(handles.pbSave,'enable','on');
  set(handles.pbRestore,'enable','on');
else
  %set(handles.pbSave,'enable','off');
  set(handles.pbRestore,'enable','off');
end

subject = handles.flac.subject;
hemi = handles.flac.hemi;
UseTal = handles.flac.UseTalairach;
if(UseTal == 0 & isempty(subject)) v = 1; end
if(~isempty(subject))
  if(strcmp(subject,'fsaverage') & strcmp(hemi,'lh')) v = 2; end
  if(strcmp(subject,'fsaverage') & strcmp(hemi,'rh')) v = 3; end
  if(strcmp(subject,'self') & strcmp(hemi,'lh')) v = 5; end
  if(strcmp(subject,'self') & strcmp(hemi,'rh')) v = 6; end
end
if(UseTal == 1) v = 4; end
set(handles.menuAnatomy,'value',v);

return;

%------------------------------------------------------
function err = checkTR(handles)
err = 1;
if(isempty(handles.flac.TR))
  errordlg('You must set the TR before you can set this value.');  
  return;
end
err = 0;
return;

%------------------------------------------------------
function flac = flacinit

flac = fast_ldflac;
flac.fsd = 'bold';
flac.runlistfile = '';
flac.inorm = 100; % Changed from 100 to 1000 4/9/09, back to 100 4/6/10
flac.TR = '';
flac.parfile = '';
flac.acfbins = 10; % set to 0 to turn off whitening
flac.stimulusdelay = 0;
flac.autostimdur = 0;
flac.ana.designtype = 'event-related';
flac.ana.nconditions = 1;
flac.ana.PolyOrder = 2;
flac.ana.extreg = '';
flac.ana.nextreg = 3;
flac.ana.firfit = 0;
flac.ana.timewindow = 26;
flac.ana.prestim = 4;
flac.ana.TER = -1;
flac.ana.gammafit = 1;
flac.ana.gamdelay = 2.25;
flac.ana.gamtau = 1.25;
flac.ana.gamexp = 2;
flac.ana.spmhrffit = 0;
flac.ana.nspmhrfderiv = 0;
flac.ana.ncycles = 4;
flac.ana.con = [];
flac.ana.nregressors = 1;

flac.ana.ConditionNames = '';
for n = 1:flac.ana.nconditions
  tmp = sprintf('Condition%02d',n);
  flac.ana.ConditionNames = strvcat(flac.ana.ConditionNames,tmp);
end


return;


% --- Executes on slider movement.
function slPolyFit_Callback(hObject, eventdata, handles)
% hObject    handle to slPolyFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.flac.ana.PolyOrder = get(hObject,'Value');
handles = setstate(handles);
guidata(hObject, handles);
return;

% ===================================================================
function lbContrast_Callback(hObject, eventdata, handles)
if(handles.editing_contrast) return; end
contents = get(hObject,'String');
connum = get(hObject,'Value');
ncontrasts = length(handles.flac.ana.con);
if(connum <= ncontrasts)
  hMkCon = mkcontrast_gui('init',handles.figure1,handles.flac,connum)
else
  hMkCon = mkcontrast_gui('init',handles.figure1,handles.flac);
end
handles.editing_contrast = 1;
guidata(hObject, handles);
uiwait(hMkCon);
ud = get(handles.figure1,'userdata');
if(~isempty(ud.cspec))
  if(~isfield(ud.cspec,'delete'))
    % Add contrast to the list
    handles.flac.ana.con(connum).cspec = ud.cspec;
  else
    % Delete contrast from the list
    handles.flac.ana.con(connum:end-1) = ...
	handles.flac.ana.con(connum+1:end);
    handles.flac.ana.con(end) = [];
  end
  ncontrasts = length(handles.flac.ana.con);
  clear tmpstr;
  for nthcon = 1:ncontrasts
    tmpstr{nthcon} = handles.flac.ana.con(nthcon).cspec.name;
  end
  tmpstr{nthcon+1} = 'Add Contrast';
  set(handles.lbContrast,'string',tmpstr);
end
handles.editing_contrast = 0;
setstate(handles);
guidata(hObject, handles);
return;

% ========================================================================
function slNConditions_Callback(hObject, eventdata, handles)
cont = DeleteContrastsDialog(handles);
if(~cont)
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end
handles.flac.ana.con = [];
Nc = round(get(hObject,'Value'));
if(Nc > handles.flac.ana.nconditions)
  tmp = sprintf('Condition%02d',Nc);
  handles.flac.ana.ConditionNames = strvcat(handles.flac.ana.ConditionNames,tmp);
end
handles.flac.ana.nconditions = Nc;
handles = setstate(handles);
guidata(hObject, handles);
return;

% ========================================================================
function slNConditions_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);
set(hObject,'Min',1);
set(hObject,'Max',20);
step = 1/(20-1);
set(hObject,'SliderStep',[step step]);
return;

% ===============================================================
function slSPMHRFNDeriv_Callback(hObject, eventdata, handles)
if(0)
  cont = DeleteContrastsDialog(handles);
  if(~cont)
    handles = setstate(handles);
    guidata(hObject, handles);
    return;
  end
  handles.flac.ana.con = [];
end
Nderiv = round(get(hObject,'Value'));
handles.flac.ana.nspmhrfderiv = round(get(hObject,'Value'));
handles = setstate(handles);
setstate(handles);
guidata(hObject, handles);
return;

% ===============================================================
function slSPMHRFNDeriv_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor',[.9 .9 .9]);
set(hObject,'BackgroundColor',[.9 .9 .9]);
set(hObject,'Min',0);
set(hObject,'Max',3);
step = 1/(3-0);
set(hObject,'SliderStep',[step step]);
set(hObject,'Value',0);
return;

% ===============================================================
function ok = OKToSave(handles,quiet)
ok = 0;
if(nargin == 1) quiet = 0; end
flac = handles.flac;
ana = flac.ana;
if(isempty(flac.funcstem))
  if(quiet) return; end
  errordlg('You must set the funcstem before you can save analysis.');  
  return;
end
if(isempty(handles.flac.TR))
  if(quiet) return; end
  errordlg('You must set the TR before you can save analysis.');  
  return;
end
if(isempty(flac.funcstem))
  if(quiet) return; end
  errordlg('You must set the funcstem before you can save analysis.');  
  return;
end
if(strcmp(ana.designtype,'event-related') | ...
   strcmp(ana.designtype,'blocked'))
  if(isempty(handles.flac.parfile))
    if(quiet) return; end
    errordlg('You must set the paradigm file before you can save analysis.');  
    return;
  end
end
ok = 1;

return;

% =============================================================
function needed = SaveNeeded(handles)
needed = 1; 
a = handles.originalflac;
b = handles.flac;

if(isempty(a)) return; end
if(~strcmp(a.name,b.name)) return; end
if(~strcmp(a.fsd,b.fsd)) return; end
if(~strcmp(a.runlistfile,b.runlistfile)) return; end
if(~strcmp(a.funcstem,b.funcstem)) return; end
if(a.inorm ~= b.inorm) return; end
if(a.TR ~= b.TR) return; end
if(~strcmp(a.ana.extreg,b.ana.extreg)) return; end
if(a.ana.PolyOrder ~= b.ana.PolyOrder) return; end
if(~strcmp(a.ana.designtype,b.ana.designtype)) return; end
if(a.acfbins ~= b.acfbins) return; end
if(a.stimulusdelay ~= b.stimulusdelay) return; end

if(~strcmp(a.ana.designtype,'event-related') & ...
   ~strcmp(a.ana.designtype,'blocked'))
    if(a.ana.ncycles ~= b.ana.ncycles) return; end  
    needed = 0;
    return;
end
  
% Only gets here if event-related or blocked

if(~strcmp(a.parfile,b.parfile)) return; end
if(a.autostimdur ~= b.autostimdur) return; end
if(a.ana.nconditions ~= b.ana.nconditions) return; end
if(a.ana.firfit ~= b.ana.firfit) return; end
if(a.ana.firfit)
  if(a.ana.timewindow ~= b.ana.timewindow) return; end
  if(a.ana.prestim ~= b.ana.prestim) return; end
  if(a.ana.TER ~= b.ana.TER) return; end
end
if(a.ana.gammafit ~= b.ana.gammafit) return; end
if(a.ana.gammafit)
  if(a.ana.gamdelay ~= b.ana.gamdelay) return; end  
  if(a.ana.gamtau ~= b.ana.gamtau) return; end  
  if(a.ana.gamexp ~= b.ana.gamexp) return; end  
end
if(a.ana.spmhrffit ~= b.ana.spmhrffit) return; end
if(a.ana.spmhrffit)
  if(a.ana.nspmhrfderiv ~= b.ana.nspmhrfderiv) return; end  
end
a_ncontrasts = length(a.ana.con);
b_ncontrasts = length(b.ana.con);
if(a_ncontrasts ~= b_ncontrasts) return; end

ncontrasts = length(a.ana.con);
for nth = 1:ncontrasts
  Ca = a.ana.con(nth).cspec.ContrastMtx_0;
  Cb = b.ana.con(nth).cspec.ContrastMtx_0;
  if(size(Ca,1) ~= size(Cb,1)) return; end
  if(size(Ca,2) ~= size(Cb,2)) return; end
  if(Ca ~= Cb) return; end
end
needed = 0;

return;

%=======================================================
function cont = DeleteContrastsDialog(handles);
% cont = 1 means to continue op, deleting any contrasts that exist
cont = 0; 

ncontrasts = length(handles.flac.ana.con);
if(ncontrasts == 0) cont = 1; return; end

qstring = 'This operation will cause the number of regressors ';
qstring = [qstring 'to change. This will invalidate your '];
qstring = [qstring 'contrasts, forcing them to be deleted. '];
qstring = [qstring 'Do you want to continue with this operation?'];
button = questdlg(qstring,'WARNING','yes','no','no');
if(strcmp(button,'no')) return; end
cont = 1;

return;





% ===========================================================
% ===========================================================
function ebFIRPreStim_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;
% --- Executes during object creation, after setting all properties.
function ebFuncStem_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;
% --- Executes during object creation, after setting all properties.
function ebFSD_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;
% --- Executes during object creation, after setting all properties.
function ebRLF_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;
% --- Executes during object creation, after setting all properties.
function ebTR_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;
% --- Executes during object creation, after setting all properties.
function ebTPExcludeFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebTPExcludeFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% ========================================================================
function ebParFile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function lbContrast_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;
% ===================================================================
function slPolyFit_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
return;
% ===================================================================
function ebAnalysisName_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;

% --- Executes during object creation, after setting all properties.
function ebNCycles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebNCycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function ebGammaAlpha_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebGammaAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function ebGammaTau_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebGammaTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function ebGammaDelta_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebGammaDelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes during object creation, after setting all properties.
function ebFIRTER_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;
% --- Executes during object creation, after setting all properties.
function ebFIRTotTimeWin_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;


















%%%%%% May be orphaned %%%%%%%%%%%%%%%%%%%
function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function FigWBD(hObject, eventdata, handles)
fprintf('WBD\n');
return;


% --- Executes on button press in cbAutoStimDur.
function cbAutoStimDur_Callback(hObject, eventdata, handles)
% hObject    handle to cbAutoStimDur (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of cbAutoStimDur
newval = get(handles.cbAutoStimDur,'value');
handles.flac.autostimdur = newval;
handles = setstate(handles);
guidata(hObject, handles);
return;


%========================================================
function ebRED_Callback(hObject, eventdata, handles)
tmp = get(handles.ebRED,'string');
if(isempty(tmp))
  errordlg('You cannot set RED to be empty');  
  handles = setstate(handles);
  guidata(hObject, handles);
  return;
end
handles.flac.RefEventDur = sscanf(tmp,'%f');
fprintf('RED %g\n',handles.flac.RefEventDur);
handles = setstate(handles);
guidata(hObject, handles);
return;

% --- Executes during object creation, after setting all properties.
function ebRED_CreateFcn(hObject, eventdata, handles)
set(hObject,'BackgroundColor','white');
return;


% --- Executes on selection change in menuAnatomy.
function menuAnatomy_Callback(hObject, eventdata, handles)
fprintf('menuAnatomy\n');
item = get(handles.menuAnatomy,'Value');
fprintf('menuAnatomy %d\n',item);
switch (item)
  case 1,
   handles.flac.subject = '';
   handles.flac.hemi = '';
   handles.flac.UseTalairach = 0;
  case 2,
   handles.flac.subject = 'fsaverage';
   handles.flac.hemi = 'lh';
   handles.flac.UseTalairach = 0;
  case 3,
   handles.flac.subject = 'fsaverage';
   handles.flac.hemi = 'rh';
   handles.flac.UseTalairach = 0;
  case 4,
   handles.flac.subject = '';
   handles.flac.hemi = '';
   handles.flac.UseTalairach = 1;
  case 5,
   handles.flac.subject = 'self';
   handles.flac.hemi = 'lh';
   handles.flac.UseTalairach = 0;
  case 6,
   handles.flac.subject = 'self';
   handles.flac.hemi = 'rh';
   handles.flac.UseTalairach = 0;
end  
handles = setstate(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function menuAnatomy_CreateFcn(hObject, eventdata, handles)
% hObject    handle to menuAnatomy (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


