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
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mkanalysis_gui

% Last Modified by GUIDE v2.5 11-Jun-2007 18:52:24

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
handles.version = '$Id: mkanalysis_gui.m,v 1.4 2007/06/11 23:23:35 greve Exp $';
handles.saveneeded = 1;
handles.flac = [];
handles.clone = '';

%handles = parse_args(handles,varargin);
if(isempty(MkAnalysisName))
  errordlg('Global Variable MkAnalysisName not set.');  
  fprintf('Global Variable MkAnalysisName not set.\n');  
  return;
end

set(handles.ebAnalysisName,'string',MkAnalysisName);
if(exist(MkAnalysisName,'dir'))
  fprintf('Loading %s\n',MkAnalysisName);
  handles.flac = fast_ldanaflac(MkAnalysisName);
  if(isempty(handles.flac))
    msg = sprintf('Could not load %s\n',MkAnalysisName);
    errordlg(msg);
    fprintf('ERROR: %s\n',msg);
    return;
  end
  handles.saveneeded = 0;    
  setstate(handles);
end

if(~isempty(MkAnalysisClone))
  if(exist(MkAnalysisClone,'dir'))
    fprintf('Cloning %s\n',MkAnalysisClone);
    handles.flac = fast_ldanaflac(MkAnalysisClone);
    if(isempty(handles.flac))
      msg = sprintf('Could not load %s\n',MkAnalysisClone);
      errordlg(msg);
      fprintf('ERROR: %s\n',msg);
      return;
    end
    handles.flac.name = MkAnalysisName;
    handles.flac.ana.analysis = MkAnalysisName;
    handles.clone = MkAnalysisClone;
    setstate(handles);
  end
end

if(isempty(handles.flac))
  fprintf('Using default parameters\n');
  handles.flac = flacinit;
  handles.flac.name = MkAnalysisName;
  handles.flac.ana.analysis = MkAnalysisName;
  setstate(handles);
end
  
handles.originalflac = handles.flac;

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
setstate(handles);
guidata(hObject, handles);
return;

% --- Executes during object creation, after setting all properties.
function ebFuncStem_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebFuncStem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ebFSD_Callback(hObject, eventdata, handles)
% hObject    handle to ebFSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebFSD as text
%        str2double(get(hObject,'String')) returns contents of ebFSD as a double


% --- Executes during object creation, after setting all properties.
function ebFSD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebFSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ebRLF_Callback(hObject, eventdata, handles)
% hObject    handle to ebRLF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebRLF as text
%        str2double(get(hObject,'String')) returns contents of ebRLF as a double


% --- Executes during object creation, after setting all properties.
function ebRLF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebRLF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%=======================================================================
function ebTR_Callback(hObject, eventdata, handles)
tmp = get(ebTR,'string');
if(isempty(tmp))
  errordlg('You cannot set TR to be empty');  
  tmp = sprintf('%f',handles.flac.TR);
  set(ebTR,'string',tmp);  
  return;
end
handles.flac.TR = sscanf(tmp,'%f');
setstate(handles);
guidata(hObject, handles);
return;

% --- Executes during object creation, after setting all properties.
function ebTR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebTR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ebTPExcludeFile_Callback(hObject, eventdata, handles)
% hObject    handle to ebTPExcludeFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebTPExcludeFile as text
%        str2double(get(hObject,'String')) returns contents of ebTPExcludeFile as a double


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


% ========================================================================
function cbERBlock_Callback(hObject, eventdata, handles)
if(get(handles.cbERBlock,'value') == 1)
  handles.flac.ana.designtype = 'event-related';
else
  handles.flac.ana.designtype = 'abblocked';  
end
setstate(handles);
guidata(hObject, handles);

% ========================================================================
function ebNConditions_Callback(hObject, eventdata, handles)

% ========================================================================
function ebNConditions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ========================================================================
function ebParFile_Callback(hObject, eventdata, handles)
handles.flac.par4file = get(handles.ebParFile,'string');
setstate(handles);
guidata(hObject, handles);

return;

% ========================================================================
function ebParFile_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% ----------------------------------------------------------
function cbHRFGamma_Callback(hObject, eventdata, handles)
newval = get(handles.cbHRFGamma,'value');
if(handles.flac.ana.gammafit & newval == 0)
  % Tried to uncheck; recheck and return.
  set(handles.cbHRFGamma,'value',1);
  return;
end
handles.flac.ana.firfit = 0;
handles.flac.ana.gammafit = 1;
handles.flac.ana.spmhrffit = 0;
setstate(handles);
guidata(hObject, handles);
return;

% ----------------------------------------------------------
function cbHRFFIR_Callback(hObject, eventdata, handles)
newval = get(handles.cbHRFFIR,'value');
if(handles.flac.ana.firfit & newval == 0)
  % Tried to uncheck; recheck and return.
  set(handles.cbHRFFIR,'value',1);
  return;
end
handles.flac.ana.firfit = 1;
handles.flac.ana.gammafit = 0;
handles.flac.ana.spmhrffit = 0;
setstate(handles);
guidata(hObject, handles);
return;

% ----------------------------------------------------------
function cbHRFSPMHRF_Callback(hObject, eventdata, handles)
newval = get(handles.cbHRFSPMHRF,'value');
if(handles.flac.ana.spmhrffit & newval == 0)
  % Tried to uncheck; recheck and return.
  set(handles.cbHRFSPMHRF,'value',1);
  return;
end
handles.flac.ana.firfit = 0;
handles.flac.ana.gammafit = 0;
handles.flac.ana.spmhrffit = 1;
setstate(handles);
guidata(hObject, handles);
return;

% --------------------------------------------------------------
function ebFIRTotTimeWin_Callback(hObject, eventdata, handles)
if(checkTR(handles))
  tmp = sprintf('%d',handles.flac.ana.timewindow);
  set(handles.ebFIRTotTimeWin,'string',tmp);
  return;
end
win = sscanf(get(handles.ebFIRTotTimeWin,'string'),'%f');
ter = handles.flac.ana.TER;
if(rem(win,ter) ~= 0)
  errordlg('Your Time Window must be an integer multiple of the TER');
  tmp = sprintf('%d',handles.flac.ana.timewindow);
  set(handles.ebFIRTotTimeWin,'string',tmp);
  return;
end
handles.flac.ana.timewindow = win;
setstate(handles);
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function ebFIRTotTimeWin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebFIRTotTimeWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ebSPMHRFNDeriv_Callback(hObject, eventdata, handles)
% hObject    handle to ebSPMHRFNDeriv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebSPMHRFNDeriv as text
%        str2double(get(hObject,'String')) returns contents of ebSPMHRFNDeriv as a double


% --- Executes during object creation, after setting all properties.
function ebSPMHRFNDeriv_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebSPMHRFNDeriv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ebFIRPreStim_Callback(hObject, eventdata, handles)
% hObject    handle to ebFIRPreStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebFIRPreStim as text
%        str2double(get(hObject,'String')) returns contents of ebFIRPreStim as a double


% --- Executes during object creation, after setting all properties.
function ebFIRPreStim_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebFIRPreStim (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ebFIRTER_Callback(hObject, eventdata, handles)
% hObject    handle to ebFIRTER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebFIRTER as text
%        str2double(get(hObject,'String')) returns contents of ebFIRTER as a double


% --- Executes during object creation, after setting all properties.
function ebFIRTER_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebFIRTER (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ebGammaDelta_Callback(hObject, eventdata, handles)
% hObject    handle to ebGammaDelta (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebGammaDelta as text
%        str2double(get(hObject,'String')) returns contents of ebGammaDelta as a double


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



function ebGammaTau_Callback(hObject, eventdata, handles)
% hObject    handle to ebGammaTau (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebGammaTau as text
%        str2double(get(hObject,'String')) returns contents of ebGammaTau as a double


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



function ebGammaAlpha_Callback(hObject, eventdata, handles)
% hObject    handle to ebGammaAlpha (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebGammaAlpha as text
%        str2double(get(hObject,'String')) returns contents of ebGammaAlpha as a double


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



function ebNCycles_Callback(hObject, eventdata, handles)
% hObject    handle to ebNCycles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebNCycles as text
%        str2double(get(hObject,'String')) returns contents of ebNCycles as a double


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


% ==============================================================
function ebAnalysisName_Callback(hObject, eventdata, handles)
tmp = get(handles.ebAnalysisName,'string');
if(isempty(tmp))
  errordlg('You cannot set Analysis Name to be empty');  
  set(handles.ebAnalysisName,'string',handles.flac.name);
  return;
end
if(~isempty(find(tmp,' ')))
  errordlg('Analysis Name cannot have blanks');  
  set(handles.ebAnalysisName,'string',handles.flac.name);
  return;
end
handles.flac.name = tmp;
setstate(handles);
guidata(hObject, handles);
return;

% --- Executes during object creation, after setting all properties.
function ebAnalysisName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebAnalysisName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% ==============================================================
function cbPeriodicDesign_Callback(hObject, eventdata, handles)
if(get(handles.cbPeriodicDesign,'value') == 1)
  handles.flac.ana.designtype = 'abblocked';  
else
  handles.flac.ana.designtype = 'event-related';
end
setstate(handles);
guidata(hObject, handles);
return;

% ==============================================================
function pbSave_Callback(hObject, eventdata, handles)
return;

% =======================================================
function cbMCExtReg_Callback(hObject, eventdata, handles)
val = get(handles.cbMCExtReg,'value');
if(val) handles.flac.extreg = 'mcextreg';
else    handles.flac.extreg = '';
end
setstate(handles);
guidata(hObject, handles);
return;

% =======================================================
function cbINorm_Callback(hObject, eventdata, handles)
val = get(cbINorm,'value');
if(val) handles.flac.inorm = 100;
else    handles.flac.inorm = 0;
end
setstate(handles);
guidata(hObject, handles);
return;

% --- Executes on button press in pbLoad.
function pbLoad_Callback(hObject, eventdata, handles)
% hObject    handle to pbLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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

%--------------------------------------------------%
function setstate(handles)

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

set(handles.ebFSD,'string',flac.fsd);
set(handles.ebRLF,'string',flac.runlistfile);
set(handles.ebTR,'string',flac.TR);
set(handles.ebTPExcludeFile,'string',flac.tpexcfile);
if(flac.inorm ~= 0) set(handles.cbINorm,'value',1);
else                set(handles.cbINorm,'value',0);
end
if(strcmp(ana.extreg,'mcextreg')) set(handles.cbMCExtReg,'value',1);
else                              set(handles.cbMCExtReg,'value',0);
end

tmp = sprintf('PolyFit %d',ana.PolyOrder);
set(handles.txPolyFit,'string',tmp);
set(handles.slPolyFit,'value',ana.PolyOrder);

set(handles.ebFIRTotTimeWin,'string',flac.ana.timewindow);
set(handles.ebFIRPreStim,'string',flac.ana.prestim);
if(flac.TR > 0 & flac.ana.TER < 0) flac.ana.TER = flac.ana.TR; end
set(handles.ebFIRTER,'string',flac.ana.TER);

if(strcmp(ana.designtype,'event-related') | ...
   strcmp(ana.designtype,'blocked'))
  EnableERBlock = 1;
  EnableABBlock = 0;
  % Enable ERBlock ----------------------
  set(handles.cbERBlock,'value',1);
  set(handles.txNConditions,'enable','on');
  set(handles.ebNConditions,'enable','on');
  set(handles.txParFile,'enable','on');
  set(handles.ebParFile,'enable','on');
  set(handles.cbHRFFIR,'enable','on');
  set(handles.cbHRFGamma,'enable','on');
  set(handles.cbHRFSPMHRF,'enable','on');
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
  set(handles.ebNConditions,'enable','off');
  set(handles.txParFile,'enable','off');
  set(handles.ebParFile,'enable','off');
  set(handles.cbHRFFIR,'enable','off');
  set(handles.cbHRFGamma,'enable','off');
  set(handles.cbHRFSPMHRF,'enable','off');
  % Enable ABBlocked ------------------
  set(handles.cbPeriodicDesign,'value',1);
  set(handles.ebNCycles,'enable','on');
  set(handles.txNCycles,'enable','on');
end

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

if(EnableERBlock && isempty(handles.flac.par4file))
  set(handles.ebParFile,'backgroundcolor','red')
else
  set(handles.ebParFile,'backgroundcolor','white')
end

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
set(handles.ebSPMHRFNDeriv,'enable',SPMHRFEnable);

if(~isempty(handles.clone))
  tmp = sprintf('Cloned from: %s\n',handles.clone);
  set(handles.txClone,'string',tmp);
end

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
flac.inorm = 100;
flac.TR = '';
flac.par4file = '';
flac.acfbins = 10; % set to 0 to turn off whitening
flac.stimulusdelay = 0;
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

return;


% --- Executes on slider movement.
function slPolyFit_Callback(hObject, eventdata, handles)
% hObject    handle to slPolyFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.flac.ana.PolyOrder = get(hObject,'Value');
setstate(handles);
guidata(hObject, handles);
return;

% --- Executes during object creation, after setting all properties.
function slPolyFit_CreateFcn(hObject, eventdata, handles)
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
return;

