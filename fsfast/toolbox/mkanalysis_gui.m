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

% Last Modified by GUIDE v2.5 10-Jun-2007 18:10:31

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
handles.version = '$Id: mkanalysis_gui.m,v 1.3 2007/06/11 03:08:43 greve Exp $';

%handles = parse_args(handles,varargin);
if(~isempty(MkAnalysisName))
  set(handles.ebAnalysisName,'string',MkAnalysisName);
  if(exist(MkAnalysisName,'dir'))
    handles.flac = fast_ldanaflac(MkAnalysisName);
    setstate(handles);
  end
end


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



function ebFuncStem_Callback(hObject, eventdata, handles)
% hObject    handle to ebFuncStem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebFuncStem as text
%        str2double(get(hObject,'String')) returns contents of ebFuncStem as a double


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



function ebTR_Callback(hObject, eventdata, handles)
% hObject    handle to ebTR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebTR as text
%        str2double(get(hObject,'String')) returns contents of ebTR as a double


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


% --- Executes on button press in cbERBlock.
function cbERBlock_Callback(hObject, eventdata, handles)
% hObject    handle to cbERBlock (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: get(hObject,'Value') returns toggle state of cbERBlock
if(get(handles.cbERBlock,'value') == 1)
  handles.flac.ana.designtype = 'event-related';
else
  handles.flac.ana.designtype = 'abblocked';  
end
setstate(handles);
guidata(hObject, handles);

function ebNConditions_Callback(hObject, eventdata, handles)
% hObject    handle to ebNConditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebNConditions as text
%        str2double(get(hObject,'String')) returns contents of ebNConditions as a double


% --- Executes during object creation, after setting all properties.
function ebNConditions_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ebNConditions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cbHRFGamma.
function cbHRFGamma_Callback(hObject, eventdata, handles)
% hObject    handle to cbHRFGamma (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbHRFGamma


% --- Executes on button press in cbHRFFIR.
function cbHRFFIR_Callback(hObject, eventdata, handles)
% hObject    handle to cbHRFFIR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbHRFFIR


% --- Executes on button press in cbHRFSPMHRF.
function cbHRFSPMHRF_Callback(hObject, eventdata, handles)
% hObject    handle to cbHRFSPMHRF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbHRFSPMHRF



function ebFIRTotTimeWin_Callback(hObject, eventdata, handles)
% hObject    handle to ebFIRTotTimeWin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebFIRTotTimeWin as text
%        str2double(get(hObject,'String')) returns contents of ebFIRTotTimeWin as a double


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



function ebAnalysisName_Callback(hObject, eventdata, handles)
% hObject    handle to ebAnalysisName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ebAnalysisName as text
%        str2double(get(hObject,'String')) returns contents of ebAnalysisName as a double


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


% --- Executes on button press in cbPeriodicDesign.
function cbPeriodicDesign_Callback(hObject, eventdata, handles)
% hObject    handle to cbPeriodicDesign (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(handles.cbPeriodicDesign,'value') == 1)
  handles.flac.ana.designtype = 'abblocked';  
else
  handles.flac.ana.designtype = 'event-related';
end
setstate(handles);
guidata(hObject, handles);


% --- Executes on button press in pbSave.
function pbSave_Callback(hObject, eventdata, handles)
% hObject    handle to pbSave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in pmPolyFit.
function pmPolyFit_Callback(hObject, eventdata, handles)
% hObject    handle to pmPolyFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns pmPolyFit contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pmPolyFit


% --- Executes during object creation, after setting all properties.
function pmPolyFit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pmPolyFit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in cbMCExtReg.
function cbMCExtReg_Callback(hObject, eventdata, handles)
% hObject    handle to cbMCExtReg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbMCExtReg


% --- Executes on button press in cbINorm.
function cbINorm_Callback(hObject, eventdata, handles)
% hObject    handle to cbINorm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of cbINorm


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
set(handles.ebFSD,'string',flac.fsd);
set(handles.ebRLF,'string',flac.runlistfile);
set(handles.ebTR,'string',flac.TR);
set(handles.ebTPExcludeFile,'string',flac.tpexcfile);
if(flac.inorm ~= 0) set(handles.cbINorm,'value',1);
else                set(handles.cbINorm,'value',0);
end
if(strcmp(ana.extreg,'mcext')) set(handles.cbMCExtReg,'value',1);
else                           set(handles.cbMCExtReg,'value',1);
end

if(strcmp(ana.designtype,'event-related') | ...
   strcmp(ana.designtype,'blocked'))
  % Enable ERBlock ----------------------
  set(handles.cbERBlock,'value',1);
  set(handles.txNConditions,'enable','on');
  set(handles.ebNConditions,'enable','on');
  set(handles.cbHRFFIR,'enable','on');
  set(handles.cbHRFGamma,'enable','on');
  set(handles.cbHRFSPMHRF,'enable','on');

  if(get(handles.cbHRFFIR,'value') == 1)
    FIRenable   = 'on';
    GammaEnable = 'off';
    SPMHRFEnable = 'off';
  end
  if(get(handles.cbHRFGamma,'value') == 1)
    FIRenable   = 'off';
    GammaEnable = 'on';
    SPMHRFEnable = 'off';
  end
  if(get(handles.cbHRFSPMHRF,'value') == 1)
    FIRenable    = 'off';
    GammaEnable  = 'off';
    SPMHRFEnable = 'on';
  end
  
  % Disable ABBlocked ------------------
  set(handles.cbPeriodicDesign,'value',0);
  set(handles.ebNCycles,'enable','off');
  set(handles.txNCycles,'enable','off');
else
  % Enable ABBlocked ------------------
  set(handles.cbPeriodicDesign,'value',1);
  set(handles.ebNCycles,'enable','on');
  set(handles.txNCycles,'enable','on');
  % Disable ERBlock ----------------------
  set(handles.cbERBlock,'value',0);
  set(handles.txNConditions,'enable','off');
  set(handles.ebNConditions,'enable','off');

  set(handles.cbHRFFIR,'enable','off');
  set(handles.txFIRTotTimeWin,'enable','off');
  set(handles.ebFIRTotTimeWin,'enable','off');
  set(handles.txFIRPreStim,'enable','off');
  set(handles.ebFIRPreStim,'enable','off');
  set(handles.txFIRTER,'enable','off');
  set(handles.ebFIRTER,'enable','off');

  set(handles.cbHRFGamma,'enable','off');
  set(handles.txGammaDelta,'enable','off');
  set(handles.ebGammaDelta,'enable','off');
  set(handles.txGammaTau,'enable','off');
  set(handles.ebGammaTau,'enable','off');
  set(handles.txGammaAlpha,'enable','off');
  set(handles.ebGammaAlpha,'enable','off');

  set(handles.cbHRFSPMHRF,'enable','off');
  set(handles.txSPMHRFNDeriv,'enable','off');
  set(handles.ebSPMHRFNDeriv,'enable','off');
  
end

set(handles.txFIRTotTimeWin,'enable','on');
set(handles.ebFIRTotTimeWin,'enable','on');
set(handles.txFIRPreStim,'enable','on');
set(handles.ebFIRPreStim,'enable','on');
set(handles.txFIRTER,'enable','on');
set(handles.ebFIRTER,'enable','on');

set(handles.txGammaDelta,'enable','on');
set(handles.ebGammaDelta,'enable','on');
set(handles.txGammaTau,'enable','on');
set(handles.ebGammaTau,'enable','on');
set(handles.txGammaAlpha,'enable','on');
set(handles.ebGammaAlpha,'enable','on');

set(handles.txSPMHRFNDeriv,'enable','on');
set(handles.ebSPMHRFNDeriv,'enable','on');






return;
