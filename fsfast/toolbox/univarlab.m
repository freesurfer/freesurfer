function varargout = univarlab(varargin)
% UNIVARLAB Application M-file for univarlab.fig
%    FIG = UNIVARLAB launch univarlab GUI.
%    UNIVARLAB('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 07-Aug-2003 02:41:50

if(nargin == 0)  
  % LAUNCH GUI
  fig = openfig(mfilename,'new');
  % Generate a structure of handles to pass to callbacks, and store it. 
  handles = guihandles(fig); % only run this once, later use h = guidata
  handles.hunivarlab = fig;
  handles = InitGUIParams(handles);
  handles = Synth(handles);
  handles = PlotRaw(handles);
guidata(fig, handles);
  if(nargout > 0) varargout{1} = fig; end
elseif(ischar(varargin{1})) 
  % INVOKE NAMED SUBFUNCTION OR CALLBACK
  try
    if(nargout) [varargout{1:nargout}] = feval(varargin{:}); 
    else feval(varargin{:}); 
    end
  catch
    disp(lasterr);
  end
end

%| ABOUT CALLBACKS:
%| GUIDE automatically appends subfunction prototypes to this file, and 
%| sets objects' callback properties to call them through the FEVAL 
%| switchyard above. This comment describes that mechanism.
%|
%| Each callback subfunction declaration has the following form:
%| <SUBFUNCTION_NAME>(H, EVENTDATA, HANDLES, VARARGIN)
%|
%| The subfunction name is composed using the object's Tag and the 
%| callback type separated by '_', e.g. 'slider2_Callback',
%| 'figure1_CloseRequestFcn', 'axis1_ButtondownFcn'.
%|
%| H is the callback object's handle (obtained using GCBO).
%|
%| EVENTDATA is empty, but reserved for future use.
%|
%| HANDLES is a structure containing handles of components in GUI using
%| tags as fieldnames, e.g. handles.figure1, handles.slider2. This
%| structure is created at GUI startup using GUIHANDLES and stored in
%| the figure's application data using GUIDATA. A copy of the structure
%| is passed to each callback.  You can store additional information in
%| this structure at GUI startup, and you can change the structure
%| during callbacks.  Call guidata(h, handles) after changing your
%| copy to replace the stored original so that subsequent callbacks see
%| the updates. Type "help guihandles" and "help guidata" for more
%| information.
%|
%| VARARGIN contains any extra arguments you have passed to the
%| callback. Specify the extra arguments by editing the callback
%| property in the inspector. By default, GUIDE sets the property to:
%| <MFILENAME>('<SUBFUNCTION_NAME>', gcbo, [], guidata(gcbo))
%| Add any extra arguments after the last argument, before the final
%| closing parenthesis.



% --------------------------------------------------------------------
function varargout = quit_pb_Callback(h, eventdata, gd, varargin)




% --------------------------------------------------------------------
function varargout = reset_pb_Callback(h, eventdata, gd, varargin)

% --------------------------------------------------------------------
function varargout = synth_pb_Callback(h, eventdata, gd, varargin)
gd = Synth(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
%univarlab('PlotRaw',gd);
%guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = GamDeltaSynth_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.GamDeltaSynth(2) | val > gd.GamDeltaSynth(3))
  set(h,'string',sprintf('%g',gd.GamDeltaSynth(1)));
  return;
end
gd.GamDeltaSynth(1) = val;
gd = MakeFirstLevelAnalysis(gd);
gd = Synth(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = GamTauSynth_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.GamTauSynth(2) | val > gd.GamTauSynth(3))
  set(h,'string',sprintf('%g',gd.GamTauSynth(1)));
  return;
end
gd.GamTauSynth(1) = val;
gd = MakeFirstLevelAnalysis(gd);
gd = Synth(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = C1AmpSynth_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.C1AmpSynth(2) | val > gd.C1AmpSynth(3))
  set(h,'string',sprintf('%g',gd.C1AmpSynth(1)));
  return;
end
gd.C1AmpSynth(1) = val;
gd = Synth(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = C2AmpSynth_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.C2AmpSynth(2) | val > gd.C2AmpSynth(3))
  set(h,'string',sprintf('%g',gd.C2AmpSynth(1)));
  return;
end
gd.C2AmpSynth(1) = val;
gd = Synth(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = NoiseStd_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.NoiseStd(2) | val > gd.NoiseStd(3))
  set(h,'string',sprintf('%g',gd.NoiseStd(1)));
  return;
end
gd.NoiseStd(1) = val;
gd = Synth(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = NoiseAR1_Callback(h, eventdata, gd, varargin)




% --------------------------------------------------------------------
function varargout = FIRPreWin_Callback(h, eventdata, gd, varargin)




% --------------------------------------------------------------------
function varargout = FIRPostWin_Callback(h, eventdata, gd, varargin)




% --------------------------------------------------------------------
function varargout = GamDeltaEst_Callback(h, eventdata, gd, varargin)




% --------------------------------------------------------------------
function varargout = GamTauEst_Callback(h, eventdata, gd, varargin)




% --------------------------------------------------------------------
function varargout = C1AmpMan_Callback(h, eventdata, gd, varargin)




% --------------------------------------------------------------------
function varargout = C2AmpMan_Callback(h, eventdata, gd, varargin)




% --------------------------------------------------------------------
function varargout = Whiten_cb_Callback(h, eventdata, gd, varargin)

% --------------------------------------------------------------------
function varargout = ViewFIRDesign_Callback(h, eventdata, gd, varargin)
if(isempty(gd.hXFIR) | ~ishandle(gd.hXFIR)) gd.hXFIR = figure; end
figure(gd.hXFIR);
imagesc(gd.XFIR);
colorbar;
title('FIR Design Matrix: nrows = ntimepoints');
figure(gd.hunivarlab);
guidata(gd.hunivarlab,gd)
return;

% --------------------------------------------------------------------
function varargout = ViewGamDesign_Callback(h, eventdata, gd, varargin)
if(isempty(gd.hXGam) | ~ishandle(gd.hXGam)) gd.hXGam = figure; end
figure(gd.hXGam);
imagesc(gd.XGam);
colorbar;
title('Gamma Design Matrix: nrows = ntimepoints');
figure(gd.hunivarlab);
guidata(gd.hunivarlab,gd)
return;

%-------------------------------------------------%
function gd = PlotRaw(gd)
if(isempty(gd.hRaw) | ~ishandle(gd.hRaw)) gd.hRaw = figure; end
figure(gd.hRaw);
plot(gd.t,gd.yobserved, gd.t,gd.ysignal,'+-', gd.t,gd.yhatfir, ...
     gd.t,gd.yhatgam);
title('Raw Time Courses');
legend('Observed','True','FIR','Gamma');
figure(gd.hunivarlab);
guidata(gd.hunivarlab,gd)
return;

% --------------------------------------------------------------------
function gd = InitGUIParams(gd)
% ---------- synthesis parameters ------------------%
gd.TR = 2;
gd.ntp = 120;
gd.GamDeltaSynth = [2.25  0.00 5.00];
gd.GamTauSynth   = [1.25  0.10 2.00];
gd.C1AmpSynth    = [1.00 -2.00 2.00];
gd.C2AmpSynth    = [1.00 -2.00 2.00];
gd.schedule      = [];   % stimulus schedule
gd.GamSynthFLA   = [];   % First level analysis struct for synth
gd.NoiseStd      = [1.00  0.00 2.00];
gd.NoiseAR1      = [0.00  0.00 1.00];

% ---------- FIR Estimation parameters ------------------%
gd.FIRPSDMin     = [-4.00 -8.00  0.00];
gd.FIRPSDMax     = [20.00  4.00 40.00];
gd.FIRTotWin     = gd.FIRPreWin(1) + gd.FIRPostWin(1);
gd.FIRRStd       = 0.0; % Estimated value of Fir resid std
gd.FIREstFLA        = [];  % First level analysis struct for est FIR

% ---------- Gamma Estimation parameters ------------------%
gd.GamDeltaEst = gd.GamDeltaSynth; % Use this val in estimation model
gd.GamTauEst   = gd.GamTauSynth;   % Use this val in estimation model
gd.C1AmpEst    = 0.0;  % Estimated value of C1
gd.C2AmpEst    = 0.0;  % Estimated value of C2
gd.GamRStd     = 0.0;  % Estimated value of Gamma resid std
gd.GamEstFLA   = [];   % First level analysis struct for est gamma

% ---------- Manual-Gamma Estimation parameters ------------------%
gd.C1AmpMan    = 0.0;  % Use this val to construct manual signal
gd.C2AmpMan    = 0.0;  % Use this val to construct manual signal
gd.ManGamRStd  = 0.0;  % Estimated value of ManGamma resid std

gd.Whiten      = 0;    % Whitening flag (0 or 1)

% ---------- Raw Synthesized Vectors -----------------------%
gd.t = gd.TR*[0:gd.ntp-1]';
gd.ypar = [];
gd.ynoise = [];
gd.yc1signal = [];
gd.yc2signal = [];
gd.ysignal = [];
gd.yobserved = [];

% ---------- Raw Estimation Vectors -----------------------%
gd.yhatfir = [];
gd.resfir = [];
gd.acffir = [];
gd.yhatgam = [];
gd.resgam = [];
gd.acfgam = [];
gd.yhatmangam = [];
gd.resmangam = [];
gd.acfmangam = [];

gd = NewStimSched(gd);
gd = MakeFirstLevelAnalysis(gd);

gd.hXFIR = [];
gd.hXGam = [];
gd.hRaw = [];

return;
% --------------------------------------------------------------------
function gd = NewStimSched(gd)

nc1 = round(gd.ntp/3);
nc2 = round(gd.ntp/3);
n0 = gd.ntp - (nc1+nc2);

r = [zeros(1,n0) ones(1,nc1) 2*ones(1,nc2)]';
gd.schedule = [gd.t r(randperm(gd.ntp))];

return;
% --------------------------------------------------------------------
function gd = MakeFirstLevelAnalysis(gd)

%-------- Est FIR ------------------------%
gd.FIREstFLA = fast_flacfg_struct;
gd.FIREstFLA.flaname = 'FIR';
gd.FIREstFLA.TR      = gd.TR;

gd.FIREstFLA.sesscfg           = fast_sesscfg_struct;
gd.FIREstFLA.sesscfg.ntp       = gd.ntp;
gd.FIREstFLA.sesscfg.evschlist(1).evsch = gd.schedule;
gd.FIREstFLA.sesscfg.runlist   = '001';
gd.FIREstFLA.nthrun            = 1;

fxline = sprintf('effect fixed cond1 fir 1 %g %g %g 0',...
		       gd.FIRPSDMin(1),gd.TR,gd.FIRPSDMax(1));
gd.FIREstFLA.fxlist(1).fx = fast_fxcfg('parseline',fxline);
fxline = sprintf('effect fixed cond2 fir 2 %g %g %g 0',...
		       gd.FIRPSDMin(1), gd.TR, gd.FIRPSDMax(1));
gd.FIREstFLA.fxlist(2).fx = fast_fxcfg('parseline',fxline);
gd.XFIR = fast_fla_desmat(gd.FIREstFLA);

%-------- Est Gamma ------------------------%
gd.GamEstFLA = gd.FIREstFLA;
fxline = sprintf('effect fixed cond1 gamma 1 0 .1 30 0 0 %g %g 0',...
		 gd.GamDeltaEst(1),gd.GamTauEst(1));
gd.GamEstFLA.fxlist(1).fx = fast_fxcfg('parseline',fxline);
fxline = sprintf('effect fixed cond2 gamma 2 0 .1 30 0 0 %g %g 0',...
		 gd.GamDeltaEst(1),gd.GamTauEst(1));
gd.GamEstFLA.fxlist(2).fx = fast_fxcfg('parseline',fxline);
gd.XGam = fast_fla_desmat(gd.GamEstFLA);

%-------- Synth Gamma ------------------------%
gd.GamSynthFLA = gd.FIREstFLA;
fxline = sprintf('effect fixed cond1 gamma 1 0 .1 30 0 0 %g %g 0',...
		 gd.GamDeltaSynth(1),gd.GamTauSynth(1));
gd.GamSynthFLA.fxlist(1).fx = fast_fxcfg('parseline',fxline);
fxline = sprintf('effect fixed cond2 gamma 2 0 .1 30 0 0 %g %g 0',...
		 gd.GamDeltaSynth(1),gd.GamTauSynth(1));
gd.GamSynthFLA.fxlist(2).fx = fast_fxcfg('parseline',fxline);
gd.XGamSynth = fast_fla_desmat(gd.GamSynthFLA);

return;
% --------------------------------------------------------------------

function gd = Synth(gd)

gd.ynoise    = gd.NoiseStd(1)*randn(gd.ntp,1);
gd.yc1signal = gd.XGamSynth(:,1) * gd.C1AmpSynth(1);
gd.yc2signal = gd.XGamSynth(:,2) * gd.C2AmpSynth(1);
gd.ysignal   = gd.yc1signal + gd.yc2signal;
gd.yobserved = gd.ysignal + gd.ynoise;

gd.FIRbeta = (inv(gd.XFIR'*gd.XFIR)*gd.XFIR')*gd.yobserved;
gd.yhatfir = gd.XFIR * gd.FIRbeta;
gd.resfir  = gd.yobserved - gd.yhatfir;

gd.Gambeta = (inv(gd.XGam'*gd.XGam)*gd.XGam')*gd.yobserved;
gd.yhatgam = gd.XGam * gd.Gambeta;
gd.resgam  = gd.yobserved - gd.yhatgam;

return;








