function varargout = univarlab(varargin)
% UNIVARLAB Application M-file for univarlab.fig
%    FIG = UNIVARLAB launch univarlab GUI.
%    UNIVARLAB('callback_name', ...) invoke the named callback.
% Last Modified by GUIDE v2.0 20-Aug-2003 23:16:15
%

% To do:
%   Closing/Quiting functions
%   Choose whether to display some traces
%   Better HDR
%   Info/Help
%   TR/Ntp
%   Signifiances
%   Slider for Manual
%   Test for Matrix Condition
%   BoxCar Width
%   Fix FixACF
%   Seed
%   Gamma First Derivative
%   Slice Timing Correction
%   Save all schedules, re-gen only when asked
%   Separate synthesis model for both conditions


if(nargin == 0)  
  % LAUNCH GUI
  fig = openfig(mfilename,'new');
  % Generate a structure of handles to pass to callbacks, and store it. 
  gd = guihandles(fig); % only run this once, later use h = guidata
  gd.hunivarlab = fig;
  gd = InitGUIParams(gd);
  gd = SynthNoise(gd);
  gd = ScaleFilterNoise(gd);
  gd = NewStimSched(gd);
  gd = MakeFirstLevelAnalysis(gd);
  gd = SynthC1Signal(gd);
  gd = SynthC2Signal(gd);
  gd = SynthSignal(gd);
  gd = SynthObserved(gd);
  gd = Estimate(gd);
  gd = PlotRaw(gd);
  gd = PlotACF(gd);
  gd = PlotHRF(gd);
  guidata(fig, gd);
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
return;
%----------------- end of main ---------------------------%
%----------------- end of main ---------------------------%
%----------------- end of main ---------------------------%

% --------------------------------------------------------------------
function gd = InitGUIParams(gd)
  % ---------- synthesis parameters ------------------%
gd.TR = 2;
gd.ntp = 120;
gd.GamDeltaSynth = [2.25  0.00 5.00];
gd.GamTauSynth   = [1.25  0.10 2.00];
gd.C1AmpSynth    = [1.00 -2.00 2.00];
gd.C2AmpSynth    = [1.00 -2.00 2.00];
gd.ScheduleType  = 'rper';
gd.schedule      = [];   % stimulus schedule, [t, id]
gd.GamSynthFLA   = [];   % First level analysis struct for synth
gd.NoiseStd      = [1.00  0.00 2.00];
gd.NoiseAR1      = [0 -0.9 0.9];
gd = PrepNoiseACF(gd);
gd.FixACF = 0;
gd.FIRFixACFMtx = [];
gd.GamFixACFMtx = [];
gd.ACFShowLagZero = 1;

  % ---------- FIR Estimation parameters ------------------%
gd.FIRPSDMin     = [-4.00 -8.00  0.00];
gd.FIRPSDMax     = [20.00  4.00 40.00];
gd.FIRTotWin     = gd.FIRPSDMax(1) - gd.FIRPSDMin(1);
gd.FIRRStd       = 0.0; % Estimated value of Fir resid std
gd.FIREstFLA        = [];  % First level analysis struct for est FIR

  % ---------- Gamma Estimation parameters ------------------%
gd.GamDeltaEst = gd.GamDeltaSynth; % Use this val in estimation model
gd.GamTauEst   = gd.GamTauSynth;   % Use this val in estimation model
gd.GamNDeriv   = 0;
gd.C1AmpEst    = 0.0;  % Estimated value of C1
gd.C2AmpEst    = 0.0;  % Estimated value of C2
gd.GamRStd     = 0.0;  % Estimated value of Gamma resid std
gd.GamEstFLA   = [];   % First level analysis struct for est gamma

  % ---------- Manual-Gamma Estimation parameters ------------------%
gd.C1AmpMan    = [1.0 -2 2];  % Use this val to construct manual signal
gd.C2AmpMan    = [1.0 -2 2];  % Use this val to construct manual signal

gd.Whiten      = 0;    % Whitening flag (0 or 1)
gd.FitMean     = 0;    % Include regressor for MeanOffset

  % ---------- Raw Synthesized Vectors -----------------------%
gd.t = gd.TR*[0:gd.ntp-1]';
gd.ypar = [];
gd.ynoise0 = []; % unscaled
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
gd.yhatmangam = []; % Manual gamma
gd.resmangam = [];  % Manual gamma

gd.hXFIR = [];
gd.hXGam = [];
gd = SetGUIParams(gd);

gd.ShowRaw = 0;
gd.ShowRawSignal = 1;
gd.ShowRawObserved = 1;
gd.ShowACF = 0;
gd.ShowHRF = 1;

return;

% --------------------------------------------------------------------
function varargout = quit_pb_Callback(h, eventdata, gd, varargin)

% --------------------------------------------------------------------
function varargout = reset_pb_Callback(h, eventdata, gd, varargin)
gd = InitGUIParams(gd);
gd = SynthNoise(gd);
gd = ScaleFilterNoise(gd);
gd = NewStimSched(gd);
gd = MakeFirstLevelAnalysis(gd);
gd = SynthC1Signal(gd);
gd = SynthC2Signal(gd);
gd = SynthSignal(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab, gd);
return;
  
% --------------------------------------------------------------------
% Synthesize a new schedule %
function varargout = Schedule_pb_Callback(h, eventdata, gd, varargin)
gd = NewStimSched(gd);
gd = MakeFirstLevelAnalysis(gd);
gd = SynthC1Signal(gd);
gd = SynthC2Signal(gd);
gd = SynthSignal(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = Schedule_pum_Callback(h, eventdata, gd, varargin)
s = get(h,'string');
v = get(h,'value');
scheduletype = lower(char(s(v,:)));
%scheduletype = lower(deblank(s(v,:)));
if(strcmp(scheduletype,gd.ScheduleType)) return; end
gd.ScheduleType = scheduletype;
gd = NewStimSched(gd);
gd = MakeFirstLevelAnalysis(gd);
gd = SynthC1Signal(gd);
gd = SynthC2Signal(gd);
gd = SynthSignal(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);

return;

% --------------------------------------------------------------------
% Synthesize new noise %
function varargout = synth_pb_Callback(h, eventdata, gd, varargin)
gd = SynthNoise(gd);
gd = ScaleFilterNoise(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);
%univarlab('PlotRaw',gd);
%guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
% Change Delta parameter used in Synthesis Gamma Model
function varargout = GamDeltaSynth_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.GamDeltaSynth(2) | val > gd.GamDeltaSynth(3))
  set(h,'string',sprintf('%g',gd.GamDeltaSynth(1)));
  return;
end
gd.GamDeltaSynth(1) = val;
gd = MakeFirstLevelAnalysis(gd);
gd = SynthC1Signal(gd);
gd = SynthC2Signal(gd);
gd = SynthSignal(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
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
gd = SynthC1Signal(gd);
gd = SynthC2Signal(gd);
gd = SynthSignal(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
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
gd = SynthC1Signal(gd);
gd = SynthSignal(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
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
gd = SynthC2Signal(gd);
gd = SynthSignal(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
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
gd = ScaleFilterNoise(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = NoiseAR1_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.NoiseAR1(2) | val > gd.NoiseAR1(3))
  set(h,'string',sprintf('%g',gd.NoiseAR1(1)));
  return;
end
gd.NoiseAR1(1) = val;
gd = PrepNoiseACF(gd);
gd = ScaleFilterNoise(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = FIRPSDMin_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.FIRPSDMin(2) | val > gd.FIRPSDMin(3))
  set(h,'string',sprintf('%g',gd.FIRPSDMin(1)));
  return;
end
val = gd.TR*round(val/gd.TR);
set(h,'string',sprintf('%g',val));
gd.FIRPSDMin(1) = val;
gd.FIRTotWin  = gd.FIRPSDMax(1) - gd.FIRPSDMin(1);
set(gd.FIRTotWinTxt,'string',sprintf('%g',gd.FIRTotWin));
gd = MakeFirstLevelAnalysis(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = FIRPSDMax_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
val = gd.TR*round(val/gd.TR);
if(val < gd.FIRPSDMax(2) | val > gd.FIRPSDMax(3))
  set(h,'string',sprintf('%g',gd.FIRPSDMax(1)));
  return;
end
set(h,'string',sprintf('%g',val));
if(abs(val-gd.FIRPSDMax(1))==0) return; end
gd.FIRPSDMax(1) = val;
gd.FIRTotWin  = gd.FIRPSDMax(1) - gd.FIRPSDMin(1);
set(gd.FIRTotWinTxt,'string',sprintf('%g',gd.FIRTotWin));
gd = MakeFirstLevelAnalysis(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = GamDeltaEst_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.GamDeltaEst(2) | val > gd.GamDeltaEst(3))
  set(h,'string',sprintf('%g',gd.GamDeltaEst(1)));
  return;
end
gd.GamDeltaEst(1) = val;
gd = MakeFirstLevelAnalysis(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);

% --------------------------------------------------------------------
function varargout = GamTauEst_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.GamTauEst(2) | val > gd.GamTauEst(3))
  set(h,'string',sprintf('%g',gd.GamTauEst(1)));
  return;
end
gd.GamTauEst(1) = val;
gd = MakeFirstLevelAnalysis(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd);

% --------------------------------------------------------------------
function varargout = C1AmpMan_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.C1AmpMan(2) | val > gd.C1AmpMan(3))
  set(h,'string',sprintf('%g',gd.C1AmpMan(1)));
  return;
end
gd.C1AmpMan(1) = val;
set(gd.C1ManAmp_sl,'value',gd.C1AmpMan(1));
gd = EstimateManual(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = C1ManAmp_sl_Callback(h, eventdata, gd, varargin)
val = get(h,'value');
gd.C1AmpMan(1) = val;
set(gd.C1AmpMan_et,'string',sprintf('%g',gd.C1AmpMan(1)));
gd = EstimateManual(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = C2AmpMan_Callback(h, eventdata, gd, varargin)
val = sscanf(get(h,'string'),'%f');
if(val < gd.C2AmpMan(2) | val > gd.C2AmpMan(3))
  set(h,'string',sprintf('%g',gd.C2AmpMan(1)));
  return;
end
gd.C2AmpMan(1) = val;
set(gd.C2ManAmp_sl,'value',gd.C1AmpMan(2));
gd = EstimateManual(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;

% --------------------------------------------------------------------
function varargout = C2ManAmp_sl_Callback(h, eventdata, gd, varargin)
val = get(h,'value');
gd.C2AmpMan(1) = val;
set(gd.C2AmpMan_et,'string',sprintf('%g',gd.C2AmpMan(1)));
gd = EstimateManual(gd);
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd);
return;


% --------------------------------------------------------------------
function varargout = Whiten_cb_Callback(h, eventdata, gd, varargin)
gd.Whiten = ~gd.Whiten;
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
return;

% --------------------------------------------------------------------
function varargout = FixACF_cb_Callback(h, eventdata, gd, varargin)
gd.FixACF = ~gd.FixACF;
gd = PlotACF(gd);
return;

% --------------------------------------------------------------------
function varargout = FitMean_cb_Callback(h, eventdata, gd, varargin)
gd.FitMean = ~gd.FitMean; 
gd = MakeFirstLevelAnalysis(gd);
gd = SynthC1Signal(gd);
gd = SynthC2Signal(gd);
gd = SynthSignal(gd);
gd = SynthObserved(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab, gd);
return;

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

% --------------------------------------------------------------------
function varargout = GamDeriv_cb_Callback(h, eventdata, gd, varargin)
gd.GamNDeriv = ~gd.GamNDeriv;
gd = MakeFirstLevelAnalysis(gd);
gd = Estimate(gd);
gd = PlotRaw(gd);
gd = PlotACF(gd);
gd = PlotHRF(gd);
guidata(gd.hunivarlab, gd);

return;


%-------------------------------------------------%
function gd = PlotHRF(gd)
if(gd.ShowHRF==0) return; end
if(~isfield(gd,'hHRF') | isempty(gd.hHRF) | ~ishandle(gd.hHRF)) 
  gd.hHRF = figure; 
  set(gcf,'DeleteFcn','univarlab(''DeleteHRF'')');
  ud.hunivarlab = gd.hunivarlab;
  set(gd.hHRF,'userdata',ud);
end
figure(gd.hHRF);

gd.GamSynthFLA.nthfx = 1;
GamSynthPSD = fast_fxcfg('irftaxis',gd.GamSynthFLA);
GamSynthHRF = fast_fxcfg('irfmatrix',gd.GamSynthFLA);
GamSynthHRF1 = GamSynthHRF * gd.C1AmpSynth(1);
GamSynthHRF2 = GamSynthHRF * gd.C2AmpSynth(1);

gd.FIREstFLA.nthfx = 1;
[FIRHRF1 FIRPSD] = fast_fla_irf(gd.FIREstFLA,gd.FIRbeta);
gd.FIREstFLA.nthfx = 2;
[FIRHRF2 FIRPSD] = fast_fla_irf(gd.FIREstFLA,gd.FIRbeta);

gd.GamEstFLA.nthfx = 1;
[GamEstHRF1 GamEstPSD] = fast_fla_irf(gd.GamEstFLA,gd.Gambeta);
gd.GamEstFLA.nthfx = 2;
[GamEstHRF2 GamEstPSD] = fast_fla_irf(gd.GamEstFLA,gd.Gambeta);

%FIRPSD = fast_fxcfg('irftaxis',gd.FIREstFLA);
%nPerFIR = length(FIRPSD);
%FIRHRF1 = gd.FIRbeta(1:nPerFIR);
%FIRHRF2 = gd.FIRbeta(nPerFIR+1:2*nPerFIR);

%gd.GamEstFLA.nthfx = 1;
%GamEstPSD = fast_fxcfg('irftaxis',gd.GamEstFLA);
%GamEstHRF = fast_fxcfg('irfmatrix',gd.GamEstFLA);
%GamEstHRF1 = GamEstHRF * gd.Gambeta(1:1+gd.GamNDeriv);
%GamEstHRF2 = GamEstHRF * gd.Gambeta(2+gd.GamNDeriv:2+2*gd.GamNDeriv);

plot(GamSynthPSD,GamSynthHRF1,'+-',...
     GamSynthPSD,GamSynthHRF2,'+-',...
     FIRPSD,FIRHRF1,FIRPSD,FIRHRF2,...
     GamEstPSD,GamEstHRF1,GamEstPSD,GamEstHRF2);
title('Hemodynamic Responses');
legend('Ideal1','Ideal2','FIR1','FIR2','Gam1','Gam2');
xlabel('post-stimulus delay (sec)');
figure(gd.hunivarlab);
guidata(gd.hunivarlab,gd)
return;

%-------------------------------------------------%
% Delete function for HRF window %
function DeleteHRF
ud = get(gcf,'userdata');
gd = guidata(ud.hunivarlab);
gd.ShowHRF = 0;
gd.hHRF = [];
guidata(gd.hunivarlab,gd)
return;

% --------------------------------------------------------------------
function varargout = ShowHRF_mi_Callback(h, eventdata, gd, varargin)
if(gd.ShowHRF)  figure(gd.hHRF); return; end
gd.ShowHRF = 1;
gd = PlotHRF(gd);
guidata(gd.hunivarlab,gd)
return;

%-------------------------------------------------%
function gd = PlotRaw(gd)
if(gd.ShowRaw == 0) return; end
if(~isfield(gd,'hRaw') | isempty(gd.hRaw) | ~ishandle(gd.hRaw)) 
  gd.hRaw = figure; 
  set(gcf,'Interruptible','off'); % smooth animations %
  set(gcf,'DoubleBuffer','on');   % smooth animations %
  set(gcf,'BusyAction','cancel'); % dont build up a lot of events %
  set(gcf,'renderer','painters'); % seems to be the best
  set(gcf,'DeleteFcn','univarlab(''DeleteRaw'')');
  ud.hunivarlab = gd.hunivarlab;
  set(gd.hRaw,'userdata',ud);
  ud.hToggleRawSignal = ...
      uicontrol('Style', 'checkbox', ...
		'String', 'Signal',...
		'Position', [1 1 50 25], ...
		'Callback', 'univarlab(''ToggleRawTrace'',''signal'')',...
		'tooltipstring','Toggle Display of Signal Raw Trace',...
		'value',gd.ShowRawSignal);
  ud.hToggleRawObserved = ...
      uicontrol('Style', 'checkbox', ...
		'String', 'Observed',...
		'Position', [1 26 50 25], ...
		'Callback', 'univarlab(''ToggleRawTrace'',''observed'')',...
		'tooltipstring','Toggle Display of Observed Trace',...
		'value',gd.ShowRawSignal);
  ud.hunivarlab = gd.hunivarlab;
  set(gd.hRaw,'userdata',ud);
end
figure(gd.hRaw);
tmp = [gd.yobserved gd.ysignal gd.yhatfir gd.yhatgam gd.yhatmangam];
ymin = min(reshape1d(tmp));
ymax = max(reshape1d(tmp));
yrange = ymax-ymin;
ind1 = find(gd.schedule(:,2)==1);
nind1 = length(ind1);
p1 = (ymin-.05*yrange)*ones(nind1,1);
ind2 = find(gd.schedule(:,2)==2);
nind2 = length(ind2);
p2 = (ymin-.05*yrange)*ones(nind2,1);

%p1 = gd.schedule(ind1,2) * (ymax-ymin)/2 + ymin;
%p2 = gd.schedule(ind2,2) * (ymax-ymin)/2 + ymin;

plot(gd.t,gd.ysignal,'+-', ...
     gd.t,gd.yobserved, 'o-', ...
     gd.t,gd.yhatfir, ...
     gd.t,gd.yhatgam, gd.t,gd.yhatmangam, ...
     gd.t(ind1),p1,'*', gd.t(ind2),p2,'d');
legend('True Signal','Observed','FIR','Gamma','Man-Gamma','C1Stim','C2Stim');
gd.hRawTrace = get(gca,'children');
gd.RawTraceId = strvcat('signal','observed','fir','gamma','man-gamma',...
			'c1stim','c2stim');
gd.RawTraceId = flipud(gd.RawTraceId);
gd = HideRaw(gd);

%if(gd.ShowRawSignal)
%else
%plot(gd.t,gd.yobserved, 'o-', ...
%     gd.t,gd.yhatfir, ...
%     gd.t,gd.yhatgam, gd.t,gd.yhatmangam, ...
%     gd.t(ind1),p1,'*', gd.t(ind2),p2,'d');
%legend('Observed','FIR','Gamma','Man-Gamma','C1Stim','C2Stim');
%end

title('Raw Time Courses');
xlabel('time (sec)');
figure(gd.hunivarlab);
guidata(gd.hunivarlab,gd)
return;

%-------------------------------------------------%
% Delete function for Raw window %
function DeleteRaw
ud = get(gcf,'userdata');
gd = guidata(ud.hunivarlab);
gd.ShowRaw = 0;
gd.hRaw = [];
guidata(gd.hunivarlab,gd)
return;

% --------------------------------------------------------------------
function varargout = ShowRaw_mi_Callback(h, eventdata, gd, varargin)
if(gd.ShowRaw) 
  figure(gd.hRaw);
  return;
end
gd.ShowRaw = 1;
gd = PlotRaw(gd);
guidata(gd.hunivarlab,gd)
return;

%-------------------------------------------------%
% Callback for Raw figure toggle switch %
function gd = ToggleRawTrace(traceid)
ud = get(gcf,'userdata');
gd = guidata(ud.hunivarlab);
switch(traceid)
 case 'signal',   gd.ShowRawSignal   = ~gd.ShowRawSignal;
 case 'observed', gd.ShowRawObserved = ~gd.ShowRawObserved;
end
gd = HideRaw(gd);
guidata(gd.hunivarlab,gd)
return;

%-------------------------------------------------%
function gd = HideRaw(gd)
i = strmatch('signal',gd.RawTraceId);
if(gd.ShowRawSignal) set(gd.hRawTrace(i),'visible','on');
else                 set(gd.hRawTrace(i),'visible','off');
end  
i = strmatch('observed',gd.RawTraceId);
if(gd.ShowRawObserved) set(gd.hRawTrace(i),'visible','on');
else                   set(gd.hRawTrace(i),'visible','off');
end  
return;

%-------------------------------------------------%
function gd = PlotACF(gd)

if(gd.ShowACF == 0) return; end
if(~isfield(gd,'hACF') | isempty(gd.hACF) | ~ishandle(gd.hACF)) 
  gd.hACF = figure; 
  set(gcf,'Interruptible','off'); % smooth animations %
  set(gcf,'DoubleBuffer','on');   % smooth animations %
  set(gcf,'BusyAction','cancel'); % dont build up a lot of events %
  set(gcf,'renderer','painters'); % seems to be the best
  set(gcf,'DeleteFcn','univarlab(''DeleteACF'')');
  ud.hShowLagZero = ...
      uicontrol('Style', 'checkbox', ...
		'String', 'Show Lag Zero',...
		'Position', [1 1 100 25], ...
		'Callback', 'univarlab(''ShowLagZero'')',...
		'tooltipstring','Toggle Display of Zeroth Lag',...
		'value',1);
  ud.hunivarlab = gd.hunivarlab;
  set(gd.hACF,'userdata',ud);
end
figure(gd.hACF);

nlag = 10;
if(gd.ACFShowLagZero) lag = 1:nlag;
else                  lag = 2:nlag;
end
z = zeros(size(lag));

if(gd.FixACF)
  if(isempty(gd.FIRFixACFMtx))
    [tmp gd.FIRFixACFMtx] = fast_yacf_kjw(gd.firacf(1:2),gd.RFIR);
  else 
    tmp = gd.FIRFixACFMtx * gd.firacf(1:2);
    tmp = tmp/tmp(1);
  end
  firacf = tmp(2).^[0:gd.ntp-1];

  if(isempty(gd.GamFixACFMtx))
    [tmp gd.GamFixACFMtx] = fast_yacf_kjw(gd.gamacf(1:2),gd.RGam);
  else 
    tmp = gd.GamFixACFMtx * gd.gamacf(1:2);
    tmp = tmp/tmp(1);
  end
  gamacf = tmp(2).^[0:gd.ntp-1];
  %firacf = fast_yacf_kjw(firacf,gd.RFIR);
  %gamacf = fast_yacf_kjw(gamacf,gd.RGam);
else
  firacf = gd.firacf(2).^[0:gd.ntp-1];
  gamacf = gd.gamacf(2).^[0:gd.ntp-1];
end

plot(lag-1,gd.NoiseACF(lag), '+-', ...
     lag-1,gd.firacf(lag), 'o-',...
     lag-1,firacf(lag), 'o-',...
     lag-1,gd.gamacf(lag), 'd-',...
     lag-1,gamacf(lag), 'd-',...
     lag-1,z,'k-..');
title('Autocorrelation Function');
xlabel('lag (time points)');
legend('True ACF','FIR-Raw','FIR-AR1','Gamma-Raw','Gamma-AR1');
guidata(gd.hunivarlab,gd)
figure(gd.hunivarlab);
return;

%-------------------------------------------------%
% Delete function for ACF window %
function DeleteACF
ud = get(gcf,'userdata');
gd = guidata(ud.hunivarlab);
gd.ShowACF = 0;
gd.hACF = [];
guidata(gd.hunivarlab,gd)
return;

% --------------------------------------------------------------------
function varargout = ShowACF_mi_Callback(h, eventdata, gd, varargin)
if(gd.ShowACF)  figure(gd.hACF); return; end
gd.ShowACF = 1;
gd = PlotACF(gd);
guidata(gd.hunivarlab,gd)

% --------------------------------------------------------------------
function varargout = ShowACF_pb_Callback(h, eventdata, gd, varargin)
gd.ShowACF = 1;
gd = PlotACF(gd);
set(gd.ShowACF_pb,'enable','off');
guidata(gd.hunivarlab,gd)
return;

%-------------------------------------------------%
% Callback for ACF figure toggle switch %
function gd = ShowLagZero
ud = get(gcf,'userdata');
gd = guidata(ud.hunivarlab);
gd.ACFShowLagZero = ~gd.ACFShowLagZero;
gd = PlotACF(gd);
guidata(gd.hunivarlab,gd)
return;


%-----------------------------------------------------%
function gd = SetGUIParams(gd)

if(strcmp(gd.ScheduleType,'rper')) v = 1; end
if(strcmp(gd.ScheduleType,'fier')) v = 2; end
if(strcmp(gd.ScheduleType,'blocked')) v = 3; end
set(gd.Schedule_pum,'value',v);

set(gd.GamDeltaSynth_et,'string',sprintf('%g',gd.GamDeltaSynth(1)));
set(gd.GamTauSynth_et,'string',sprintf('%g',gd.GamTauSynth(1)));
set(gd.C1AmpSynth_et,'string',sprintf('%g',gd.C1AmpSynth(1)));
set(gd.C2AmpSynth_et,'string',sprintf('%g',gd.C2AmpSynth(1)));
set(gd.C1AmpMan_et,'string',sprintf('%g',gd.C1AmpMan(1)));
set(gd.C2AmpMan_et,'string',sprintf('%g',gd.C2AmpMan(1)));
set(gd.NoiseStd_et,'string',sprintf('%g',gd.NoiseStd(1)));
set(gd.NoiseAR1_et,'string',sprintf('%g',gd.NoiseAR1(1)));
set(gd.FIRPSDMin_et,'string',sprintf('%g',gd.FIRPSDMin(1)));
set(gd.FIRTotWinTxt,'string',sprintf('%g',gd.FIRTotWin));


return;
% --------------------------------------------------------------------
function gd = NewStimSched(gd)

switch(gd.ScheduleType)
 case {'rper'}
  % randomize sequence and timing
  nc1 = round(gd.ntp/3);
  nc2 = round(gd.ntp/3);
  n0 = gd.ntp - (nc1+nc2);
  r = [zeros(1,n0) ones(1,nc1) 2*ones(1,nc2)]';
  gd.schedule = [gd.t r(randperm(gd.ntp))];
 case {'fier'}
  % randomize sequence, separate by 20 sec
  ntpsep = round(20/gd.TR); 
  nc = round((gd.ntp/ntpsep)/2);
  r = [ones(1,nc) 2*ones(1,nc)];
  r = r(randperm(2*nc));
  r = [r; zeros(ntpsep-1,2*nc)];
  r = reshape1d(r);
  gd.schedule = [gd.t r];
 case {'blocked'}
  % randomize sequence, include null, block length = 20sec
  % May be a problem when using BCW
  ntpblock = round(20/gd.TR);
  nblocks = round(gd.ntp/ntpblock);
  nc1blocks = round(nblocks/4);
  nc2blocks = round(nblocks/4);
  nc0blocks = nblocks-(nc1blocks+nc2blocks);
  seq = [ones(1,nc1blocks) 2*ones(1,nc2blocks)];
  seq = seq(randperm(length(seq)));
  seq = reshape1d([zeros(1,length(seq)); seq])';
  r = reshape1d(repmat(seq,[ntpblock 1]));
  gd.schedule = [gd.t r];
  
end

return;
% --------------------------------------------------------------------
function gd = MakeFirstLevelAnalysis(gd)

set(gd.hunivarlab,'pointer','watch');

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
if(gd.FitMean)
  gd.FIREstFLA.fxlist(3).fx = fast_fxcfg('parseline',...
   [sprintf('effect random drift polynomial %d',0)]);
end
[gd.XFIR gd.FIREstFLA] = fast_fla_desmat(gd.FIREstFLA);
gd.RFIR = eye(gd.ntp) - gd.XFIR*inv(gd.XFIR'*gd.XFIR)*gd.XFIR'; 
gd.FIRDOF = size(gd.XFIR,1)-size(gd.XFIR,2);
gd.FIRDOF = size(gd.XFIR,1)-size(gd.XFIR,2);

  %-------- Est Gamma ------------------------%
gd.GamEstFLA = gd.FIREstFLA;
fxline = sprintf('effect fixed cond1 gamma 1 0 .1 30 0 0 %g %g %d',...
		 gd.GamDeltaEst(1),gd.GamTauEst(1),gd.GamNDeriv);
gd.GamEstFLA.fxlist(1).fx = fast_fxcfg('parseline',fxline);
fxline = sprintf('effect fixed cond2 gamma 2 0 .1 30 0 0 %g %g %d',...
		 gd.GamDeltaEst(1),gd.GamTauEst(1),gd.GamNDeriv);
gd.GamEstFLA.fxlist(2).fx = fast_fxcfg('parseline',fxline);
[gd.XGam gd.GamEstFLA] = fast_fla_desmat(gd.GamEstFLA);
gd.RGam = eye(gd.ntp) - gd.XGam*inv(gd.XGam'*gd.XGam)*gd.XGam'; 
gd.GamDOF = size(gd.XGam,1)-size(gd.XGam,2);

  %-------- Synth Gamma ------------------------%
gd.GamSynthFLA = gd.FIREstFLA;
fxline = sprintf('effect fixed cond1 gamma 1 0 .1 30 0 0 %g %g 0',...
		 gd.GamDeltaSynth(1),gd.GamTauSynth(1));
gd.GamSynthFLA.fxlist(1).fx = fast_fxcfg('parseline',fxline);
fxline = sprintf('effect fixed cond2 gamma 2 0 .1 30 0 0 %g %g 0',...
		 gd.GamDeltaSynth(1),gd.GamTauSynth(1));
gd.GamSynthFLA.fxlist(2).fx = fast_fxcfg('parseline',fxline);
if(gd.FitMean) clear gd.GamSynthFLA.fxlist(3); end
[gd.XGamSynth gd.GamSynthFLA] = fast_fla_desmat(gd.GamSynthFLA);
set(gd.hunivarlab,'pointer','arrow');
return;

% --------------------------------------------------------------------
function gd = SynthNoise(gd)
gd.ynoise0    = randn(gd.ntp,1);
return;

% --------------------------------------------------------------------
function gd = PrepNoiseACF(gd)
gd.NoiseACF = (gd.NoiseAR1(1)).^[0:gd.ntp-1];
gd.NoiseF   = chol(toeplitz(gd.NoiseACF));
return;

% --------------------------------------------------------------------
function gd = ScaleFilterNoise(gd)
gd.ynoise    = gd.NoiseStd(1)*gd.ynoise0;
if(gd.NoiseAR1(1) ~= 0) 
  gd.ynoise  = gd.NoiseF * gd.ynoise;
end
return;

% --------------------------------------------------------------------
function gd = SynthC1Signal(gd)
gd.yc1signal = gd.XGamSynth(:,1) * gd.C1AmpSynth(1);
return;

% --------------------------------------------------------------------
function gd = SynthC2Signal(gd)
gd.yc2signal = gd.XGamSynth(:,2) * gd.C2AmpSynth(1);
return;

% --------------------------------------------------------------------
function gd = SynthSignal(gd)
gd.ysignal   = gd.yc1signal + gd.yc2signal;
return;

% --------------------------------------------------------------------
function gd = SynthObserved(gd)
gd.yobserved = gd.ysignal + gd.ynoise;
gd.cnr = sum(gd.ysignal.^2)/sum(gd.ynoise.^2);
s = sprintf('CNR %6.2f',gd.cnr);
set(gd.CNR_tx,'string',s);
return;

% --------------------------------------------------------------------
function gd = Estimate(gd)
gd.FIRbeta = (inv(gd.XFIR'*gd.XFIR)*gd.XFIR')*gd.yobserved;
gd.yhatfir = gd.XFIR * gd.FIRbeta;
gd.resfir  = gd.yobserved - gd.yhatfir;
gd.resfirvar = sum(gd.resfir.^2)/gd.FIRDOF;
gd.resfirstd = sqrt(gd.resfirvar);
gd.firacf = fast_acorr(gd.resfir);
if(gd.Whiten)
  acf = gd.firacf(2) .^ [0:gd.ntp-1];
  W = inv(chol(toeplitz(acf)));
  y = W*gd.yobserved;
  X = W*gd.XFIR;
  gd.FIRbeta = (inv(X'*X)*X')*y;
  gd.yhatfir = gd.XFIR * gd.FIRbeta;
  gd.resfir  = W*(gd.yobserved - gd.yhatfir);
  gd.resfirvar = sum(gd.resfir.^2)/gd.FIRDOF;
  gd.resfirstd = sqrt(gd.resfirvar);
  gd.firacf = fast_acorr(gd.resfir);
end

gd.Gambeta = (inv(gd.XGam'*gd.XGam)*gd.XGam')*gd.yobserved;
gd.yhatgam = gd.XGam * gd.Gambeta;
gd.resgam  = gd.yobserved - gd.yhatgam;
gd.resgamvar = sum(gd.resgam.^2)/gd.GamDOF;
gd.resgamstd = sqrt(gd.resgamvar);
gd.gamacf = fast_acorr(gd.resgam);
if(gd.Whiten)
  acf = gd.gamacf(2) .^ [0:gd.ntp-1];
  W = inv(chol(toeplitz(acf)));
  y = W*gd.yobserved;
  X = W*gd.XGam;
  gd.Gambeta = (inv(X'*X)*X')*y;
  gd.yhatgam = gd.XGam * gd.Gambeta;
  gd.resgam  = W*(gd.yobserved - gd.yhatgam);
  gd.resgamvar = sum(gd.resgam.^2)/gd.GamDOF;
  gd.resgamstd = sqrt(gd.resgamvar);
  gd.gamacf = fast_acorr(gd.resgam);
end

set(gd.FIRRStdTxt,'string',sprintf('%3.2f',gd.resfirstd));
set(gd.GamRStdTxt,'string',sprintf('%3.2f',gd.resgamstd));
set(gd.GamRStdTxt,'string',sprintf('%3.2f',gd.resgamstd));
set(gd.C1AmpEstTxt,'string',sprintf('%3.2f',gd.Gambeta(1)));
set(gd.C2AmpEstTxt,'string',sprintf('%3.2f',gd.Gambeta(2+gd.GamNDeriv)));

gd.FIRFixACFMtx = []; % Will need new matrix
gd.GamFixACFMtx = []; % Will need new matrix
gd = EstimateManual(gd);

return;

%----------------------------------------------------------%
function gd = EstimateManual(gd);

b = [gd.C1AmpMan(1) gd.C2AmpMan(1); zeros(gd.GamNDeriv,2)];
b = reshape1d(b);
if(gd.FitMean) b = [b; 0]; end
gd.yhatmangam = gd.XGam * b;
gd.resmangam  = gd.yobserved - gd.yhatmangam;
gd.resmangamvar = sum(gd.resmangam.^2)/gd.GamDOF;
gd.resmangamstd = sqrt(gd.resmangamvar);
set(gd.ManGamRStdTxt,'string',sprintf('%3.2f',gd.resmangamstd));
return;


























