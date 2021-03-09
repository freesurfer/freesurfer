function varargout = rawplot(varargin)
% rawplot Application M-file for rawplot.fig
% fig = rawplot launch rawplot GUI.
%    rawplot('callback_name', ...) invoke the named callback.


%
% rawplot.m
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
  % Launch new GUI
  fig = openfig(mfilename,'new');
  set(fig,'name','FSFAST Raw Plot');

  % Use system color scheme for figure:
  set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

  % Generate a structure of handles to pass to callbacks, and store it. 
  gd = guihandles(fig); % only run once
  gd.hfig = fig;
  gd = InitGUIParams(gd);
  gd = CheckVolumes(gd);
  guidata(fig, gd);

  if(nargout > 0) varargout{1} = fig; end

elseif(ischar(varargin{1}))
  % INVOKE NAMED SUBFUNCTION OR CALLBACK
  try
    if(nargout) [varargout{1:nargout}] = feval(varargin{:}); 
    else feval(varargin{:});
    end
  catch disp(lasterr);
  end
end

return;
%------------------------------------------------------%
%------------------------------------------------------%
%------------------------------------------------------%

%------------------------------------------------------%
function gd = InitGUIParams(gd)
gd.volstemlist = '';
%gd.volstemlist = 'mgh-tdr-e2a-1/bold/004/f';
%gd.volstemlist = strvcat(gd.volstemlist,'mgh-tdr-e2b-1/bold/004/f');
gd.polyfitorder = 1;
gd.CurPoint = [1 1 1]; % row,col,slice 1-based

% Get these from volume
gd.volsize = []; % rows columsn slices
gd.mristruct = []; 
gd.fft = 0;
gd.SumTC = 0;
gd.PolyFit = 1;

return;

%------------------------------------------------------%
function AddVolumes(stemlist,hfig)
nstems = size(stemlist,1);
if(nstems == 0) return; end
gd = guidata(hfig);

for n = 1:nstems
  stem = deblank(stemlist(n,:));
  gd.volstemlist = strvcat(gd.volstemlist,stem);
end
gd = CheckVolumes(gd);
gd = LoadTimeCourses(gd);
gd = PlotCurPoint(gd);
guidata(hfig,gd);
return;

%------------------------------------------------------%
function gd = CheckVolumes(gd)

nvols = size(gd.volstemlist,1);
if(nvols == 0) return; end

for n = 1:nvols
  stem = deblank(gd.volstemlist(n,:));
  [nslices nrows ncols nframes] = fmri_bvoldim(stem);
  if(nslices==0)
    fprintf('ERROR loading %s\n',stem);
    return;
  end
end

gd.volsize = [nrows ncols nslices];
gd.mristruct = fast_ldbhdr(stem);

set(gd.TR_st,'string',sprintf('TR = %g',gd.mristruct.tr));

return

%------------------------------------------------------%
function gd = LoadTimeCourses(gd)

nvols = size(gd.volstemlist,1);
if(nvols == 0) return; end

r = gd.CurPoint(1);
c = gd.CurPoint(2);
s = gd.CurPoint(3);

gd.tc = [];
for n = 1:nvols
  stem = deblank(gd.volstemlist(n,:));
  tc = fast_ldbvoxel(stem,c,r,s,1);
  if(isempty(tc))
    fprintf('ERROR loading %s\n',stem);
    return;
  end
  gd.tc = [gd.tc tc];
end

return

%------------------------------------------------------%
function SetCurPoint(r,c,s,hfig)
gd = guidata(hfig);
if(r < 1 | r > gd.volsize(1)) return; end
if(c < 1 | c > gd.volsize(2)) return; end
if(s < 1 | s > gd.volsize(3)) return; end

gd.CurPoint = [r c s];
gd = LoadTimeCourses(gd);
gd = PlotCurPoint(gd);
set(gd.CurRCS_st,'string',sprintf('CurRCS = %d %d %d',c,r,s));
guidata(hfig,gd);
return;

%------------------------------------------------------%
function gd = PlotCurPoint(gd)

if(isempty(gd.tc)) return; end
figure(gd.hfig);
nframes = size(gd.tc,1);
tc = gd.tc;
if(gd.SumTC) tc = sum(tc,2); end
if(gd.PolyFit >= 0)
  Xdt = fast_polytrendmtx(1,nframes,1,gd.PolyFit);
  tc = (eye(nframes)-Xdt*inv(Xdt'*Xdt)*Xdt')*tc;
end
if(gd.fft) 
  afft = abs(fft(tc,[],1));
  afft(1,:) = 0;
  afft = afft(1:round(nframes/2),:);
  fftaxis = fast_fftaxis(nframes,gd.mristruct.tr)';
  fftaxis = fftaxis(1:round(nframes/2),:);
  gd.hplot = plot(fftaxis,afft);
  xlabel('frequency (Hz)');
else
  gd.t = gd.mristruct.tr * [0:nframes-1];
  gd.hplot = plot(gd.t,tc);
  xlabel('time (sec)');
end

return
% --------------------------------------------------------------------
function varargout = FFT_cb_Callback(h, eventdata, gd, varargin)
gd.fft = ~gd.fft;
gd = PlotCurPoint(gd);
guidata(gd.hfig,gd);
return;

% --------------------------------------------------------------------
function varargout = Detrend_pm_Callback(h, eventdata, gd, varargin)
v = get(h,'value');
%fprintf('Detrend value = %d\n',v);
gd.PolyFit = v-2;
gd = PlotCurPoint(gd);
guidata(gd.hfig,gd);
return;



% --------------------------------------------------------------------
function varargout = SumTC_cb_Callback(h, eventdata, gd, varargin)
gd.SumTC = ~gd.SumTC ;
gd = PlotCurPoint(gd);
guidata(gd.hfig,gd);
return;
