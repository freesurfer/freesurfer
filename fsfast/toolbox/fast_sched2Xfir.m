function X = fast_sched2Xfir(tPres,ntrs,TR,IRFPSD,BCW,tDelay,PerEvW)
%
% X = fast_sched2Xfir(tPres,ntrs,TR,IRFPSD,BCW,tDelay,PerEvW)
%
% Creates a design matrix (aka stimulus convolution matrix) modeling
% the hemodynamic response as an FIR with adjustable tap weights. The
% matrix is for a single event type whose schedule is passed by tPres.
%
% tPres - list of presentation times (ie, the schedule) in seconds 
% for one event type. Time=tDelay is defined as the time the first 
% stored image was collected (ie, ignore discarded acquisitions). An 
% event before t=0 will have an effect on the matrix if its Post 
% Stimulus Window encompasses t=0; otherwise it is ignored. Events 
% found after the termination of data collection are ignored (a warning 
% is printed). If tPres is empty, a matrix of zeros is returned. 
%
% ntrs - number of functional volumes collected (do not include 
% discarded acquisitions or prescan time).
%
% TR - TR of the exmperiment (ie, time between acquisitions of 
% functional volumes) in seconds.
%
% IRFPSD - three component vector indicating the Post Stimulus
% Delay window of the Impulse Response Function. The three
% components are: PSDMin dPSD PSDMax. dPSD has also been called
% the TER. Note: this models the Impulse Response only. The
% Event Response is the Impulse Response convolved with 
% a Box Car (see BCW)
%
% BCW - Box Car Width (in sec). The Event Response is the Impulse
% Response convolved with a Box Car of width BCW. The Event
% Response PSD Window is that of the IRF except that 
% ERFPSDMax = IRFPSDMax + BCW.
%
% tDelay - set the definition of t=0. This can be handy for
% taking the slice acquisition delay into account. Ie, creating
% a different X matrix for each slice. In this case, set tDelay
% equal to the time between when the first slice is acquired 
% and the nth slice.
%
% PerEvW - presentation weighting. Length should be NPresentations. The
% matrix entry for each presentation is given the value W(n) instead
% of 1. Ignored if W=[].
%
% X will have size: Ntp by (ERFPSDMax-IRFPSDMin)/dIRFPSD
% 
% $Id: fast_sched2Xfir.m,v 1.3 2003/03/17 06:29:23 greve Exp $ 

X = [];

if(nargin < 4 & nargin > 7)
  msg = 'X = fast_sched2Xfir(tPres,ntrs,TR,IRFPSD,<BCW,tDelay,PerEvW>)';
  fprintf('%s\n',msg);
  return;
end

if(exist('BCW') ~= 1) BCW = 0; end
if(exist('tDelay') ~= 1) tDelay = 0; end

% Compute Event Response Function Post Stimulus Window %
PSDMin = IRFPSD(1);
dPSD   = IRFPSD(2);
PSDMax = IRFPSD(3)+BCW;

% Check that the TR is an integer mult of the dPSD %
if(rem(TR,dPSD) ~= 0)
  msg = sprintf('ERROR: TR (%g) must be a multiple of dPSD (%g)',...
		TR,dPSD);
  fprintf('%s\n',msg);
  return;
end

% Check that the PSDMin is an integer mult of the dPSD %
if(rem(PSDMin,dPSD) ~= 0)
  msg = sprintf('ERROR: PSDMin (%g) must be a multiple of dPSD (%g)',...
		PSDMin,dPSD);
  fprintf('%s\n',msg);
  return;
end

% Check that the BCW is an integer mult of the dPSD %
if(rem(BCW,dPSD) ~= 0)
  msg = sprintf('ERROR: BCW (%g) must be a multiple of dPSD (%g)',...
		BCW,dPSD);
  fprintf('%s\n',msg);
  return;
end

% Check that the PSDMax is an integer mult of the dPSD %
if(rem(PSDMax,dPSD) ~= 0)
  msg = sprintf('ERROR: PSDMax (%g) must be a multiple of dPSD (%g)',...
		PSDMax,dPSD);
  fprintf('%s\n',msg);
  return;
end

TimeWindow = PSDMax - PSDMin;

% Compute time of last acq
tmax = TR*(ntrs - 1);
% Compute the resampling rate
Rss = TR/dPSD;
% Compute the number of Estimates
Nh = round(TimeWindow/dPSD);
% Compute the Post Stimulus Delay at the start of the last 
% point in the window
TPostStim = PSDMax - dPSD;

% Number of presentations
Npres = length(tPres);
if(Npres == 0) 
  X = zeros(ntrs,Nh);
  return; 
end
tPres = reshape(tPres,[Npres 1]); % Just to make sure %

if(exist('PerEvW') ~= 1) PerEvW = []; end
if(isempty(PerEvW)) PerEvW = ones(Npres,1); end
if(length(PerEvW) ~= Npres)
  fprintf('ERROR: length(PerEvW) (%d) does not equal Npres (%d)\n',...
	  length(PerEvW),Npres);
  return;
end
PerEvW = reshape(PerEvW,[Npres 1]); % Just to make sure %

% Subtract tDelay from tPres
tPres = tPres - tDelay;

X = zeros(Rss*ntrs,Nh);
h = 1;
for d = PSDMin:dPSD:TPostStim,
   td = tPres+d;
   iok = find(td >= 0 & td <= tmax);
   td = td(iok);
   iPres = round(td/dPSD)+1;
   X(iPres,h) = PerEvW(iok);
   h = h + 1;
end

% SubSample %
if(Rss ~= 1)
  X = X(1:Rss:Rss*ntrs,:);
end

if(BCW ~= 0)
  Nb = BCW/dPSD;
  Na = (IRFPSD(3)-PSDMin)/dPSD;
  P = [zeros(Nb-1,Na); eye(Na); ];
  b = [ones(1,Nb) zeros(1,Na-1)];
  c = zeros(Na+Nb,1);
  c(1) = 1;
  B = toeplitz(c,b);
  X = X*B*P;
end

return;
