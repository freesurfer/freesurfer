function X = fast_sched2Xfir(tPres,ntrs,TR,TER,TPreStim,TimeWindow,W)
%
% X = fast_sched2Xfir(tPres,ntrs,TR,TER,TPreStim,TimeWindow,W)
%
% Creates a design matrix (aka stimulus convolution matrix) modeling
% the hemodynamic response as an FIR with adjustable tap weights. The
% matrix is for a single event type whose schedule is passed by tPres.
% The size of the matrix will be ntrs by Nh, Nh = TimeWindow/TER.
%
% tPres - list of presentation times (ie, the schedule) in seconds 
% for one event type. Time=0 is defined as the time the first stored 
% image was collected (ie, ignore discarded acquisitions). An event
% before t=0 will have an effect on the matrix if its Post Stimulus 
% Window encompasses t=0; otherwise it is ignored. Events found
% after the termination of data collection are ignored (a warning 
% is printed). If tPres is empty, a matrix of zeros is returned. 
%
% ntrs - number of functional volumes collected (do not include 
% discarded acquisitions or prescan time).
%
% TR - TR of the exmperiment (ie, time between acquisitions of 
% functional volumes) in seconds.
%
% TER - temporal estimation resolution (in seconds). This is the
% resolution at which the hemodynamic response is reconstructed.
% The TER must be an integer divisor of the TR. Note: arbitrary
% schedules cannot be used to perform sub TR estimation. The 
% schudule must have been designed for that purpose.
%
% TPreStim - prestimulus window (in seconds). This creates a matrix
% which will allow averaging to begin before the stimulus onset.
% Must be a integer multiple of the TER.
%
% TimeWindow - total time window (in seconds) in which the hemodynamic
% response will be modeled, including the TPreStim. Must be a integer 
% multiple of the TER. The Post Stimulus Window = TimeWindow - TPreStim.
%
% W - presentation weighting. Length should be NPresentations. The
% matrix entry for each presentation is given the value W(n) instead
% of 1. Ignored if W=[].
%
% $Id: fast_sched2Xfir.m,v 1.1 2003/03/04 20:47:38 greve Exp $ 

X = [];

if(nargin ~= 7)
  msg = 'USAGE: X = fast_sched2Xfir(tPres,ntrs,TR,TER,TPreStim,TimeWindow,W)';
  fprintf('%s\n',msg);
  return;
end

% Check that the TR is an integer mult of the TER %
if(rem(TR,TER) ~= 0)
  msg = sprintf('ERROR: TR (%g) must be a multiple of TER (%g)',...
		TR,TER);
  fprintf('%s\n',msg);
  return;
end

% Check that the TW is an interger mult of the TER %
if(rem(TimeWindow,TER) ~= 0)
  msg = sprintf('ERROR: TimeWindow (%g) must be a multiple of TER (%g)',...
		TimeWindow,TER);
  fprintf('%s\n',msg);
  return;
end

% Check that the TPreStim is an interger mult of the TER %
if(rem(TPreStim,TER) ~= 0)
  msg = sprintf('ERROR: TPreStim (%g) must be a multiple of TER (%g)',...
		TPreStim,TER);
  fprintf('%s\n',msg);
  return;
end

% Check that no presentation exceeds the maximum time %
tmax = TR*(ntrs - 1);
ind = find(tPres > tmax);
if(~isempty(ind))
  fprintf('WARNING: presentation time %g exceeds max %g ... ignoring\n',...
          max(tPres),tmax);
  ind = find(tPres <= tmax);
  tPres = tPres(ind);
end

Npres = length(tPres);
if(isempty(W)) W = ones(Npres,1); end
if(length(W) ~= Npres)
  fprintf('ERROR: length(W) (%d) does not equal Npres (%d)\n',length(W),Npres);
  return;
end

% Compute the resampling rate
Rss = TR/TER;

% Compute the number of Estimates
Nh = round(TimeWindow/TER);

% Compute the Post Stimulus Window
TPostStim = TimeWindow - TPreStim - TER;

X = zeros(Rss*ntrs,Nh);

h = 1;
for d = -TPreStim:TER:TPostStim,
   td = tPres+d;
   iok = find(td >= 0 & td <= tmax);
   td = td(iok);
   iPres = round(td/TER)+1;
   X(iPres,h) = W(iok);
   h = h + 1;
end

% SubSample %
if(Rss ~= 1)
  X = X(1:Rss:Rss*ntrs,:);
end


return;
