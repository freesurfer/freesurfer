function X = fast_sched2Xfir(tPres,ntrs,TR,psdwin,tDelay,PerEvW)
%
% X = fast_sched2Xfir(tPres,ntrs,TR,psdwin,tDelay,PerEvW)
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
% is printed). If tPres is empty, a matrix of the correct size
% filled with zeros is returned. 
%
% ntrs - number of functional volumes collected (do not include 
% discarded acquisitions or prescan time).
%
% TR - TR of the exmperiment (ie, time between acquisitions of 
% functional volumes) in seconds.
%
% psdwin - three/four component vector indicating the Post Stimulus
% Delay window. The components are: psdmin dpsd psdmax bcw. dpsd 
% has also been called the TER. bcw is the Box Car Width. If bcw
% is not present, it is assumed to be 0. See fast_psdwin for more 
% details
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
% X will have size: Ntp by (psdmax+bcw-psdmin-dpsd)/dpsd
% 
% See also: fast_psdwin
%
%


%
% fast_sched2Xfir.m
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

X = [];

if(nargin < 4 & nargin > 6)
  msg = 'X = fast_sched2Xfir(tPres,ntrs,TR,psdwin,<tDelay,PerEvW>)';
  fprintf('%s\n',msg);
  return;
end

if(exist('tDelay') ~= 1) tDelay = []; end
if(isempty(tDelay)) tDelay = 0; end

% Compute number of columns of X
Nh = fast_psdwin(psdwin,'npsdwin');
if(isempty(Nh)) return; end

psdmin = psdwin(1);
dpsd   = psdwin(2);
psdmax = psdwin(3);
if(length(psdwin) == 3) bcw = 0;
else bcw = psdwin(4);
end

TimeWindow = psdmax - psdmin;

% Compute time of last acq
tmax = TR*(ntrs - 1);
% Compute the resampling rate
Rss = round(TR/dpsd);
% Compute the Post Stimulus Delay at each point in the window
psdlist = fast_psdwin(psdwin,'erftaxis');

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
for d = psdlist'
   td = tPres+d;
   iok = find(td >= 0 & td <= tmax);
   td = td(iok);
   iPres = round(td/dpsd)+1;
   X(iPres,h) = PerEvW(iok);
   h = h + 1;
end

% SubSample %
if(Rss ~= 1)
  X = X(1:Rss:Rss*ntrs,:);
end

if(bcw ~= 0)
  B = fast_boxcarmat(psdwin);
  X = X*B;
end

return;
