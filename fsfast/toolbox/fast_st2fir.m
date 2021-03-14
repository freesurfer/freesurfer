function Xfir = fast_st2fir(st,ntp,TR,psdwin,usew)
% Xfir = fast_st2fir(st,ntp,TR,psdwin,<usew>)
%
% Computes the FIR design matrix for the given schedule of one
% event type.
%
% st = [tonset duration <weight>]
% ntp = number of time points
% TR = time between time points
% psdwin = [psdmin psdmax dpsd];
% usew:  0, [], not spec = dont use weight, 1 = use weight
%
% Notes:
%  1. Number of rows in st is the number of presentations.
%  2. Presenations may not overlap --
%     use different event types in this case.
%  3. Set your psdmax to be long enough for IRF and max duration
%  4. If st does not have weights, weights=1
%  5. Does not force dpsd to be an integer divisor of TR,
%     but it is a good idea.
%  6. If two stimuli fall within the sam TR bin, the Xfir
%     matrix will have a 2 instead of a 1.


%
% fast_st2fir.m
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

Xfir = [];
a = [];
Xfirss = [];

if(nargin < 4 | nargin > 6)
  fprintf('Xfir = fast_st2fir(st,ntp,TR,psdwin,<usew>)\n');
  return;
end
if(~exist('usew','var')) usew = []; end
if(isempty(usew)) usew = 0; end

psdmin  = psdwin(1);  % start of PSD window
psdmax  = psdwin(2);  % end of PSD window
dpsd    = psdwin(3);  % increment of PSD window
npsdwin = round((psdmax-psdmin)/dpsd);

% Empty st means that the condition is not present. This
% can only happen when the user has specified
% flac.AllowMissingCond. Return matrix of all 0s
if(isempty(st))
  Xfir = zeros(ntp,npsdwin);
  return;
end

nstim = size(st,1);

% Presentation duration - use dpsd if not specified
if(size(st,2) < 2)
  stimdur = dpsd*ones(nstim,1);
else
  stimdur = st(:,2); 
end

% Presentation weights - if not specified, set them to 1
if(size(st,2) < 3 | ~usew) 
  stimweight = ones(nstim,1);
else
  stimweight = st(:,3);
end

% The following two pieces of code prevent the case where stimuli
% of the same type are presented in a overlapping manner. It is,
% of course, ok to have the presentations of different stimulus
% types overlap (or even be simultaneous).
st = sortrows(st); % sort in order of onset time
d = diff(st(:,1));
ind = find(d == 0);
if(~isempty(ind))
  fprintf('ERROR: fast_st2fir: two or more presentations are simultaneous\n');
  return;
end

% Presentation onset times
stimonset = st(:,1); 

% SubSampling
ssr = round(TR/dpsd); % Rate

Xfirss = zeros(ntp*ssr,npsdwin);
for nthstim = 1:nstim
  tonset   = st(nthstim,1); % stimulus onset time
  wstim    = stimweight(nthstim);   % stimulus weight
  stimdur  = st(nthstim,2); % stimulus duration (sec)
  nstimdur = round(stimdur/dpsd); % stimulus duration (indices)
  for nthdur = 1:nstimdur
    % Loop through each delay (column index in X)
    for nthd = 1:npsdwin
      % time of the nth delay of the nth duration
      td = tonset + psdmin + (nthd-1)*dpsd + dpsd*nthdur;
      % Row index in X
      ind = round(td/dpsd);
      if(ind < 1 | ind > ssr*ntp) continue; end
      Xfirss(ind,nthd) = wstim;
    end
  end
end

% Subsample 
Xfir = Xfirss(1:ssr:end,:);

return;
