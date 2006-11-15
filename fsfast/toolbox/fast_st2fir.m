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
%
% $Id: fast_st2fir.m,v 1.7 2006/11/15 23:08:00 greve Exp $

Xfir = [];

if(nargin < 4 | nargin > 5)
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

% If weights are not specified, set them to 1
if(size(st,2) < 3 | ~usew) st(:,3) = 1; end

npres   = size(st,1); % number of presentations

% The following two pieces of code prevent the case where stimuli
% of the same type are presented in a overlapping manner. It is,
% of course, ok to have the presentations of different stimulus
% types overlap (or even be simultaneous).

% Make sure that presentations are not simultaneous
st = sortrows(st); % sort in order of onset time
d = diff(st(:,1));
ind = find(d == 0);
if(~isempty(ind))
  fprintf('ERROR: fast_st2fir: two or more presentations are simultaneous\n');
  return;
end

% Check whether the offset of one stimulus is past the onset of the
% next. Allows them to overlay by less than dpsd/2.
tonset  = st(:,1);
toffset = st(:,1) + st(:,2);
ind = find(toffset(1:end-1) > (tonset(2:end) + dpsd/2));
if(~isempty(ind))
  fprintf('ERROR: fast_st2fir: two or more presentations overlap\n');
  return;
end

% Alloc and set to 0
Xfir = zeros(ntp,npsdwin);

% Go through each presentation
for nthpres = 1:npres
  tonset0  = st(nthpres,1);
  duration = st(nthpres,2);
  weight   = st(nthpres,3);
  if(duration ==0) nduration = 1;
  else             nduration = duration/dpsd;
  end

  % Go through each increment in the duration
  for nthduration = 1:nduration
    tonset = tonset0 + (nthduration-1)*dpsd;
  
    % Rows in the design matrix (0-based)
    r1 = round((tonset+psdmin)/TR);
    r2 = round((tonset+psdmax)/TR)-1;
    r = r1:r2;
    
    % Columns in the design matrix (0-based)
    c = round((TR*r-tonset)/dpsd);
    
    % Convert to 1-based
    r = r + 1;
    c = c + 1;
    
    % Only keep the ones that are in bounds
    indok = find(r > 0 & r <= ntp & c > 0 & c <= npsdwin);
    if(isempty(indok)) continue; end
    r = r(indok);
    c = c(indok);

    % Compute the indicies in the design matrix
    ind = sub2ind(size(Xfir),r,c);
    % Set the components in the design matrix to the weight
    Xfir(ind) = Xfir(ind) + weight;
  end
  
end

return
