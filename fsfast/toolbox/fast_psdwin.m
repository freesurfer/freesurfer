function rt = fast_psdwin(psdwin,DoWhat)
% rt = fast_psdwin(psdwin,DoWhat)
%
% Manages the Post-Stimulus Delay Window
%
% psdwin = [psdmin dpsd psdmax <bcw>]
%
% DoWhat = 'check', 'npsdwin', 'erfpsdwin'
%  If DoWhat is null or does not exist, then 'check'
% 
% psdmin is the time (sec) at which the window starts
%   This is similar to the "prestim" time, except that
%   psdmin = -prestim.
%
% dpsd is the temporal resolution inside the window.
%   This is the same as TER. In the past TER was forced
%   to apply to all event types.
%
% psdmax is the maximum post-stimulus delay of the Impulse 
%   Response Window. This is related to "TimeWindow" :
%   TimeWindow = psdmax - psdmin - dpsd
%   The '-dpsd' is needed because it is assumed that each
%   point in the window spans dpsd seconds, so the last
%   point ends at psdmax.
%
% bcw - length of box car (in seconds) to convolve the
%   impulse response. This will lengthen the length of the
%   window. If not present, then bcw = 0.
%
% DoWhat can be:
%   not present - same as 'check'
%   empty - same as 'check'
%   'check' - checks window, returns 1 if ok, 0 otherwise
%   'erfpsdmax' - max PSD of the Event Response Function,
%      computed as erfpsdmax = psdmax + bcw, ie the IRF max
%      after convolving with a boxcar of length bcw.
%   'npsdwin' - the number of dpsds in the window, computed
%   DoWhat is not case sensitive.
%  
% The psd window is checked regardless of what DoWhat is. If
% there is an error and DoWhat is not 'check', then an empty
% matrix is returned.
%  
% $Id: fast_psdwin.m,v 1.1 2003/03/18 06:21:01 greve Exp $

rt = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('rt = fast_psdwin(psdwin,DoWhat)\n');
  return;
end
if(exist('DoWhat') ~= 1) DoWhat = ''; end
if(isempty(DoWhat)) DoWhat = 'check'; end
DoWhat = lower(DoWhat);

lenpsdwin = length(psdwin);
if(lenpsdwin ~= 3 & lenpsdwin ~= 4)
  fprintf('ERROR: psdwin has %d components, must be 3 or 4\n',lenpsdwin);
  return;
end

psdmin = psdwin(1);
dpsd   = psdwin(2);
psdmax = psdwin(3);
if(lenpsdwin == 3) bcw = 0;
else               bcw = psdwin(4);
end

% Check PSD no matter what %
errmsg = '';
if(dpsd <= 0)
  errmsg = strvcat(errmsg,sprintf('ERROR: dpsd = %g, must be > 0\n',dpsd));
end
if(psdmax <= psdmin)
  errmsg = strvcat(errmsg,...
    sprintf('ERROR: psdmax = %g <= dpsdmin = %g\n',psdmax,psdmin));
end
if(bcw < 0)
  errmsg = strvcat(errmsg,sprintf('ERROR: bcw = %g, must be >= 0\n',bcw));
end
if(rem(bcw,dpsd) ~= 0)
  errmsg = strvcat(errmsg,...
     sprintf('ERROR: bcw=%g not int mult of dpsd=%g\n',bcw,dpsd));
end
if(rem(psdmin,dpsd) ~= 0)
  errmsg = strvcat(errmsg,...
     sprintf('ERROR: psdmin=%g not int mult of dpsd=%g\n',psdmin,dpsd));
end
if(rem(psdmax,dpsd) ~= 0)
  errmsg = strvcat(errmsg,...
     sprintf('ERROR: psdmax=%g not int mult of dpsd=%g\n',psdmax,dpsd));
end

if(~isempty(errmsg))
  for n = 1:size(errmsg,1)
    fprintf('%s\n',deblank(errmsg(n,:)));
  end
  if(strcmp(DoWhat,'check')) rt = 0; end
  return;
else
  ok = 1;
end

erfpsdmax = psdmax + bcw;
npsdwin = (erfpsdmax-psdmin-dpsd)/dpsd;

switch(DoWhat)
  case 'check'
   rt = ok;
  case 'npsdwin'
   rt = npsdwin;
  case 'erfpsdmax'
   rt = erfpsdmax;
 otherwise
  fprintf('ERROR: DoWhat=%s, unrecognized\n',DoWhat);
  fprintf(' Legal options are: check, erfpsdmax, and npsdwin\n');
  return;
end

return;
