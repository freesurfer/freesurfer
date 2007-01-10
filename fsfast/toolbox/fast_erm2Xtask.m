function [X, Nc] = fast_par2Xtask(par,ntrs,TR,ERM,W)
% [X, Nc] = fast_par2Xtask(par,ntrs,TR,ERM,W)
%
% ERM - vector of event-response model structures. The
% structure has (at least) the following fields:
%   name - string with the name of the model
%   params - vector with parameter list
%
% If the length of ERM is 1, then the same ERM is applied
% to all event types.
%
% Possible ERMs and their parameters:
%  fir         - tprestim, ter, timewindow
%  gamma       - delay, dispersion, boxcarwidth
%  gamma+deriv - delay, dispersion, boxcarwidth
%
% See also:
%   fast_par2nconds, fast_sched2Xfir, fast_sched2Xgamma, 
%   fast_sched2Xgammaderiv


%
% fast_erm2Xtask.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:30 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

X = [];
Nc = [];

if(nargin ~= 5)
  fprintf('[X, Nc] = fast_par2Xtask(par,ntrs,TR,ERM,W)\n');
  return;
end

% Get the number of conditions %
[Nc CondList Holes] = fast_par2nconds(par);
if(Holes)
  fprintf('ERROR: holes in the paradigm\n');
  return;
end

if(length(ERM) ~= Nc & length(ERM) ~= 1)
  fprintf('ERROR: Number of ERMs does not equal number of conditions\n');
  return;
end

X0 = [];
for c = 1: Nc
  indc = find(par(:,2)==c);

  tPres = par(indc,1);

  if(~isempty(W)) Wc = W(indc);
  else Wc = [];
  end

  if(length(ERM) ~= 1) ERMc = ERM(c);
  else                 ERMc = ERM;
  end

  switch(lower(ERMc.name))

     case 'fir',
       TPreStim   = ERMc.params(1);
       TER        = ERMc.params(2);
       TimeWindow = ERMc.params(3);
       Xc = fast_sched2Xfir(tPres,ntrs,TR,TER,TPreStim,TimeWindow,Wc);
       if(isempty(Xc)) return; end

     case 'gamma',
       Delay       = ERMc.params(1);
       Dispersion  = ERMc.params(2);
       BoxCarWidth = ERMc.params(3);
       Xc = fast_sched2Xgamma(tPres,ntrs,TR,Delay,Dispersion,BoxCarWidth,Wc);
       if(isempty(Xc)) return; end

     case 'gamma+deriv',
       Delay       = ERMc.params(1);
       Dispersion  = ERMc.params(2);
       BoxCarWidth = ERMc.params(3);
       Xc = fast_sched2Xgammaderiv(tPres,ntrs,TR,Delay,Dispersion,BoxCarWidth,Wc);
       if(isempty(Xc)) return; end

     otherwise,
       fprintf('ERROR: Event Response Model %s unrecognized\n',ERMc.name);
       return;

  end

  X0 = [X0 Xc];
end

X = X0;

return;
