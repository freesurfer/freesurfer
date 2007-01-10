function [Xfir, Nc] = fast_par2Xfir(par,ntrs,TR,TER,TPreStim,TimeWindow,W)
% [Xfir Nc] = fast_par2Xfir(par,ntrs,TR,TER,TPreStim,TimeWindow,W)


%
% fast_par2Xfir.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

Xfir = [];
Nc = [];

if(nargin ~= 7)
fprintf('[Xfir Nc] = fast_par2Xfir(par,ntrs,TR,TER,TPreStim,TimeWindow,W)\n');
return;
end

% Get the number of conditions %
Nc = fast_par2nconds(par);

Xfir = [];
for c = 1: Nc
  indc = find(par(:,2)==c);
  tPres = par(indc,1);
  if(~isempty(W)) Wc = W(indc);
  else Wc = [];
  end
  Xfirc = fast_sched2Xfir(tPres,ntrs,TR,TER,TPreStim,TimeWindow,Wc);
  Xfir = [Xfir Xfirc];
end

return;
