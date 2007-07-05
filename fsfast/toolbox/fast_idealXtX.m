function [XtX C vrf] = fast_idealXtX(EventProb,Nh)
% [XtX C vrf] = fast_idealXtX(EventProb,Nh)
%
% Computes "ideal" FIR XtX, where "ideal" is the asymptotic XtX for an
% infinite length run. Assumes each event type lasts 1
% TER. Multiply by Ntp to get the ideal for a certain number of
% time points.
%
% EventProb - list of probabilities of a non-null event types. The
%   sum of EventProb < 1.
% Nh - number of FIR taps per event type
%
% C = ones(1,Nstim*Nh)/(Nstim*Nh);
% vrf = 1/(C*inv(XtX)*C');
%
% There is a simulation at the end of this file (after return) to
% check the accuracy.
%
% $Id: fast_idealXtX.m,v 1.1 2007/07/05 16:42:08 greve Exp $

%
% fast_idealXtX
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/07/05 16:42:08 $
%    $Revision: 1.1 $
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

XtX = [];
if(nargin ~= 2)
  fprintf('[XtX C vrf] = fast_idealXtX(EventProb,Nh)\n');
  return;
end

ptot = sum(EventProb);
if(ptot >= 1)
  fprintf('ERROR: fast_idealXtX: ptot = %g >= 1\n',ptot);
  return;
end

NeventTypes = length(EventProb);

nXtX = NeventTypes*Nh;
XtX = zeros(nXtX,nXtX);

i1 = 1;
for ithEventType = 1:NeventTypes
  i2 = i1+Nh-1;
  j1 = 1;
  Pi = EventProb(ithEventType);
  for jthEventType = 1:NeventTypes
    j2 = j1+Nh-1;
    Pj = EventProb(jthEventType);
    XtXij = Pi*Pj*ones(Nh,Nh);
    if(ithEventType == jthEventType)
      % Replace diag with Pi=Pj
      XtXij = XtXij - diag(diag(XtXij)) + Pi*eye(Nh,Nh);
    else
      % Replace diag with 0, ie, can't have two events
      % at the same time.
      XtXij = XtXij - diag(diag(XtXij));
    end
    XtX(i1:i2,j1:j2) = XtXij;
    j1 = j2 + 1;
  end
  i1 = i2 + 1;
end

C = ones(1,NeventTypes*Nh)/(NeventTypes);
vrf = 1/(C*inv(XtX)*C');

return;
%-----------------------------------------------

% Simulation to check %
TR = 1;
Ntrs = 100;
Trun = TR*Ntrs;
Nper = round(EventProb*Ntrs);
Tper = TR*ones(size(EventProb));
par = fmri_synthpar3(Nper,Tper,1,Trun,TR,0);
psdwin = [0 Nh*TR TR];

X = [];
for nthEventType = 1:NeventTypes
  ind = find(par(:,2) == nthEventType);
  st = par(ind,1);
  Xfir = fast_st2fir(st,Ntrs,TR,psdwin);
  X = [X Xfir];
end
XtX0 = (X'*X);
aXtX = XtX0/Ntrs;
rms = sqrt(mean((XtX(:)-aXtX(:)).^2));

avrf  = 1/(C*inv(aXtX)*C');
avrf0 = 1/(C*inv(XtX0)*C');
fprintf('N=%d, vrf = %g, avrf = %g, avrf0 = %g, rms = %g\n',...
	Ntrs,vrf,avrf,avrf0,rms);

figure(1);
imagesc(XtX);
colorbar;
showfigxy;

figure(2);
imagesc(aXtX);
colorbar;
showfigxy;

figure(3);
imagesc(abs(XtX-aXtX));
colorbar;
showfigxy;




