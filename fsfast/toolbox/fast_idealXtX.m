function [XtX C vrf] = fast_idealXtX(NrepsPer,TR,Ntrs,psdwin)
% [XtX C vrf] = fast_idealXtX(NrepsPer,TR,Ntrs,psdwin)
%
% Computes "ideal" FIR XtX, where "ideal" is the asymptotic XtX for an
% infinite length run with the given stimulus density (ie, NrepsPer/Ntrs).
% 
%
% NrepsPer - list of the number of repetitions of a non-null event types.
% psdwin = [psdmin psdmax dpsd];
%
% C = ones(1,Nstim*Nh)/(Nstim*Nh);
% vrf = 1/(C*inv(XtX)*C');
%
% There is a simulation at the end of this file (after return) to
% check the accuracy.
%

%
% fast_idealXtX
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

XtX = [];
if(nargin ~= 4)
  fprintf('[XtX C vrf] = fast_idealXtX(NrepsPer,TR,Ntrs,psdwin)\n');
  return;
end

psdmin  = psdwin(1);  % start of PSD window
psdmax  = psdwin(2);  % end of PSD window
dpsd    = psdwin(3);  % increment of PSD window
Nh      = round((psdmax-psdmin)/dpsd);
Rss     = round(TR/dpsd);
EventProb = NrepsPer/(Rss*Ntrs);

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

XtX = Ntrs*XtX;

C = dpsd*ones(1,NeventTypes*Nh)/(NeventTypes);
vrf = 1/(C*inv(XtX)*C');


return;
%-----------------------------------------------

% Simulation to check %
fprintf('Simulating\n');

Trun = TR*Ntrs;
Tper = TR*ones(size(EventProb));
par = fmri_synthpar3(NrepsPer,Tper,1,Trun,dpsd,0);

X = [];
for nthEventType = 1:NeventTypes
  ind = find(par(:,2) == nthEventType);
  st = par(ind,1);
  Xfir = fast_st2fir(st,Ntrs,TR,psdwin);
  X = [X Xfir];
end
XtX0 = (X'*X);
aXtX = XtX0;

rms = sqrt(mean((XtX(:)-aXtX(:)).^2));
avrf  = 1/(C*inv(aXtX)*C');
avrf0 = 1/(C*inv(XtX0)*C');
fprintf('N=%d, vrf = %g, avrf = %g, avrf0 = %g, rms = %g\n',...
	Ntrs,vrf,avrf,avrf0,rms);

figure(1);
imagesc(XtX);
colorbar;
showfigxy;
title('Ideal');

figure(2);
imagesc(aXtX);
colorbar;
showfigxy;
title('Actual');

figure(3);
imagesc(abs(XtX-aXtX));
colorbar;
showfigxy;
title('Diff');




