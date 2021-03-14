function [Xfir, Nc] = fast_par2Xfir(par,ntrs,TR,TER,TPreStim,TimeWindow,W)
% [Xfir Nc] = fast_par2Xfir(par,ntrs,TR,TER,TPreStim,TimeWindow,W)


%
% fast_par2Xfir.m
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
  psdwin = [TPreStim TER TimeWindow-TPreStim];
  Xfirc = fast_sched2Xfir(tPres,ntrs,TR,psdwin,Wc);
  Xfir = [Xfir Xfirc];
end

return;
