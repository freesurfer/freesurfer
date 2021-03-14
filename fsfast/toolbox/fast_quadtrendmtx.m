function Xqtrend = fast_quadtrendmtx(run,ntrs,nruns)
% Xqtrend = fast_quadtrendmtx(run,ntrs,nruns)
%
% Quadratic trend - centered half-way into the run; mean zero.
% This vector will be orthogonal to those produced by
% fast_baselinemtx fast_trendmtx 


%
% fast_quadtrendmtx.m
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

if(nargin ~= 3)
  msg = 'USAGE: Xqtrend = fast_quadtrendmtx(run,ntrs,nruns)';
  qoe(msg);error(msg);
end

t = [0:ntrs-1]'; %'
t = t - mean(t);
v = t.^2;
v = v - mean(v);
v = v./sqrt(sum(v.^2));

Xqtrend        = zeros(ntrs,nruns);
Xqtrend(:,run) = v;

return;
