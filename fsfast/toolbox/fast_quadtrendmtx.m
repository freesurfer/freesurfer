function Xqtrend = fast_quadtrendmtx(run,ntrs,nruns)
% Xqtrend = fast_quadtrendmtx(run,ntrs,nruns)
%
% Quadratic trend - centered half-way into the run; mean zero.
% This vector will be orthogonal to those produced by
% fast_baselinemtx fast_trendmtx 

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
