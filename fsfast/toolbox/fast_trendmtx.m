function Xtrend = fast_trendmtx(run,ntrs,nruns)
% Xbasline = fast_trendmtx(run,ntrs,nruns)

if(nargin ~= 3)
  msg = 'USAGE: Xbasline = fast_trendmtx(run,ntrs,nruns)';
  qoe(msg);error(msg);
end

v = [0:ntrs-1]'; %'
v = v - mean(v);
v = v./sqrt(sum(v.^2));

Xtrend        = zeros(ntrs,nruns);
Xtrend(:,run) = v;

return;
