function Xbaseline = fast_baselinemtx(run,ntrs,nruns)
% Xbaseline = fast_baselinemtx(run,ntrs,nruns)

if(nargin ~= 3)
  msg = 'USAGE: Xbaseline = fast_baselinemtx(run,ntrs,nruns)';
  qoe(msg);error(msg);
end

v = ones(ntrs,1);

Xbaseline        = zeros(ntrs,nruns);
Xbaseline(:,run) = v;

return;
