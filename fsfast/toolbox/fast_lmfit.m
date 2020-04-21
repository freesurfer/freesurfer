function [beta, exitflag, cost, beta2] = fast_lmfit(y,X,L)
% [beta exitflag cost beta2] = fast_lmfit(y,X,L)

if(nargin ~= 3)
  fprintf('[beta exitflag cost beta2] = fast_lmfit(y,X,L)');
  return;
end

%Using this inline was *slow*
%lmcost = inline('sum(abs(y - X*beta).^L)','beta','y','X','L');

%fprintf('Computing L2 solution\n');
beta2 = inv(X'*X)*X'*y;

vdof = size(X,1) - size(X,2);

nbeta = size(X,2);
nvox = size(y,2);
beta = zeros(nbeta,nvox);
exitflag = zeros(1,nvox);
cost = zeros(1,nvox);
n10 = round(nvox/100);
%fprintf('Starting voxel-wise iteration (%d)\n',nvox);
%keyboard
tic;
for nthvox = 1:nvox
  %fprintf('%4d \n',nthvox);
  if(nthvox == 0 | rem(nthvox,n10) == 0)
    pctdone = 100*nthvox/nvox;
    ttogo = (toc/nthvox)*(nvox-nthvox);
    fprintf('%4.1f %5d/%5d %8.2f %8.2f\n',...
	    pctdone,nthvox,nvox,toc/60,ttogo/60);
  end
  yv = y(:,nthvox);
  beta2v = beta2(:,nthvox);
  [betav costv exitflagv] = fminsearch(@lmcost,beta2v,[],yv,X,L);
  beta(:,nthvox) = betav;
  exitflag(1,nthvox) = exitflagv;
  cost(1,nthvox) = costv;
end

return;

%-------------------------------------------------%
function cost = lmcost(beta,y,X,L)
cost = sum(abs(y - X*beta).^L);
return;
