function [rmse, racfest] = fast_rphacfrmse(params,racf,R,taper)
% [rmse, racfest] = fast_rphacfrmse(params,racf,R,<taper>)

if(nargin < 1 | nargin > 4)
  fprintf('[rmse racfest]= fast_rphacfrmse(params,racf,R,<taper>)\n');
  return;
end

alpha = params(1);
rho   = params(2);
T     = params(3);

nf = length(racf);
tau = [0:nf-1]';

nacf = fast_rphacf(alpha,rho,T,nf);
racfest = fast_cvm2acor(R*toeplitz(nacf)*R,1);

if(nargin ~= 4)
  rmse = mean( (racf-racfest).^2 );
  rmseb = mean( (racf-nacf).^2 );
  rmse = rmse + rmseb/4;
else
  rmse = mean( ((racf-racfest).*taper).^2 );
end

return;