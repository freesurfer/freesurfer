function acf = fast_rphacf(alpha,rho,T,nf)
% acf = fast_rphacf(alpha,rho,T,nf)

if(nargin < 1 | nargin > 4)
  fprintf('acf = fast_rphacf(alpha,rho,T,nf)\n');
  return;
end

tau = [0:nf-1]';
acf = alpha * (rho.^tau) .* cos(2*pi*tau/T);
acf(1) = 1;

return;


