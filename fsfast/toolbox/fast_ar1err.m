function [err, Lmin, rhobest] = fast_ar1err(rho,R,Mr)
% [err Lmin] = fast_ar1err(rho,R,Mr)
%
% rho is a list of numbers between 0 and 1
% R is the resdiual error forming matrix before whitening
% Mr = r*r'/nv; %' where r is the residual error before whitening
% 
%  acf = r.^nnf;
%  L = toeplitz(acf);
%  err(n) = sum(reshape1d(Mr-R*L*R).^2);  
%
%  Lmin is the best from the list.

err=[];
Lmin=[];

if(nargin ~= 3)
  fprintf('[err Lmin] = fast_ar1err(rho,R,Mr)\n');
  return;
end

nf = size(R,1);
nnf = 0:nf-1;

if(size(Mr,1) ~= nf)
  fprintf('ERROR: size(Mr,1) ~= nf)\n');
  return;
end

n=1;
for r = rho
  acf = r.^nnf;
  L = toeplitz(acf);
  err(n) = sum(reshape1d(Mr-R*L*R).^2);  
  if(n==1)
    errmin = err(n);
    Lmin = L;
    rhobest = r;
  else
    if(errmin > err(n))
      Lmin = L;
      errmin = err(n);
      rhobest = r;
    end
  end
  n = n+1;
end
