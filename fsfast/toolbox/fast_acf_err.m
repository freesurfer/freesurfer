function acferr = fast_acf_err(racf,R,yacf)
% acferr = fast_acf_err(racf,R,yacf)

if(nargin ~= 3)
  fprintf('acferr = fast_acf_err(racf,R,yacf)\n');
  return;
end

nf = length(racf);
Vy = toeplitz(yacf);

% Compute expected racf %
for l = 1:nf
  Dl = diag(ones(nf-l+1,1),l-1);  
  racfexp(l) = trace(R*Dl*R*Vy);
end
racfexp = racfexp/racfexp(1);




