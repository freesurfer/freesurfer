function [yacf, M] = fast_yacf_kjw(racf,R)
% yacf = fast_yacf_kjw(racf,R)

if(nargin ~= 2)
  fprintf('yacf = fast_yacf_kjw(racf,R)\n');
  return;
end

[nf nv] = size(racf);
M = fast_kjw_mtx(R,nf);

yacf = inv(M)*racf;
yacf = yacf./repmat(yacf(1,:),[nf 1]);


toc


return;
