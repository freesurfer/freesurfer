function [yacf, M] = fast_yacf_kjw(racf,R)
% yacf = fast_yacf_kjw(racf,R)

if(nargin ~= 2)
  fprintf('yacf = fast_yacf_kjw(racf,R)\n');
  return;
end

nf = length(racf);
M = zeros(nf,nf);

tic;
Dl = eye(nf);
for l = 1:nf
  %Dl = diag(ones(nf-(l-1),1),l-1);  
  D(:,:,l) = Dl;
  DpDt(:,:,l) = Dl+Dl'; %'
  Dl = fast_mshift(Dl,[0 1],0);
end


M = zeros(nf,nf);
for l = 1:nf
  fprintf('l = %d, %g\n',l,toc);
  %Dl = diag(ones(nf-l,1),l);
  Dl = diag(ones(nf-(l-1),1),l-1);  
  %RDl = fast_mshift(R,[0 l-1],0);
  RDlR = R*Dl*R;
  M(l,1) = trace(RDlR);
  for j = 2:nf
    M(l,j) = trace(RDlR*DpDt(:,:,j)); 
  end
end

%tmp = racf .* [nf:-1:1]'; %'
%yacf = inv(M)*tmp;
yacf = inv(M)*racf;
yacf = yacf/yacf(1);
toc


return;
