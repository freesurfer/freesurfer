function  M = fast_kjw_mtx(R,p)
% M = fast_kjw_mtx(R, <p>)

if(nargin ~= 1 & nargin ~= 2)
  fprintf('M = fast_kjw_mtx(R, <p>)\n');
  return;
end

nf = size(R,1);
if(~exist('p')) p = nf; end

M = zeros(p,p);

tic;
fprintf('KJW Matrix: Stage 1\n');
Dl = eye(nf);
for l = 1:p
  fprintf('l = %d, %g\n',l,toc);
  %Dl = diag(ones(nf-(l-1),1),l-1);  
  D(:,:,l) = Dl;
  DpDt(:,:,l) = Dl+Dl'; %'
  Dl = fast_mshift(Dl,[0 1],0);
end

fprintf('KJW Matrix: Stage 2\n');
M = zeros(p,p);
for l = 1:p
  fprintf('l = %d, %g\n',l,toc);
  %Dl = diag(ones(p-l,1),l);
  Dl = D(:,:,l);
  %RDl = fast_mshift(R,[0 l-1],0);
  RDlR = R*Dl*R;
  M(l,1) = trace(RDlR);
  for j = 2:p
    M(l,j) = trace(RDlR*DpDt(:,:,j)); 
  end
end



return;
