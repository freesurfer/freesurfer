function [rho, nperslice, refslice] = fast_asorder(func,mask)
% [rho nperslice refslice] = fast_asorder(func,mask)

rho = [];
nperslice = [];
refslice = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('[rho nperslice refslice] = fast_asorder(func,mask)\n');
  return;
end

if(ischar(func))
  f = MRIread(func);
  if(isempty(f)) return; end
  func = f.vol;
  clear f;
end

[nr nc nslices nf] = size(func);

if(nargin == 1) mask = []; end
if(isempty(mask))  mask = ones(nr,nc,nslices); end

for nthslice = 1:nslices
  fs = fast_vol2mat(func(:,:,nthslice,:));
  ind = find(mask(:,:,nthslice));
  nperslice(nthslice) = length(ind);
  if(nperslice(nthslice) ~= 0) 
    % mean timecourse across nth slice
    fmn(:,nthslice) = mean(fs(:,ind),2); 
  else
    fmn(:,nthslice) = randn(nf,1); % hack
  end
end

X = fast_polytrendmtx(1,nf,1,2);
R = eye(nf) - X*inv(X'*X)*X';
rmn = R*fmn;
rmsse = sqrt(sum(rmn.^2));
rmnnorm = rmn./repmat(rmsse,[nf 1]);

refslice = round(nslices/2);
tcref = repmat(rmnnorm(:,refslice),[1 nslices]);
rho = sum(rmnnorm.*tcref);

return;




