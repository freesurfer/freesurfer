function [gcomp, theta, thetahat, Fsig, beta, Rrow] = tdr_ghostcomp(kpcn,epipar,order)
% [gcomp theta thetahat Fsig] = tdr_ghostcomp(kpcn,epipar,<order>)
%
% kpcn is (3,ncols,n1,n2,n3) with the rows NOT LR reversed
% gcomp is (ncols n1 n2 n3) 
% order is polynomial order (def is 2)
%
% revenfixed = reven .* gcomp (careful with transpose)
%

%order = 2;   % poly order
rthresh = .2; % keep the top rthresh

gcomp = [];
if(nargin < 2 | nargin > 3)
  fprintf('[gkcomp theta thetahat Fsig] = tdr_ghostcomp(kpcn,epipar,<order>)\n');
  return;
end

if(~exist('order','var')) order = []; end
if(isempty(order)) order = 2; end

[npcn ncols n1 n2 n3] = size(kpcn);
if(npcn ~= 3) 
  fprintf('ERROR: npcn = %d, should be 3)\n',npcn);
  return;
end
nsamples = n1*n2*n3;

if(isstruct(epipar))
  [kvec gvec] = kspacevector2(epipar);
  % Compute the Ideal col and row DFT reconstruction matrices
  Frow = fast_dftmtx(kvec);
  Frow = fast_svdregpct(Frow,90);
  Rrow = inv(Frow); % Dont do transpose, mult vects on right
else
  Rrow = epipar;
  Frow = inv(Rrow'*Rrow)*Rrow';
end


% Flip even rows left-right
kpcn(2,:,:,:) = flipdim(kpcn(2,:,:,:),2);

% kpcn is (3,ncols,n1,n2,n3)
kpcn = reshape(kpcn,[3 ncols n1*n2*n3]);
kpcn = permute(kpcn,[2 3 1]);

% Reconstruct the rows
pcn1 = Rrow*kpcn(:,:,1);
pcn2 = Rrow*kpcn(:,:,2);
pcn3 = Rrow*kpcn(:,:,3);

pcnref = (pcn1+pcn3)/2;
pcnmov = pcn2;

for nth = 1:nsamples
  [thetahat(:,nth) theta(:,nth) Fsig(:,nth) beta] = ...
      tdr_phdist_est0(pcnref(:,nth),pcnmov(:,nth),rthresh,order);
end

% Use the same for everyone
tmp = mean(thetahat,2);
thetahat = repmat(tmp,[1 nsamples]);

gcomp = exp(+i*thetahat);

% Apply compensation to the pcnmov
pcnmovcomp = pcnmov .* gcomp;

% Convert it back to kspace
kpcnmovcomp = Frow*pcnmovcomp;

% The compensated should look very close to the ref
kpcnref = (kpcn(:,:,1) + kpcn(:,:,3))/2;
kpcnmov = kpcn(:,:,2);
abskpcnref = (abs(kpcn(:,:,1)) + abs(kpcn(:,:,3)))/2;
abskpcnmov = abs(kpcn(:,:,2));
nn = 1:ncols;
%plot(nn,abskpcnref(:,1),'+-',nn,abskpcnmov(:,1),nn,abs(kpcnmovcomp(:,1)));
%legend('ref','mov','movcomp')

gcomp = reshape(gcomp,[ncols n1 n2 n3]);

gcompmov = exp(+i*thetahat/2);
gcompref = exp(-i*thetahat/2);


return;

%--------------------------------------------------------------------------
function [thetahat, theta, Fsig, beta] = tdr_phdist_est0(Ref,Mov,rthresh,order)

thetahat = [];
theta = [];
Fsig = [];
beta  = [];

if(nargin < 2 | nargin > 4)
  fprintf('[thetahat theta Fsig beta] = tdr_phdist_est(Ref,Mov,rthresh,order)\n');
  return;
end

if(exist('rthresh') ~= 1) rthresh = []; end 
if(isempty(rthresh)) rthresh = 2; end 
if(exist('order') ~= 1) order = []; end 
if(isempty(order)) order = 2; end 
  
if(max(abs(size(Ref)-size(Mov))) ~= 0)
  fprintf('ERROR: size mismatch between Ref and Mov\n');
  return;
end
[nf nv] = size(Ref);

% Average image space complex
RefMn = mean(abs(Ref),2);
MovMn = mean(abs(Mov),2);
RefMovMn = (RefMn+MovMn)/2;

% Compute the absolute segmentation threshold %
%mn = mean(RefMovMn);
%athresh = rthresh * mn; 

s = sort(RefMovMn);
athresh = s(round((1-rthresh)*nf)); % Keep top rthresh percent

% Segment - all columns get the same segmentation %
indover = find(RefMovMn > athresh);
nover = length(indover);
if(nover == 0) 
  fprintf('ERROR: no voxels over threshold\n');
  return;
end
indunder = find(RefMovMn <= athresh);

% Compute the phase distortion %
%theta = unwrap(angle(Ref./Mov));
theta = angle(Ref./Mov);

% Create a polynomial design matrix %
X = fast_polytrendmtx(1,nf,1,order); 

% Use only the suprathreshold frames for estimation %
thetaover = unwrap(theta(indover,:));
Xover = X(indover,:);

% Perform the estimation on suprathreshold frames %
[beta rvar vdof] = fast_glmfit(thetaover,Xover);

% Compute the estimate of the phase distortion across all frames
thetahat = X*beta;

% Make thetahat pass thru 0 at first sample
% thetahat = thetahat - repmat(thetahat(1,:),[nf 1]);

% Do a test to make sure things are ok
C = zeros(1,order+1);
C(2) = 1;
[F Fsig ces] = fast_fratio(beta,X,rvar,C);
Fsig = -log10(Fsig);
if(min(Fsig) < 200)
  fprintf('WARNING: tdr_phdist_est(): min(Fsig) = %g < 200  ',min(Fsig));
  fprintf('(threshold is %g)\n',rthresh);
end

%plot(1:nf,theta,1:nf,thetahat)
%keyboard

return;
