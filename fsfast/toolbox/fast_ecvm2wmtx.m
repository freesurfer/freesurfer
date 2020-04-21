function [W AutoCor AutoCor0 ncnd] = fast_ecvm2wmtx(ErrCovMtx);
% [W AutoCor AutoCor0 ncnd] = fast_ecvm2wmtx(ErrCovMtx);
%   ErrCovMtx = (e*e')/size(e,2);
% W is the whitening matrix
% AutoCor is the conditioned autocorrelation function
% AutoCor0 is the unconditioned autocorrelation function
% ncnd is the number of conditioning steps.
%
% This implements the old selxavg-style method for computing the
% whitening matrix. This method regularizes by recursively tukey
% tapering the autocorrelation function until the covariance 
% matrix has a positive min eigenvalue. Note that it does not use
% taumax or pct (variables which are present in fast_selxavg.m
% but not actually used there either).
% 
% Note: this produced a non-causal whitening matrix, but that's
% probably unimportant in the end as W*S*W' is the same.
%
% Stanard way to do it:
%  A = toeplitz(AutoCor);
%  WA = inv(chol(A)');
%  Such that: WA*A*WA' is the identity.
%

if(nargin ~= 1)
  fprintf('[W AutoCor AutoCor0 ncnd] = fast_ecvm2wmtx(ErrCovMtx);\n');
  return;
end

fprintf('Using old selxavg-style whitening method\n');

nf = size(ErrCovMtx,1);
AutoCor0 = fast_cvm2acor(ErrCovMtx,1);
AutoCor = AutoCor0;
[cnd mineig S] = fast_acfcond(AutoCor);
ncnd = 1;
while(mineig < 0 | cnd > 100)
  fprintf('  ncnd = %d, mineig = %g, cnd = %g\n',ncnd,mineig,cnd);
  AutoCor = AutoCor .* tukeytaper(nf);
  [cnd mineig S] = fast_acfcond(AutoCor);
  ncnd = ncnd + 1;
end
fprintf('      ncnd = %d, mineig = %g, cnd = %g\n',ncnd,mineig,cnd);
[utmp stmp vtmp] = svd(S);
W = utmp * inv(sqrt(stmp)) * vtmp';
  


return;

