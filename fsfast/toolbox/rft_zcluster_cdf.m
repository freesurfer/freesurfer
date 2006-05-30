function cdf = rft_zcluster_cdf(csize,vthresh,fwhm,ssize,D)
% cdf = rft_zcluster_cdf(csize,vthresh,fwhm,ssize,D)
%
% Prob that a cluster >= csize defined by threshold vthresh will be
% found in a D-dim z-field with fwhm smoothness of ssize search space
% 

cdf = [];
if(nargin ~= 5)
  fprintf('cdf = rft_zcluster_cdf(csize,vthresh,fwhm,ssize,D)\n');
  return;
end

W = fwhm/sqrt(4*log(2));
dLh = W.^(-D);
%dLh = 0.493971;

Em = ssize .* ((2*pi).^(-(D+1)/2)) .* dLh .* (vthresh.^(D-1)) .* ...
     exp(-(vthresh.^2)/2); 

beta = ( ((gamma(D/2+1) .* Em)) ./ (ssize.*mri_zcdf(-vthresh))).^(2./D);

cdf = exp(-beta.*(csize.^(2./D)));

%keyboard

return;


