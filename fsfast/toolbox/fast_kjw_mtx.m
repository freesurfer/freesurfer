function  M = fast_kjw_mtx(R,p)
% M = fast_kjw_mtx(R, <p>)
%
% Computes the matrix that remove bias from residual autocorrelation
% function using Keith Worsley's fmristats method. The bias is induced
% by projecting out the task components to form the residual.
%
% R is the residual forming matrix
% p - can't remember what this does, but it's optional
%
% Notes:
%  1. As the number of frames increases, the computation is 
%     more intense and the correction matrix is more singular.
%     There is some happy medium where the accuracy of the 
%     ACF is best, but it depends upon the design among other
%     things.
%
% See also: fast_acorr, fast_yacf_kjw.
%
%


%
% fast_kjw_mtx.m
%
% Original Author: Doug Greve
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if(nargin ~= 1 & nargin ~= 2)
  fprintf('M = fast_kjw_mtx(R, <p>)\n');
  return;
end

nf = size(R,1);
if(~exist('p')) p = nf; end

%fprintf('KJW Matrix: Stage 1\n');
Dl = eye(nf);
for l = 1:p
  %if(mod(l,5)==0 | l==1) fprintf('l = %d, %g\n',l,toc); end
  D(:,:,l) = Dl;
  DpDt(:,:,l) = Dl+Dl'; %'
  Dl = fast_mshift(Dl,[0 1],0);
end

%fprintf('KJW Matrix: Stage 2\n');
M = zeros(p,p);
for l = 1:p
  %if(mod(l,5)==0 | l==1) fprintf('l = %d, %g\n',l,toc); end
  Dl = D(:,:,l);
  RDlR = R*Dl*R;
  M(l,1) = trace(RDlR);
  for j = 2:p
    M(l,j) = trace(RDlR*DpDt(:,:,j)); 
  end
end
%fprintf('KJW Matrix: cond = %g  (%g)\n',cond(M),toc);


return;
