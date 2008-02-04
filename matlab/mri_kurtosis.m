function g2 = mri_kurtosis(x,unbiasedflag,dim)
% g2 = mri_kurtosis(x,<unbiasedflag>,<dim>)
% Kurtosis estimator. 
% Biased by default. For unbiased set unbiasedflag=1.
% To make estimator zero mean, subtract 3.
%
% $Id: mri_kurtosis.m,v 1.3 2008/02/04 22:48:46 greve Exp $

%
% mri_kurtosis.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2008/02/04 22:48:46 $
%    $Revision: 1.3 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%


g2 = [];
if(nargin < 1 & nargin > 3)
  fprintf('g2 = mri_kurtosis(x,<unbiasedflag>,<dim>)\n');
  return;
end

if(~exist('unbiasedflag','var')) unbiasedflag = []; end
if(isempty(unbiasedflag)) unbiasedflag = 0; end
if(~exist('dim','var')) dim = []; end
if(isempty(dim)) dim = 1; end

szx = size(x);
n = szx(dim);

xmean = mean(x,dim);
tmp = ones(1,length(szx));
tmp(dim) = n;
xmeanRep = repmat(xmean,tmp);

m4 = n*sum((x - xmeanRep).^4,dim);
m2 = sum((x - xmeanRep).^2,dim);

if(unbiasedflag)
  b1 = (n+1)*(n-1)/((n-2)*(n-3));
  b2 = ((n-1).^2)/((n-2)*(n-3));
  g2 = b1*(m4./(m2.^2)) - 3*b2 + 3;  
else
  % Biased
  g2 = m4./(m2.^2);
end

%G2 = b1*(m4./(m2.^2)) + 3*b2;

return;
