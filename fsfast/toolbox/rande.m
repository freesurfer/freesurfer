function r = rande(rdim,order)
%
% r = rande(dim,<order>)
%
% Generates matrix of given dimension whose elements 
% have an erlang distribution. The expected value
% is 1. The variance is 1/r. The stddev is 1/sqrt(r).
%
% If the order is unspecified, defalts to 1.
% An order of 1 is an exponential distribution.
%
% r = rande([5 50],3);
%
% See also erlang.m
%
%


%
% rande.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:07 $
%    $Revision: 1.5 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

if(exist('rdim')  ~= 1) rdim  = 1; end
if(exist('order') ~= 1) order = 1; end

r = 0;
for n = 1:order
  r = r + -log(rand(prod(rdim),1));
end
r = r/order; % make avg=1
r = reshape(r, [rdim 1]);

return;
