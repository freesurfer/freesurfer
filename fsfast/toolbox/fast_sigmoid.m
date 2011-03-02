function y = fast_sigmoid(x,dx,s);
% y = fast_sigmoid(x,<dx>,<s>);
%
% y passes thru 0.5 at x=dx. Default dx=0.
% s controls the slope. Default s=1.
% dx and s can either be scalars or same size as x.
%
% y = 0 at x = -inf 
% y = 1 at x = +inf 
% y will be the same size as x. 
%
% $Id: fast_sigmoid.m,v 1.2 2011/03/02 00:04:05 nicks Exp $

%
% fast_sigmoid.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:05 $
%    $Revision: 1.2 $
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

y = [];
if(nargin < 0 | nargin > 3)
  fprintf('y = fast_sigmoid(x,<dx>,<s>);\n');
  return;
end

if(~exist('dx','var')) dx = []; end
if(isempty(dx)) dx = 0; end

if(~exist('s','var')) s = []; end
if(isempty(s)) s = 1; end

y = 1./(1 + exp(-s.*(x-dx)));

return;
