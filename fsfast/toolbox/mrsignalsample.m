function [mrs, mrsstdexp] = mrsignalsample(S0mn,S0std,R2smn,R2sstd,Rawstd,TE,alpha,sz)
% [mrs, mrsstdexp] = mrsignalsample(S0mn,S0std,R2smn,R2sstd,TE,alpha,sz)


%
% mrsignalsample.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:34 $
%    $Revision: 1.2 $
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

if(~exist('sz','var')) sz = []; end
if(isempty(sz)) sz = 1; end

if(length(sz) == 1) sz = [sz 1]; end

S0  = S0mn  + S0std*randn(sz);
S0  = abs(S0);
R2s = R2smn + R2sstd*randn(sz);
R2s = abs(R2s);

nraw = Rawstd*randn(sz);

mrs = sin(alpha)*S0.*exp(-R2s*TE) + nraw;

mrsstdexp = sqrt( (S0std*exp(-TE*R2smn))^2 + ...
		  (R2sstd*TE*S0mn*exp(-TE*R2smn))^2 + ...
		  Rawstd);

return;



