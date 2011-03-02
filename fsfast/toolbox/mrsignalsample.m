function [mrs, mrsstdexp] = mrsignalsample(S0mn,S0std,R2smn,R2sstd,Rawstd,TE,alpha,sz)
% [mrs, mrsstdexp] = mrsignalsample(S0mn,S0std,R2smn,R2sstd,TE,alpha,sz)


%
% mrsignalsample.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:07 $
%    $Revision: 1.3 $
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



