function sratio = fast_sratio(num,den)
% sratio = fast_sratio(num,den)
% Computes a "signed" ratio, ie, 
%  if num > den then sratio = +num/den
%  if num < den then sratio = -den/num
%

%
% fast_fast_sratio.m
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

if(nargin ~= 2)
  fprintf('sratio = fast_sratio(num,den)\n');
  return;
end

mri = [];
if(isfield(num,'vol'))  
  num = num.vol; 
  mri = num;
end
if(isfield(den,'vol'))  
  mri = den;
  den = den.vol; 
end
 
sratio = zeros(size(num));
ind = find(num > den & den ~= 0);
sratio(ind) = num(ind)./den(ind);
ind = find(num < den & num ~= 0);
sratio(ind) = -den(ind)./num(ind);

if(~isempty(mri))
  mri.vol = sratio;
  sratio = mri;
end


return;


