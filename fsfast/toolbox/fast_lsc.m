function lsc = fast_lsc(y)
% lsc = fast_lsc(y)


%
% fast_lsc.m
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

lsc = [];
if(nargin ~= 1)
  fprintf('lsc = fast_lsc(y)\n');
  return;
end

[nr nc ns nf] = size(y);

% normalize
yss = sqrt(sum(y.^2,4));
y = y./repmat(yss,[1 1 1 nf]);

lsc = zeros(nr,nc,ns,6);
%tic;
for r = 2:nr-1
  %fprintf('r = %d, %g\n',r,toc);
  for c = 2:nc-1
    for s = 2:ns-1
      y0 = squeeze(y(r,c,s,:));

      yn = squeeze(y(r-1,c,s,:));
      lsc(r,c,s,1) = sum(y0.*yn);
      yn = squeeze(y(r+1,c,s,:));
      lsc(r,c,s,2) = sum(y0.*yn);

      yn = squeeze(y(r,c-1,s,:));
      lsc(r,c,s,3) = sum(y0.*yn);
      yn = squeeze(y(r,c+1,s,:));
      lsc(r,c,s,4) = sum(y0.*yn);

      yn = squeeze(y(r,c,s-1,:));
      lsc(r,c,s,5) = sum(y0.*yn);
      yn = squeeze(y(r,c,s+1,:));
      lsc(r,c,s,6) = sum(y0.*yn);
    
    end % ds
  end % dc
end % dr


return;



