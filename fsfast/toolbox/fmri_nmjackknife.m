function a = fmri_nmjackknife(v)
% Jackknife averaging with normalization
%
% a = fmri_nmjackknife(v)
%
% size(v) = (Nh,Nvox,Nsamples)


%
% fmri_nmjackknife.m
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

[Nh Nvox Nsamples] = size(v);

for s = 1:Nsamples
  jk = find([1:Nsamples] ~= s);
  vtemplate  = fmri_norm(mean(v(:,:,jk),3),2);
  %vtemplate  = fmri_norm(randn(size(v(:,:,1))),2);

  if(size(v,1) > 1)
    a(s,:) = sum(v(:,:,s).*vtemplate);
  else
    a(s,:) = v(:,:,s).*vtemplate;
  end

end

return;
