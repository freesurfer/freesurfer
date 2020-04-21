function [roibaselist, hemi, croilist] = striphemi(roilist,ncut)
% [roibaselist hemi] = striphemi(roilist,<ncut>)
% Removes hemifield designations from list of roi names
% hemi is a list equal to the roilist where 1=lh, 2=rh, 0=other
% croilist is a list of the contralateral roi names
% If ncut is nonempty and greater than 0 then removes ncut letters
%   from the end of each roiname. This can be helpful when removing
%   things like '_thickness'.

if(~exist('ncut','var')) ncut = []; end
if(isempty(ncut)) ncut = 0; end

nrois = size(roilist,1);

hemi = zeros(nrois,1);
roibaselist = '';
croilist = '';
for nth = 1:nrois
  s = deblank(roilist(nth,:));
  [hemi(nth) hemistr sbase croi] = roiname2hemi(s,ncut); 
  roibaselist = strvcat(roibaselist,sbase);
  croilist = strvcat(croilist,croi);
end

return

