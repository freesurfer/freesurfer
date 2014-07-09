function roilist2 = striphemi(roilist,ncut)
% roilist2 = striphemi(roilist,<ncut>)
% Removes hemifield designations from list of roi names
% If ncut is nonempty and greater than 0 then removes ncut letters
%   from the end of each roiname. This can be helpful when removing
%   things like '_thickness'.
% $Id: striphemi.m,v 1.2 2014/07/09 23:17:13 greve Exp $

if(~exist('ncut','var')) ncut = []; end
if(isempty(ncut)) ncut = 0; end

nrois = size(roilist,1);

hemistringlist = strvcat('Left-','Right-','ctx-lh-','ctx-rh-',...
			 'wm-lh-','wm-rh-','rh_','lh_');
nhemistrings = size(hemistringlist,1);

roilist2 = '';
for nth = 1:nrois
  s = deblank(roilist(nth,:));
  s = s(1:end-ncut);
  hit = 0;
  for nthhs = 1:nhemistrings
    hs = deblank(hemistringlist(nthhs,:));
    k = strfind(s,hs);
    if(~isempty(k))
      nh = length(hs);
      s = s(nh+1:end);
      roilist2 = strvcat(roilist2,s);
      hit = 1;
      break;
    end
  end
  if(hit == 0) roilist2 = strvcat(roilist2,s); end
end

return


