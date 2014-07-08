function roilist2 = striphemi(roilist)
% roilist2 = striphemi(roilist)
% Removes hemifield designations from list of roi names
% $Id: striphemi.m,v 1.1 2014/07/08 16:03:25 greve Exp $

nrois = size(roilist,1);

hemistringlist = strvcat('Left-','Right-','ctx-lh-','ctx-rh-','wm-lh-','wm-rh-');
nhemistrings = size(hemistringlist,1);

roilist2 = '';
for nth = 1:nrois
  s = deblank(roilist(nth,:));
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


