function [roilist2, hemi] = striphemi(roilist,ncut)
% [roilist2 hemi] = striphemi(roilist,<ncut>)
% Removes hemifield designations from list of roi names
% hemi is a list equal to the roilist where 1=lh, 2=rh, 0=other
% If ncut is nonempty and greater than 0 then removes ncut letters
%   from the end of each roiname. This can be helpful when removing
%   things like '_thickness'.
% $Id: striphemi.m,v 1.4 2014/10/01 01:53:31 greve Exp $

if(~exist('ncut','var')) ncut = []; end
if(isempty(ncut)) ncut = 0; end

nrois = size(roilist,1);


lhhemistringlist = strvcat('Left-','ctx-lh-','wm-lh-','lh_');
nlhhemistrings = size(lhhemistringlist,1);
rhhemistringlist = strvcat('Right-','ctx-rh-','wm-rh-','rh_');
nrhhemistrings = size(rhhemistringlist,1);

lhhemistringlistpost = strvcat('_Left','_L'); % order important
nlhhemistringspost = size(lhhemistringlistpost,1);
rhhemistringlistpost = strvcat('_Right','_R'); % order important
nrhhemistringspost = size(rhhemistringlistpost,1);

hemi = zeros(nrois,1);
roilist2 = '';
for nth = 1:nrois
  s = deblank(roilist(nth,:));
  s = s(1:end-ncut);
  hit = 0;
  for nthhs = 1:nlhhemistrings
    hs = deblank(lhhemistringlist(nthhs,:));
    k = strfind(s,hs);
    if(~isempty(k))
      nh = length(hs);
      s = s(nh+1:end);
      roilist2 = strvcat(roilist2,s);
      hemi(nth) = 1;
      hit = 1;
      break;
    end
  end
  if(hit) continue; end
  for nthhs = 1:nrhhemistrings
    hs = deblank(rhhemistringlist(nthhs,:));
    k = strfind(s,hs);
    if(~isempty(k))
      nh = length(hs);
      s = s(nh+1:end);
      roilist2 = strvcat(roilist2,s);
      hemi(nth) = 2;
      hit = 1;
      break;
    end
  end
  if(hit) continue; end
  for nthhs = 1:nlhhemistringspost
    hs = deblank(lhhemistringlistpost(nthhs,:));
    k = strfind(s,hs);
    if(~isempty(k))
      nh = length(hs);
      s = s(1:end-nh);
      roilist2 = strvcat(roilist2,s);
      hit = 1;
      hemi(nth) = 1;
      break;
    end
  end
  if(hit) continue; end
  for nthhs = 1:nrhhemistringspost
    hs = deblank(rhhemistringlistpost(nthhs,:));
    k = strfind(s,hs);
    if(~isempty(k))
      nh = length(hs);
      s = s(1:end-nh);
      roilist2 = strvcat(roilist2,s);
      hit = 1;
      hemi(nth) = 2;
      break;
    end
  end
  if(hit) continue; end

  if(hit == 0) roilist2 = strvcat(roilist2,s); end
end

return


