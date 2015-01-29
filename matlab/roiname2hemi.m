function [hemi, hemistr, roibasename, roicontraname] = roiname2hemi(roiname,ncut)
% [hemi, hemistr, roibasename, roicontraname] = roiname2hemi(roiname,<ncut>)
% hemi = 1 for left, hemistr = lh
% hemi = 2 for right, hemistr = rh
% hemi = 0 for neither, hemistr = na (contra=roiname)
% roibasename is roiname with hemisphere designation removed
% roicontraname is name of the contralateral roi (no cut)
% ncut - remove ncut characters from the end of roiname before
%   determining hemi
% If roiname is a string array, then it is run on each member

if(nargin ~= 1 & nargin ~= 2) 
  fprintf('[hemi, hemistr, roibasename, roicontraname] = roiname2hemi(roiname,<ncut>)\n');
  return;
end

if(~exist('ncut','var')) ncut = []; end
if(isempty(ncut)) ncut = 0; end

if(size(roiname,1) > 1)
  hemi = zeros(size(roiname,1),1);
  hemistr = '';
  roibasename = '';
  roicontraname = '';
  for n = 1:size(roiname,1)
    [hemin hemistrn roibasenamen roicontranamen] = roiname2hemi(roiname(n,:),ncut);
    hemi(n) = hemin; 
    hemistr = strvcat(hemistr,hemistrn);
    roibasename = strvcat(roibasename,roibasenamen);
    roicontraname = strvcat(roicontraname,roicontranamen);
  end
  return;
end

% order of lhhemilist must be same as rhhemilist
lhhemistringlist = strvcat('Left-','Left_','ctx-lh-','wm-lh-','lh_');
nlhhemistrings = size(lhhemistringlist,1);
rhhemistringlist = strvcat('Right-','Right_','ctx-rh-','wm-rh-','rh_');
nrhhemistrings = size(rhhemistringlist,1);

lhhemistringlistpost = strvcat('_Left','_L'); % order important
nlhhemistringspost = size(lhhemistringlistpost,1);
rhhemistringlistpost = strvcat('_Right','_R'); % order important
nrhhemistringspost = size(rhhemistringlistpost,1);

s = roiname(1:end-ncut);

for nthhs = 1:nlhhemistrings
  hs = deblank(lhhemistringlist(nthhs,:));
  k = strfind(s,hs);
  if(~isempty(k))
    nh = length(hs);
    roibasename = s(nh+1:end);
    hemi = 1;
    hemistr = 'lh';
    chs = deblank(rhhemistringlist(nthhs,:));
    roicontraname = sprintf('%s%s',chs,roiname(nh+1:end));
    return;
  end
end

for nthhs = 1:nrhhemistrings
  hs = deblank(rhhemistringlist(nthhs,:));
  k = strfind(s,hs);
  if(~isempty(k))
    nh = length(hs);
    roibasename = s(nh+1:end);
    hemi = 2;
    hemistr = 'rh';
    chs = deblank(lhhemistringlist(nthhs,:));
    roicontraname = sprintf('%s%s',chs,roiname(nh+1:end));
    return;
  end
end

for nthhs = 1:nlhhemistringspost
  hs = deblank(lhhemistringlistpost(nthhs,:));
  k = strfind(s,hs);
  if(~isempty(k))
    nh = length(hs);
    roibasename = s(1:end-nh);
    hemi = 1;
    hemistr = 'lh';
    chs = deblank(rhhemistringlistpost(nthhs,:));
    roicontraname = sprintf('%s%s',roiname(1:end-nh),chs);
    return;
  end
end

for nthhs = 1:nrhhemistringspost
  hs = deblank(rhhemistringlistpost(nthhs,:));
  k = strfind(s,hs);
  if(~isempty(k))
    nh = length(hs);
    roibasename = s(1:end-nh);
    hemi = 2;
    hemistr = 'rh';
    chs = deblank(lhhemistringlistpost(nthhs,:));
    roicontraname = sprintf('%s%s',roiname(1:end-nh),chs);
    return;
  end
end

roibasename = s;
roicontraname = roiname;
hemi = 0;
hemistr = 'na';


