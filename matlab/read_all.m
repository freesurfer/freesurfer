

%
% read_all.m
%
% Original Author: Bruce Fischl
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


names = str2mat('lh.thickness_alex.asc', 'lh.thickness_anders.asc', 'lh.thickness_arthur.asc','lh.thickness_buckner.asc','lh.thickness_chris.asc','lh.thickness_cindy.asc','lh.thickness_comite.asc','lh.thickness_gray.asc','lh.thickness_holmes.asc','lh.thickness_janine.asc','lh.thickness_jennifer.asc','lh.thickness_jonas.asc','lh.thickness_julieg.asc','lh.thickness_margaret.asc','lh.thickness_menell.asc','lh.thickness_moore.asc','lh.thickness_nancy.asc','lh.thickness_nava.asc','lh.thickness_nkarp.asc','lh.thickness_nouchine.asc','lh.thickness_paula.asc','lh.thickness_reppas.asc','lh.thickness_rhodes.asc','lh.thickness_roger.asc','lh.thickness_seppo.asc','lh.thickness_sherri.asc','lh.thickness_somers.asc','lh.thickness_stephan.asc','lh.thickness_talavage.asc','lh.thickness_tamily.asc','lh.thickness_tony.asc','lh.thickness_wim.asc') ; 


[rows,cols] = size(names);

for i=1:rows
		sprintf('reading curvature file %s\n', names(i,:)) 
		curv(:,i) = read_curv(names(i,:)) ;
end

ata = curv'*curv ;
[v,s2,vt] = svd(ata) ;
s = sqrt(s2) ;
u = curv*v*inv(s) ;
