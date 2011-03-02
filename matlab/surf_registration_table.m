

%
% surf_registration_table.m
%
% Original Author: Laurence Wastiaux
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
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



[Dsubj, nf, r90, thdemean, M] = surf_registration_stats('/space/neo/2/recon/buckner')
%load('/space/okapi/3/data/laurence/ADF/surf_registration/residual_var90.mat') %loads r90
%load('/space/okapi/3/data/laurence/ADF/surf_registration/Subj_Names.mat') %loads nf
%files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/surfreg_table.txt', 'w');
fprintf(fp,'# FSAFD SurfRegistrationCheck 1\n# date 20050516\n# $Id: surf_registration_table.m,v 1.3 2011/03/02 00:04:13 nicks Exp $\n# Buckner data set\n# hemi lh\n# PercentRemovedVar 90\n# ncols 1\n# nrows %d\n# label_col 1 PercentResidualVariance\n', length(r90));
nn=1;
for i=1:length(r90)
    %nb=I(i);
    if(r90(i)>0)
        fprintf(fp, '# label_row %d %s\n', nn, char(nf(i).name));
        nn=nn+1;
    end
    
end
for j=1:length(r90)
    if(r90(j)>0)
        fprintf(fp, '%f\n', r90(j));
    end
end
fclose(fp)
