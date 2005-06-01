[Dsubj, nf, r90, thdemean, M] = surf_registration_stats('/space/neo/2/recon/buckner')
%load('/space/okapi/3/data/laurence/ADF/surf_registration/residual_var90.mat') %loads r90
%load('/space/okapi/3/data/laurence/ADF/surf_registration/Subj_Names.mat') %loads nf
%files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/surfreg_table.txt', 'w');
fprintf(fp,'# FSAFD SurfRegistrationCheck 1\n# date 20050516\n# $Id: surf_registration_table.m,v 1.1 2005/06/01 17:57:43 wastiaux Exp $\n# Buckner data set\n# hemi lh\n# PercentRemovedVar 90\n# ncols 1\n# nrows %d\n# label_col 1 PercentResidualVariance\n', length(r90));
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