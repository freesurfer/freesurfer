[P, I, nf]=talairaching_dir_afd('/space/neo/2/recon/buckner', 0.01);
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_probas.mat') %loads y computed for the Buckner data set
%P=y;
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_names.mat') %loads nf (Buckner data set)
%files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/talairach_table.txt', 'w');
fprintf(fp,'# FSAFD TalairachingCheck 1\n# date 20050516\n# $Id: talairaching_table.m,v 1.1 2005/06/01 14:04:23 wastiaux Exp $\n# Buckner data set\n# ncols 1\n# nrows %d\n# label_col 1 xfmProbability\n', length(P));
for i=1:length(P)
    %u=I(i);
    %fprintf(fp, '# label_row %d %s\n', i, char(files(u).name));
    fprintf(fp, '# label_row %d %s\n', i, char(nf(i).name));
end
for j=1:length(P)
    fprintf(fp, '%f\n', P(j));
end
fclose(fp)