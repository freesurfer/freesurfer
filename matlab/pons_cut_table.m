[Dpons, Ipons]=pons_cut_dir_afd('/space/neo/2/recon/buckner',0.01)
% load('/space/okapi/3/data/laurence/ADF/cutting_planes/Dice_pons.mat') %loads Dpons
% load('/space/okapi/3/data/laurence/ADF/cutting_planes/Subj_nb_pons.mat') %loads Ipons
files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/pons_table.txt', 'w');
fprintf(fp,'# FSAFD HorizontalCuttingPlaneCheck 1\n# date 20050516\n# $Id: pons_cut_table.m,v 1.1 2005/06/01 17:40:13 wastiaux Exp $\n# Buckner data set\n# ncols 1\n# nrows %d\n# label_col 1 DiceCoefficient\n', length(Dpons));
for i=1:length(Dpons)
    nb=Ipons(i);
    fprintf(fp, '# label_row %d %s\n', i, char(files(nb).name));
end
for j=1:length(Dpons)
    fprintf(fp, '%f\n', Dpons(j));
end
fclose(fp)