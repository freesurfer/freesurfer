[D, pvals, I]=wm_seg_dir_afd('/space/neo/2/recon/buckner', 0.01);
%load('/space/okapi/3/data/laurence/ADF/wm_segmentation/WM_dice_buckner.mat') %loads D
%load('/space/okapi/3/data/laurence/ADF/wm_segmentation/subject_nb.mat') %loads I
%files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/wm_table.txt', 'w');
fprintf(fp,'# FSAFD WMsegmentationCheck 1\n# date 20050516\n# $Id: wm_seg_table.m,v 1.1 2005/06/01 14:41:54 wastiaux Exp $\n# Buckner data set\n# ncols 1\n# nrows %d\n# label_col 1 DiceCoefficient\n', length(D));
for i=1:length(D)
    nb=I(i);
    fprintf(fp, '# label_row %d %s\n', i, char(files(nb).name));
end
for j=1:length(D)
    fprintf(fp, '%f\n', D(j));
end
fclose(fp)