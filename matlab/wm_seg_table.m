[D, pvals, I]=wm_seg_dir_afd('/space/neo/2/recon/buckner', 0.01);


%
% wm_seg_table.m
%
% Original Author: Laurence Wastiaux
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


%load('/space/okapi/3/data/laurence/ADF/wm_segmentation/WM_dice_buckner.mat') %loads D
%load('/space/okapi/3/data/laurence/ADF/wm_segmentation/subject_nb.mat') %loads I
%files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/wm_table.txt', 'w');
fprintf(fp,'# FSAFD WMsegmentationCheck 1\n# date 20050516\n# $Id: wm_seg_table.m,v 1.3 2011/03/02 00:04:13 nicks Exp $\n# Buckner data set\n# ncols 1\n# nrows %d\n# label_col 1 DiceCoefficient\n', length(D));
for i=1:length(D)
    nb=I(i);
    fprintf(fp, '# label_row %d %s\n', i, char(files(nb).name));
end
for j=1:length(D)
    fprintf(fp, '%f\n', D(j));
end
fclose(fp)
