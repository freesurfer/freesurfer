
%
% cc_cut_table.m
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



[Dr, Dl, I]=cc_cut_dir_afd('/space/neo/2/recon/buckner', 0.01)
% load('/space/okapi/3/data/laurence/ADF/cutting_planes/Dice_cc_lrh_lta.mat') %loads Dr Dl
% load('/space/okapi/3/data/laurence/ADF/cutting_planes/Subj_nb_cc.mat') %loads I
files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/cc_lh_table.txt', 'w');
fprintf(fp,'# FSAFD VerticalCuttingPlaneCheck 1\n# date 20050516\n# $Id: cc_cut_table.m,v 1.3 2011/03/02 00:04:12 nicks Exp $\n# Buckner data set\n# hemi lh\n# ncols 1\n# nrows %d\n# label_col 1 DiceCoefficient\n', length(Dl));
for i=1:length(Dl)
    nb=I(i);
    fprintf(fp, '# label_row %d %s\n', i, char(files(nb).name));
end
for j=1:length(Dl)
    fprintf(fp, '%f\n', Dl(j));
end
fclose(fp)
fp2=fopen('~/cc_rh_table.txt', 'w');
fprintf(fp2,'# FSAFD VerticalCuttingPlaneCheck 1\n# date 20050516\n# $Id: cc_cut_table.m,v 1.3 2011/03/02 00:04:12 nicks Exp $\n# Buckner data set\n# hemi rh\n# ncols 1\n# nrows %d\n# label_col 1 DiceCoefficient\n', length(Dr));
for i=1:length(Dr)
    nb=I(i);
    fprintf(fp2, '# label_row %d %s\n', i, char(files(nb).name));
end
for j=1:length(Dr)
    fprintf(fp2, '%f\n', Dr(j));
end
fclose(fp2)
