

%
% pons_cut_table.m
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



[Dpons, Ipons]=pons_cut_dir_afd('/space/neo/2/recon/buckner',0.01)
% load('/space/okapi/3/data/laurence/ADF/cutting_planes/Dice_pons.mat') %loads Dpons
% load('/space/okapi/3/data/laurence/ADF/cutting_planes/Subj_nb_pons.mat') %loads Ipons
files=dir('/space/neo/2/recon/buckner');
fp=fopen('~/pons_table.txt', 'w');
fprintf(fp,'# FSAFD HorizontalCuttingPlaneCheck 1\n# date 20050516\n# $Id: pons_cut_table.m,v 1.3 2011/03/02 00:04:12 nicks Exp $\n# Buckner data set\n# ncols 1\n# nrows %d\n# label_col 1 DiceCoefficient\n', length(Dpons));
for i=1:length(Dpons)
    nb=Ipons(i);
    fprintf(fp, '# label_row %d %s\n', i, char(files(nb).name));
end
for j=1:length(Dpons)
    fprintf(fp, '%f\n', Dpons(j));
end
fclose(fp)
