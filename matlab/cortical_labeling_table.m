


%
% cortical_labeling_table.m
%
% Original Author: Laurence Wastiaux
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.6 $
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



[Dl, Dr,I]=cortical_labeling_dir_afd('/space/neo/2/recon/buckner', 0.01)
% load('/space/okapi/3/data/laurence/ADF/cortical_labeling/PercentArea_labels_parctxt.mat') %loads D2_rh D2_lh
% Dl=D2_lh;
% Dr=D2_rh;
% load('/space/okapi/3/data/laurence/ADF/cortical_labeling/Subj_nb_parctxt.mat') %loads I
files=dir('/space/neo/2/recon/buckner');
labmapfile=('/space/lyon/1/fsdev/freesurfer_dev/Simple_surface_labels2002.txt');
[label name val1 val2 val3 val4]=textread(labmapfile,'%d %s %d %d %d %d',85);
fp=fopen('~/cortlab_table_lh.txt', 'w');
sz=size(Dl);
fprintf(fp,'# FSAFD CorticalLabelingCheck 1\n# date 20050516\n# $Id: cortical_labeling_table.m,v 1.6 2011/03/02 00:04:12 nicks Exp $\n# Buckner data set\n# info_file /label/lh.parc.txt\n# hemi lh\n# ncols 84\n# nrows %d\n', sz(1));
for u=1:(length(label)-1)
    fprintf(fp, '# label_col %d %s\n', u, char(name(u+1)));
end
for i=1:length(Dl)
    nb=I(i);
    fprintf(fp, '# label_row %d %s\n', i, char(files(nb).name));
end
for j=1:sz(1)
    for k=1:84
        fprintf(fp, '%f  ', Dl(j,k));
    end
    fprintf(fp, '\n');
end
fclose(fp)
fp=fopen('cortlab_table_rh.txt', 'w');
sz=size(Dr);
fprintf(fp2,'# FSAFD CorticalLabelingCheck 1\n# date 20050516\n# $Id: cortical_labeling_table.m,v 1.6 2011/03/02 00:04:12 nicks Exp $\n# Buckner data set\n# info_file /label/rh.parc.txt\n# hemi rh\n# ncols 84\n# nrows %d\n', sz(1));
for u=1:(length(label)-1)
    fprintf(fp, '# label_col %d %s\n', u, char(name(u+1)));
end
for i=1:length(Dr)
    nb=I(i);
    fprintf(fp, '# label_row %d %s\n', i, char(files(nb).name));
end
for j=1:sz(1)
    for k=1:84
        fprintf(fp2, '%f  ', Dr(j,k));
    end
    fprintf(fp2, '\n');
end
fclose(fp2)
