#!/bin/tcsh -ef

set echo=1
set w=2
./mri_train_autoencoder   -s 1 -dt 1e-2  -m 0.5 -x 40 160 40 160 40 160 -tol 1e-6 -w $w ~/links/subjects/Long/bruce/mri/norm.mgz ./bruce.w${w}.ae
#./mri_train_autoencoder -b 1 .2   -s 1 -dt 1e-2  -m 0.5 -x 100 150 100 150 100 150 -tol 1e-6 -w $w ~/links/subjects/Long/bruce/mri/norm.mgz ./bruce.w${w}.ae
#./mri_train_autoencoder  -s 1 -dt 1e-2  -m 0.1 -b 1 .5 -c -x 100 150 100 150 100 150 -tol 1e-6 -w $w ~/links/subjects/Long/bruce/mri/norm.mgz ./bruce.w${w}.ae
#./mri_train_autoencoder -w $w ~/links/subjects/Long/bruce/mri/norm.mgz ./bruce.w${w}.ae
