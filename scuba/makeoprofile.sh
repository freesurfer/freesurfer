#! /bin/sh

mv -f profile-samples.txt profile-samples2.txt
mv -f profile-annotated.txt profile-annotated2.txt

sudo opreport -l /home/kteich/fsdev/trunk/dev/scuba/scuba > profile-samples.txt
sudo opannotate --source /home/kteich/fsdev/trunk/dev/scuba/scuba > profile-annotated.txt
