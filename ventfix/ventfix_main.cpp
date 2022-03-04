#include <stdio.h>
#include <stdlib.h>

#include "ventfix.h"

// positioned parameters:
// argv[1] - subject
// argv[2] - aseg.presurf.mgz 
// argv[3] - segids (4,43)
// example: 
// setenv FREESURFER_SEED 2
// setenv SUBJECTS_DIR /autofs/space/sulc_001/users/ventfix/bigvent.subjects
// ./ventfix 0619 aseg.presurf.mgz 4,43
int main(int argc, char *argv[])
{
  char *subject = argv[1];
  char *asegps  = argv[2];
  char *segids  = argv[3];
  
  char mdir[256] = {'\0'}, asegpspath[256] = {'\0'}, brainmask[256] = {'\0'};
  char *subjdir = getenv("SUBJECTS_DIR");
  sprintf(mdir, "%s/%s/mri", subjdir, subject);
  sprintf(asegpspath, "%s/%s", mdir, asegps);
  sprintf(brainmask, "%s/brainmask.mgz", mdir);

  printf("read %s ...\n", asegpspath);
  MRI *asegVol = MRIread(asegpspath);
  if (asegVol == NULL)
  {
    printf("Error reading %s\n", asegpspath);
    return 1;
  }

  printf("read brainmask %s ...\n", brainmask);
  MRI *maskVol = MRIread(brainmask);
  if (maskVol == NULL)
  {
    printf("Error reading brainmask %s\n", brainmask);
    return 1;
  }

  MRI *newsegVol = VentFix::fixasegps(asegVol, maskVol, segids, 0.5, -1, 10000, 1);

  char newseg[256] = {'\0'};
  sprintf(newseg, "newseg.%s.mgz", subject);
  printf("Wrting to %s\n", newseg);
  MRIwrite(newsegVol, newseg);

  MRIfree(&asegVol);
  MRIfree(&maskVol);

  return 0;

#if 0
  char *subject = argv[1];
  char *asegps  = argv[2];
  float threshmin = atof(argv[3]);
  char *segids = argv[4];

  int niters = atoi(argv[5]);
  int nmax   = atoi(argv[6]);
  int topo   = atoi(argv[7]);

  char *newseg = argv[8];

  VentFix *vf = new VentFix();

  vf->init(NULL, subject, asegps, threshmin);
  vf->binarize();
  vf->apply_mask();
  vf->create_volcluster("ocn.sum.dat", "ocn.mgh");

  char *restsegids = (char*)segids;
  char *nextsegid = strtok_r(segids, ",", &restsegids);
  while (nextsegid != NULL)
  {
    int segid = atoi(nextsegid);
    nextsegid = strtok_r(NULL, ",", &restsegids);
    printf("\nexpandSegIndices for segid %d\n", segid);
    if (nextsegid != NULL)
      vf->expandSegIndices(segid,  niters, nmax, topo, "ventfill.lh.dat",  NULL,  "ventfill.lh.json");
    else
      vf->expandSegIndices(segid, niters, nmax, topo, "ventfill.rh.dat", newseg, "ventfill.rh.json");
  }
  vf->reset();
#endif

  return 0;
}
