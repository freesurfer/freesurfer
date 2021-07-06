#include <iostream>
#include "Spline.h"
#include "mri.h"
#include "label.h"

void set_voxel_value(MRI* mri, double dx, double dy, double dz, double val)
{
  int x = rint(dx), y = rint(dy), z = rint(dz);
  switch ( mri->type )
  {
  case MRI_UCHAR:
    MRIseq_vox( mri, x, y, z, 0 ) = (unsigned char)val;
    break;
  case MRI_INT:
    MRIIseq_vox( mri, x, y, z, 0 ) = (int)val;
    break;
  case MRI_LONG:
    MRIIseq_vox( mri, x, y, z, 0 ) = (int)val;
    break;
  case MRI_FLOAT:
    MRIIseq_vox( mri, x, y, z, 0 ) = (int)val;
    break;
  case MRI_SHORT:
    MRIIseq_vox( mri, x, y, z, 0 ) = (int)val;
    break;
  default:
    break;
  }
}

int main(int argc, const char **argv)
{
  if (argc < 5) {
    std::cout << "pointset2label <waypoint file> <input volume file> <label value> <output volume label file> [-clear]\n\n"
                 "  Example: pointset2label wp.label input.mgz 3 label_out.mgz\n\n"
                 "  If -clear option is given, input volume will be cleared first. \n"
                 "  If output is the same as input, input volume file will be overwritten.\n\n";
    return -1;
  }

  char* out_fn = (char*)(argv[4]);

  bool bClear = false;
  if (argc > 5)
  {
    bClear = (strcasecmp(argv[5], "-clear") == 0);
  }

  MRI* mri;
  if (bClear)
  {
    MRI* mri_temp = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_temp)
    {
      std::cerr << "Could not read volume info from " << argv[2] << "\n";
      return -1;
    }

    mri = MRIallocSequence( mri_temp->width,
                                mri_temp->height,
                                mri_temp->depth,
                                MRI_INT, 1);
    MRIcopyHeader( mri_temp, mri );
    MRIfree(&mri_temp);
  }
  else
  {
    mri = MRIread(argv[2]);
  }

  LABEL* label = LabelRead(NULL, argv[1]);
  if (!label)
  {
    std::cerr << "Could not load point set from " << argv[1] << "\n";
    return -1;
  }
  if (label->n_points < 2)
  {
    std::cerr << "Label has less than 2 points. Cannot go on.\n";
    return -1;
  }

  int nval = atoi(argv[3]);
  if (nval <= 0)
  {
    std::cerr << "Please enter a valid label value\n";
    return -1;
  }

  if (label->coords != LABEL_COORDS_SCANNER_RAS )
  {
    LabelToScannerRAS(label, mri, label);
    std::cout << "Label coordinates are converted to scanner ras for " << argv[1] << std::endl;
  }

  float* pts = new float[label->n_points*3];
  for (int i = 0; i < label->n_points; i++)
  {
    pts[i*3] = label->lv[i].x;
    pts[i*3+1] = label->lv[i].y;
    pts[i*3+2] = label->lv[i].z;
  }

  float step_length = fmin(mri->xsize, fmin(mri->ysize, mri->zsize));
  Spline spline(pts, label->n_points, 3);
  int num = spline.GetCurveCount(step_length);
  float* new_pts = new float[num*3];
  num = spline.GetCurve(new_pts, step_length);
  double x0, y0, z0, x1, y1, z1;
  for (int i = 0; i < num-1; i++)
  {
    MRIworldToVoxel(mri, new_pts[i*3], new_pts[i*3+1], new_pts[i*3+2], &x0, &y0, &z0);
    MRIworldToVoxel(mri, new_pts[(i+1)*3], new_pts[(i+1)*3+1], new_pts[(i+1)*3+2], &x1, &y1, &z1);
    set_voxel_value(mri, x0, y0, z0, nval);
    set_voxel_value(mri, x1, y1, z1, nval);
    if (fabs(rint(x0)-rint(x1)) > 1 || fabs(rint(y0)-rint(y1)) > 1 || fabs(rint(z0)-rint(z1)) > 1)
      set_voxel_value(mri, (x0+x1)/2, (y0+y1)/2, (z0+z1)/2, nval);
  }

  int err = MRIwrite( mri, out_fn );
  if (err != 0)
    std::cerr << "Failed to write voluem to " << out_fn << std::endl;

  delete[] pts;
  delete[] new_pts;
  LabelFree(&label);
  MRIfree(&mri);
  return 0;
}

