#include <iostream>
#include "Spline.h"
#include "mri.h"
#include "label.h"

int main(int argc, const char **argv)
{
  if (argc < 5)
  {
    std::cout << "pointset2label <waypoint file> <volume label file> <label value> [output volume file]\n\n"
                 "Output file can be the same as the input volume file. In that case input volume file will be updated\n\n"
                 "  Example: pointset2label wp.label foo.mgz 3 foo_out.mgz\n\n";
    return -1;
  }

  char* out_fn = (char*)argv[2];
  if (argc >= 5)
    out_fn = (char*)(argv[4]);

  MRI* mri = MRIread(argv[2]);
  if (!mri)
  {
    std::cerr << "Could not load volume from " << argv[2] << "\n";
    return -1;
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

  float step_length = fmin(mri->xsize, fmin(mri->ysize, mri->zsize))/2;
  Spline spline(pts, label->n_points, 3);
  int num = spline.GetCurveCount(step_length);
  float* new_pts = new float[num*3];
  num = spline.GetCurve(new_pts, step_length);
  int x0, y0, z0, x1, y1, z1;
  for (int i = 0; i < num-1; i++)
  {
    MRIworldToVoxelIndex(mri, new_pts[i*3], new_pts[i*3+1], new_pts[i*3+2], &x0, &y0, &z0);
    MRIworldToVoxelIndex(mri, new_pts[(i+1)*3], new_pts[(i+1)*3+1], new_pts[(i+1)*3+2], &x1, &y1, &z1);
    switch ( mri->type )
    {
    case MRI_UCHAR:
      MRIseq_vox( mri, x0, y0, z0, 0 ) = (unsigned char)nval;
      MRIseq_vox( mri, x1, y1, z1, 0 ) = (unsigned char)nval;
      break;
    case MRI_INT:
      MRIIseq_vox( mri, x0, y0, z0, 0 ) = (int)nval;
      MRIIseq_vox( mri, x1, y1, z1, 0 ) = (int)nval;
      break;
    case MRI_LONG:
      MRIIseq_vox( mri, x0, y0, z0, 0 ) = (int)nval;
      MRIIseq_vox( mri, x1, y1, z1, 0 ) = (int)nval;
      break;
    case MRI_FLOAT:
      MRIIseq_vox( mri, x0, y0, z0, 0 ) = (int)nval;
      MRIIseq_vox( mri, x1, y1, z1, 0 ) = (int)nval;
      break;
    case MRI_SHORT:
      MRIIseq_vox( mri, x0, y0, z0, 0 ) = (int)nval;
      MRIIseq_vox( mri, x1, y1, z1, 0 ) = (int)nval;
      break;
    default:
      break;
    }
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

