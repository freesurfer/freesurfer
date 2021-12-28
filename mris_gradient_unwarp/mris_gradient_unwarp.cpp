#include <stdio.h>
#include <stdlib.h>

#include "mri.h"
#include "mriBSpline.h"
#include "GradUnwarp.h"

/* examples:
 * 1. do gradient unwarp
 * mris_gradient_unwarp/mris_gradient_unwarp $WKDIR/fs_test/gradunwarp/example/coeff_Sonata.grad $WKDIR/fs_test/gradunwarp/example/orig.mgz unwarped.mgz
 * 2. duplicate orig.mgz without gradient unwarping
 * mris_gradient_unwarp/mris_gradient_unwarp $WKDIR/fs_test/gradunwarp/example/coeff_Sonata.grad $WKDIR/fs_test/gradunwarp/example/orig.mgz orig.5.mgz 1
 */
/****** BIG PICTURES ******/
/*
 * MATRIX *vox2ras = MRIxfmCRS2XYZ(const MRI *mri, 0)
 *
 * vol = MRIread('vol.mgz');
 * [Ax Ay Az Bx By Bz] = ReadYourCoef('a.coef');
 * for col = 0; col < vol->width; col++
 *   for row
 *     for slice
 *       RAS = vox2ras*[C R S 1]'
 *       DeltaRAS = YourFunction(Ax,Ay,Az,Bx,By,Bz,RAS,R0);
 *       DistortedRAS = RAS + DeltaRAS;
 *       DistortedCRS = inv(vox2ras)*DistortedRAS; // floating point
 *       NewVol(col,row,slice) =
 *           MRIsampleVol(vol,distcol,distrow,distslice, interpolationcode);
 *
 * surf = MRISread('lh.white');
 * for n = 0:surf->nvertices-1
 *  VERTEX *v = surf->vertices[n]
 *  v->x, v->y, v->z
 */
int main(int argc, const char *argv[])
{
#if 0 // test GradUnwarp class

  GradUnwarp *gradUnwarp = new GradUnwarp();
  gradUnwarp->read_siemens_coeff(argv[1]);
  //mrisGradientNonlin->printCoeff();

  gradUnwarp->initSiemensLegendreNormfact();

  double x = atof(argv[2]);
  double y = atof(argv[3]);
  double z = atof(argv[4]);

  float Dx, Dy, Dz;
  gradUnwarp->spharm_evaluate(x, y, z, &Dx, &Dy, &Dz);

  printf(" x = %.6lf,  Dx = %.6lf\n", x, Dx);
  printf(" y = %.6lf,  Dy = %.6lf\n", y, Dy);
  printf(" z = %.6lf,  Dz = %.6lf\n", z, Dz);

  delete gradUnwarp;

  return 0;

#else

  // ??? these two variables to be set by user as command line option ???
  int InterpCode = SAMPLE_NEAREST;
  int sinchw = nint(0); 

  int (*nintfunc)( double );
  nintfunc = &nint;

  const char *gradfile = argv[1];
  const char *origmgz  = argv[2];
  const char *unwarpedmgz = argv[3];

  int duporigvol = 0;
  if (argc > 4)
    duporigvol = atoi(argv[4]);
  if (duporigvol)
    printf("\n!!! Duplicate original volume as %s ...\n\n", unwarpedmgz); 

  GradUnwarp *gradUnwarp = new GradUnwarp();
  gradUnwarp->read_siemens_coeff(gradfile);
  gradUnwarp->initSiemensLegendreNormfact();

  MRI *origvol = MRIread(origmgz);
  MATRIX *vox2ras_orig = MRIxfmCRS2XYZ(origvol, 0);

  MRI *unwarpedvol = MRIalloc(origvol->width, origvol->height, origvol->depth, origvol->type);
  unwarpedvol->x_r = origvol->x_r;
  unwarpedvol->x_a = origvol->x_a;
  unwarpedvol->x_s = origvol->x_s;

  unwarpedvol->y_r = origvol->y_r;
  unwarpedvol->y_a = origvol->y_a;
  unwarpedvol->y_s = origvol->y_s;

  unwarpedvol->z_r = origvol->z_r;
  unwarpedvol->z_a = origvol->z_a;
  unwarpedvol->z_s = origvol->z_s;

  unwarpedvol->c_r = origvol->c_r;
  unwarpedvol->c_a = origvol->c_a;
  unwarpedvol->c_s = origvol->c_s;

  unwarpedvol->xsize = origvol->xsize;
  unwarpedvol->ysize = origvol->ysize;
  unwarpedvol->zsize = origvol->zsize;

  MATRIX *CRS = MatrixAlloc(4, 1, MATRIX_REAL);
  MATRIX *RAS = MatrixAlloc(4, 1, MATRIX_REAL);;
  MATRIX *DeltaRAS = MatrixAlloc(4, 1, MATRIX_REAL);
  MATRIX *DistortedRAS = MatrixAlloc(4, 1, MATRIX_REAL);
  MATRIX *DistortedCRS = MatrixAlloc(4, 1, MATRIX_REAL);

  MRI_BSPLINE *bspline = NULL;
  if (InterpCode == SAMPLE_CUBIC_BSPLINE)
    bspline = MRItoBSpline(origvol, NULL, 3);

#if USE_VALVECTS
  float *valvects[_MAX_FS_THREADS];
  valvects[0] = (float *)calloc(sizeof(float), origvol->nframes);
#endif

  printf("width=%4d height=%4d depth=%4d nframes=%4d\n", origvol->width, origvol->height, origvol->depth, origvol->nframes);

  int c, r, s, f;
  for (c = 0; c < origvol->width; c++)
  {
    for (r = 0; r < origvol->height; r++)
    {
      for (s = 0; s < origvol->depth; s++)
      {
        if (!duporigvol)
        {
          // clear CRS, RAS, DeltaRAS, DistortedRAS, DistortedCRS
	  MatrixClear(CRS);
          MatrixClear(RAS);
          MatrixClear(DeltaRAS);
          MatrixClear(DistortedRAS);
          MatrixClear(DistortedCRS);

          CRS->rptr[1][1] = c;
          CRS->rptr[2][1] = r;
          CRS->rptr[3][1] = s;
          CRS->rptr[4][1] = 1;

          // Convert the CRS to RAS
          RAS->rptr[4][1] = 1;
          RAS = MatrixMultiply(vox2ras_orig, CRS, RAS);

          float x = RAS->rptr[1][1];
          float y = RAS->rptr[2][1];
          float z = RAS->rptr[3][1];
          float Dx, Dy, Dz;
          gradUnwarp->spharm_evaluate(x, y, z, &Dx, &Dy, &Dz);
          printf(" x=%.6lf  y=%.6lf  z=%.6lf\n", x, y, z);
          printf("Dx=%.6lf Dy=%.6lf Dz=%.6lf\n", Dx, Dy, Dz);

          DeltaRAS->rptr[1][1] = Dx;
          DeltaRAS->rptr[2][1] = Dy;
          DeltaRAS->rptr[3][1] = Dz;
        
          DistortedRAS = MatrixAdd(RAS, DeltaRAS, DistortedRAS);
          DistortedCRS = MatrixMultiply(MatrixInverse(vox2ras_orig, NULL), DistortedRAS, DistortedCRS);
        } 
 
        float fcs, frs, fss, *valvect;
        int ics, irs, iss;
        double rval = 0;

        if (!duporigvol)
	{
          fcs = DistortedCRS->rptr[1][1];
          frs = DistortedCRS->rptr[2][1];
          fss = DistortedCRS->rptr[3][1];

          ics =  nintfunc(fcs);
          irs =  nintfunc(frs);
          iss =  nintfunc(fss);

          printf("fcs = %lf (%d), frs = %lf (%d), fss = %lf (%d)\n", fcs, ics, frs, irs, fss, iss);
	}
        else
	{
          // this should output the same vol as orig
          ics =  c;
          irs =  r;
          iss =  s;
	}

#if USE_VALVECTS
        valvect = valvects[0];
#else
        valvect = new float[origvol->nframes]; 
#endif

        /* Assign output volume values */
        if (InterpCode == SAMPLE_TRILINEAR)
          MRIsampleSeqVolume(origvol, fcs, frs, fss, valvect, 0, origvol->nframes - 1);
        else {
          for (f = 0; f < origvol->nframes; f++) {
            switch (InterpCode) {
              case SAMPLE_NEAREST:
                //valvect[f] = MRIgetVoxVal(origvol, ics, irs, iss, f);
                rval = MRIgetVoxVal2(origvol, ics, irs, iss, f);
                break;
              case SAMPLE_CUBIC_BSPLINE:
                MRIsampleBSpline(bspline, fcs, frs, fss, f, &rval);
                //valvect[f] = rval;
                break;
              case SAMPLE_SINC: /* no multi-frame */
                MRIsincSampleVolume(origvol, fcs, frs, fss, sinchw, &rval);
                //valvect[f] = rval;
                break;
              default:
                printf("ERROR: MR: interpolation method %i unknown\n", InterpCode);
                exit(1);
            } // switch

            MRIsetVoxVal2(unwarpedvol, c, r, s, f, rval);
            //printf("c=%4d r=%4d s=%4d f=%4d\n", c, r, s, f);
          }
        } // f

      }   // s
    }     // r
  }       // c

  MRIwrite(unwarpedvol, unwarpedmgz);

  delete gradUnwarp;

  MatrixClear(CRS);
  MatrixClear(RAS);
  MatrixClear(DeltaRAS);
  MatrixClear(DistortedRAS);
  MatrixClear(DistortedCRS);

  MRIfree(&origvol);
  MRIfree(&unwarpedvol);

#if USE_VALVECTS
  free(valvects[0]);
#endif
  if (bspline)
    MRIfreeBSpline(&bspline);

  return 0;
#endif
}

