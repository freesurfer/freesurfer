/**
 * @file  unwarpGradientNonlinearity.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:55 $
 *    $Revision: 1.6 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */


/*

   Notes on conventions followed here.
   0) CRS means (column,row,slice)
   1) voxel row, column, slice indices (RR, CC, SS) are 0-based
   2) matrices are 1-based
   3) RAS coords versus Siemens XYZ coords. (flip the sign on x and z)

   (gdb) run --unwarp_gradient_nonlinearity allegra fullUnwarp noJacobianCorrection linear 0 testImas/621-4-100.ima COR

*/

#include <stdio.h>
#include <stdlib.h>
#include "mri.h"
#include "matrix.h"
#include "bfileio.h"

#define unwarpGradientNonlinearity_SRC
#include "unwarpGradientNonlinearity.h"

MRI *unwarpGradientNonlinearity(MRI *mri,
                                char *unwarp_gradientType,
                                char *unwarp_partialUnwarp,
                                char *unwarp_jacobianCorrection,
                                char *unwarp_interpType,
                                int unwarp_sincInterpHW)
{
  MRI *mri_unwarped;
  int RR, CC, SS, FF;
  MATRIX *M_tmp, *M_CRS_2_XYZ, *M_XYZ_2_CRS;
  MATRIX *M_XYZ_2_beadIJK, *M_CRS_2_beadIJK;
  MATRIX *voxelCRS, *tmpVec1;
  float beadI, beadJ, beadK;
  float voxel_X, voxel_Y, voxel_Z;
  float voxel_dX, voxel_dY, voxel_dZ;
  float voxel_CC_plus_dCC, voxel_RR_plus_dRR, voxel_SS_plus_dSS;
  float *bead_dX, *bead_dY, *bead_dZ;
  int maxBeadI, maxBeadJ, maxBeadK;
  double tmpVal;
  int index;
  float *voxelCornerCC_slab1, *voxelCornerCC_slab2;
  float *voxelCornerRR_slab1, *voxelCornerRR_slab2;
  float *voxelCornerSS_slab1, *voxelCornerSS_slab2;
  float *tmpPtr;
  int tt, nA, nB, nC;
  float cubeVol;
  float cv[8][3];
  int cubeFaces[12][3] = {
                           {1,     2,     3},
                           {1,     3,     4},
                           {1,     4,     5},
                           {5,     4,     8},
                           {6,     5,     8},
                           {7,     6,     8},
                           {2,     6,     7},
                           {3,     2,     7},
                           {2,     1,     5},
                           {2,     5,     6},
                           {4,     3,     8},
                           {3,     7,     8}
                         };

  /* Allocate memory and provide header for mri_unwarped */
  mri_unwarped = MRIallocSequence(mri->width,
                                  mri->height,
                                  mri->depth,
                                  MRI_FLOAT,
                                  mri->nframes);
  MRIcopyHeader(mri, mri_unwarped);

  /* Read in scanner specific distortion data */
  M_XYZ_2_beadIJK = MatrixAlloc(4, 4, MATRIX_REAL);
  uGN_loadGradientData(unwarp_gradientType,
                       M_XYZ_2_beadIJK,
                       &bead_dX,
                       &bead_dY,
                       &bead_dZ,
                       &maxBeadI,
                       &maxBeadJ,
                       &maxBeadK);

  /* Read in the matrix that maps CRS {column, row, slice} to Siemens XYZ
     scanner coordinates. Note that this involves putting minus signs
     in front of the first and third rows of the matrix that maps
     {column, row, slice} to RAS coordinates */
  M_tmp = MatrixIdentity(4, NULL);
  *(MATRIX_RELT(M_tmp, 1, 1)) = -1.0;
  *(MATRIX_RELT(M_tmp, 3, 3)) = -1.0;

  M_CRS_2_XYZ = MRIgetVoxelToRasXform(mri);
  MatrixMultiply(M_tmp, M_CRS_2_XYZ, M_CRS_2_XYZ);

  /* Invert for the matrix that takes you from Siemens scanner XYZ
     coordinates to CRS indices/coords */
  M_XYZ_2_CRS = MatrixInverse(M_CRS_2_XYZ, NULL);

  /* Compute the matrix that takes CRS to beadIJK */
  M_CRS_2_beadIJK = MatrixMultiply(M_XYZ_2_beadIJK, M_CRS_2_XYZ, NULL);

  /* voxelCRS is a matrix struct for the vector [C;R;S;1] */
  voxelCRS = VectorAlloc(4, MATRIX_REAL);
  *(MATRIX_RELT(voxelCRS, 1, 1)) = 0;
  *(MATRIX_RELT(voxelCRS, 2, 1)) = 0;
  *(MATRIX_RELT(voxelCRS, 3, 1)) = 0;
  *(MATRIX_RELT(voxelCRS, 4, 1)) = 1;

  /* tmpVec1 is used to store intermediate values */
  tmpVec1 = VectorAlloc(4, MATRIX_REAL);

  /* loop through each voxel */
  for (SS = 0; SS < mri->depth; SS++)
  {
    *(MATRIX_RELT(voxelCRS, 3, 1)) = SS;
    for (RR = 0; RR < mri->height; RR++)
    {
      *(MATRIX_RELT(voxelCRS, 2, 1)) = RR;
      for (CC = 0; CC < mri->width; CC++)
      {
        *(MATRIX_RELT(voxelCRS, 1, 1)) = CC;

        /* compute Siemens scanner XYZ for this voxel */
        MatrixMultiply(M_CRS_2_XYZ, voxelCRS, tmpVec1);
        voxel_X = *(MATRIX_RELT(tmpVec1, 1, 1));
        voxel_Y = *(MATRIX_RELT(tmpVec1, 2, 1));
        voxel_Z = *(MATRIX_RELT(tmpVec1, 3, 1));

        /* compute beadIJK for this voxel */
        MatrixMultiply(M_CRS_2_beadIJK, voxelCRS, tmpVec1);
        beadI = *(MATRIX_RELT(tmpVec1, 1, 1));
        beadJ = *(MATRIX_RELT(tmpVec1, 2, 1));
        beadK = *(MATRIX_RELT(tmpVec1, 3, 1));

        /* linearly interpolate to find [dX dY dZ] for this voxel */
        uGN_linInterp(bead_dX, bead_dY, bead_dZ,
                      beadI, beadJ, beadK,
                      maxBeadI, maxBeadJ, maxBeadK,
                      &voxel_dX, &voxel_dY, &voxel_dZ);

        /* map X+dX, etc. to CC+dCC, etc. */
        *(MATRIX_RELT(tmpVec1, 1, 1)) = voxel_X + voxel_dX;
        *(MATRIX_RELT(tmpVec1, 2, 1)) = voxel_Y + voxel_dY;
        *(MATRIX_RELT(tmpVec1, 3, 1)) = voxel_Z + voxel_dZ;
        MatrixMultiply(M_XYZ_2_CRS, tmpVec1, tmpVec1);
        voxel_CC_plus_dCC = *(MATRIX_RELT(tmpVec1, 1, 1));
        voxel_RR_plus_dRR = *(MATRIX_RELT(tmpVec1, 2, 1));
        voxel_SS_plus_dSS = *(MATRIX_RELT(tmpVec1, 3, 1));

        /* interpolate to find the value at CC+dCC, etc. */
        for (FF = 0; FF < mri->nframes; FF++)
        {
          if (strcmp(unwarp_interpType, "linear") == 0)
            MRIsampleVolumeFrame(mri,
                                 ((double) voxel_CC_plus_dCC),
                                 ((double) voxel_RR_plus_dRR),
                                 ((double) voxel_SS_plus_dSS),
                                 FF,
                                 &tmpVal);
          else if (strcmp(unwarp_interpType, "sinc") == 0)
          {
            /* Note: after twitzel writes
               MRIsincSampleVolumeFrame, uncomment the line
               below, and cut out the lines between %^& */

            /*
            MRIsincSampleVolumeFrame(mri,
            ((double) voxel_CC_plus_dCC),
            ((double) voxel_RR_plus_dRR),
            ((double) voxel_SS_plus_dSS),
            FF,
            unwarp_sincInterpHW,
            &tmpVal);
            */

            /* %^& -- start cut here --- */
            if (mri->nframes != 1)
            {
              printf("\nAt present, sinc interpolation ");
              printf("not supported for 4-D data.\n");
              printf("Use linear interpolation instead\n\n");
              exit(1);
            }
            else
              MRIsincSampleVolume(mri,
                                  ((double) voxel_CC_plus_dCC),
                                  ((double) voxel_RR_plus_dRR),
                                  ((double) voxel_SS_plus_dSS),
                                  unwarp_sincInterpHW,
                                  &tmpVal);
            /* %^& -- end cut here --- */
          }

          /* write interpolated value into mri_unwarped */
          (MRIFseq_vox(mri_unwarped, CC, RR, SS, FF)) = tmpVal;

        } /* FF */

      } /* CC */

    } /* RR */

  } /* SS */

  /* if Jacobian correction */
  if (strcmp(unwarp_jacobianCorrection, "JacobianCorrection") == 0)
  {

    voxelCornerCC_slab1 = malloc(sizeof(float)*((mri->width)+1)*((mri->height)+1));
    voxelCornerCC_slab2 = malloc(sizeof(float)*((mri->width)+1)*((mri->height)+1));
    voxelCornerRR_slab1 = malloc(sizeof(float)*((mri->width)+1)*((mri->height)+1));
    voxelCornerRR_slab2 = malloc(sizeof(float)*((mri->width)+1)*((mri->height)+1));
    voxelCornerSS_slab1 = malloc(sizeof(float)*((mri->width)+1)*((mri->height)+1));
    voxelCornerSS_slab2 = malloc(sizeof(float)*((mri->width)+1)*((mri->height)+1));

    *(MATRIX_RELT(voxelCRS, 3, 1)) = -0.5;
    for (RR = 0; RR <= mri->height; RR++)
    {
      *(MATRIX_RELT(voxelCRS, 2, 1)) = RR - 0.5;
      for (CC = 0; CC <= mri->width; CC++)
      {
        *(MATRIX_RELT(voxelCRS, 1, 1)) = CC - 0.5;

        /* index into the voxelCorner slabs */
        index = CC + RR*((mri->width) + 1);

        /* compute Siemens scanner XYZ for this point */
        MatrixMultiply(M_CRS_2_XYZ, voxelCRS, tmpVec1);
        voxel_X = *(MATRIX_RELT(tmpVec1, 1, 1));
        voxel_Y = *(MATRIX_RELT(tmpVec1, 2, 1));
        voxel_Z = *(MATRIX_RELT(tmpVec1, 3, 1));

        /* compute beadIJK for this point */
        MatrixMultiply(M_CRS_2_beadIJK, voxelCRS, tmpVec1);
        beadI = *(MATRIX_RELT(tmpVec1, 1, 1));
        beadJ = *(MATRIX_RELT(tmpVec1, 2, 1));
        beadK = *(MATRIX_RELT(tmpVec1, 3, 1));

        /* linearly interpolate to find [dX dY dZ] for this point */
        uGN_linInterp(bead_dX, bead_dY, bead_dZ,
                      beadI, beadJ, beadK,
                      maxBeadI, maxBeadJ, maxBeadK,
                      &voxel_dX, &voxel_dY, &voxel_dZ);

        /* map X+dX, etc. to CC+dCC, etc. */
        *(MATRIX_RELT(tmpVec1, 1, 1)) = voxel_X + voxel_dX;
        *(MATRIX_RELT(tmpVec1, 2, 1)) = voxel_Y + voxel_dY;
        *(MATRIX_RELT(tmpVec1, 3, 1)) = voxel_Z + voxel_dZ;
        MatrixMultiply(M_XYZ_2_CRS, tmpVec1, tmpVec1);
        voxelCornerCC_slab2[index] = *(MATRIX_RELT(tmpVec1, 1, 1));
        voxelCornerRR_slab2[index] = *(MATRIX_RELT(tmpVec1, 2, 1));
        voxelCornerSS_slab2[index] = *(MATRIX_RELT(tmpVec1, 3, 1));

      } /* CC */

    } /* RR */

    /* loop through each voxel */
    for (SS = 0; SS < mri->depth; SS++)
    {
      *(MATRIX_RELT(voxelCRS, 3, 1)) = SS + 0.5;

      /* swap slab1 and slab2 */
      tmpPtr = voxelCornerCC_slab1;
      voxelCornerCC_slab1 = voxelCornerCC_slab2;
      voxelCornerCC_slab2 = tmpPtr;

      tmpPtr = voxelCornerRR_slab1;
      voxelCornerRR_slab1 = voxelCornerRR_slab2;
      voxelCornerRR_slab2 = tmpPtr;

      tmpPtr = voxelCornerSS_slab1;
      voxelCornerSS_slab1 = voxelCornerSS_slab2;
      voxelCornerSS_slab2 = tmpPtr;

      for (RR = 0; RR <= mri->height; RR++)
      {
        *(MATRIX_RELT(voxelCRS, 2, 1)) = RR - 0.5;
        for (CC = 0; CC <= mri->width; CC++)
        {
          *(MATRIX_RELT(voxelCRS, 1, 1)) = CC - 0.5;

          /* index into the voxelCorner slabs */
          index = CC + RR*((mri->width) + 1);

          /* compute Siemens scanner XYZ for this voxel */
          MatrixMultiply(M_CRS_2_XYZ, voxelCRS, tmpVec1);
          voxel_X = *(MATRIX_RELT(tmpVec1, 1, 1));
          voxel_Y = *(MATRIX_RELT(tmpVec1, 2, 1));
          voxel_Z = *(MATRIX_RELT(tmpVec1, 3, 1));

          /* compute beadIJK for this voxel */
          MatrixMultiply(M_CRS_2_beadIJK, voxelCRS, tmpVec1);
          beadI = *(MATRIX_RELT(tmpVec1, 1, 1));
          beadJ = *(MATRIX_RELT(tmpVec1, 2, 1));
          beadK = *(MATRIX_RELT(tmpVec1, 3, 1));

          /* linearly interpolate to find [dX dY dZ] for this voxel */
          uGN_linInterp(bead_dX, bead_dY, bead_dZ,
                        beadI, beadJ, beadK,
                        maxBeadI, maxBeadJ, maxBeadK,
                        &voxel_dX, &voxel_dY, &voxel_dZ);

          /* map X+dX, etc. to CC+dCC, etc. */
          *(MATRIX_RELT(tmpVec1, 1, 1)) = voxel_X + voxel_dX;
          *(MATRIX_RELT(tmpVec1, 2, 1)) = voxel_Y + voxel_dY;
          *(MATRIX_RELT(tmpVec1, 3, 1)) = voxel_Z + voxel_dZ;
          MatrixMultiply(M_XYZ_2_CRS, tmpVec1, tmpVec1);
          voxelCornerCC_slab2[index] = *(MATRIX_RELT(tmpVec1, 1, 1));
          voxelCornerRR_slab2[index] = *(MATRIX_RELT(tmpVec1, 2, 1));
          voxelCornerSS_slab2[index] = *(MATRIX_RELT(tmpVec1, 3, 1));

        } /* CC */

      } /* RR */

      for (RR = 0; RR < mri->height; RR++)
      {
        for (CC = 0; CC < mri->width; CC++)
        {
          /* index into the voxelCorner slabs */
          index = CC + RR*((mri->width) + 1);

          cv[0][0] = voxelCornerRR_slab1[index];
          cv[0][1] = voxelCornerCC_slab1[index];
          cv[0][2] = voxelCornerSS_slab1[index];

          cv[4][0] = voxelCornerRR_slab2[index];
          cv[4][1] = voxelCornerCC_slab2[index];
          cv[4][2] = voxelCornerSS_slab2[index];

          /* index into the voxelCorner slabs */
          index = CC + (RR+1)*((mri->width) + 1);

          cv[1][0] = voxelCornerRR_slab1[index];
          cv[1][1] = voxelCornerCC_slab1[index];
          cv[1][2] = voxelCornerSS_slab1[index];

          cv[5][0] = voxelCornerRR_slab2[index];
          cv[5][1] = voxelCornerCC_slab2[index];
          cv[5][2] = voxelCornerSS_slab2[index];

          /* index into the voxelCorner slabs */
          index = (CC+1) + (RR+1)*((mri->width) + 1);

          cv[2][0] = voxelCornerRR_slab1[index];
          cv[2][1] = voxelCornerCC_slab1[index];
          cv[2][2] = voxelCornerSS_slab1[index];

          cv[6][0] = voxelCornerRR_slab2[index];
          cv[6][1] = voxelCornerCC_slab2[index];
          cv[6][2] = voxelCornerSS_slab2[index];

          /* index into the voxelCorner slabs */
          index = (CC+1) + RR*((mri->width) + 1);

          cv[3][0] = voxelCornerRR_slab1[index];
          cv[3][1] = voxelCornerCC_slab1[index];
          cv[3][2] = voxelCornerSS_slab1[index];

          cv[7][0] = voxelCornerRR_slab2[index];
          cv[7][1] = voxelCornerCC_slab2[index];
          cv[7][2] = voxelCornerSS_slab2[index];

          cubeVol = 0;
          for (tt = 0; tt < 12; tt++)
          {
            nA = cubeFaces[tt][0] - 1;
            nB = cubeFaces[tt][1] - 1;
            nC = cubeFaces[tt][2] - 1;

            cubeVol +=
              (+ cv[nA][0]*(cv[nB][1]*cv[nC][2] - cv[nB][2]*cv[nC][1])
               - cv[nA][1]*(cv[nB][0]*cv[nC][2] - cv[nB][2]*cv[nC][0])
               + cv[nA][2]*(cv[nB][0]*cv[nC][1] - cv[nB][1]*cv[nC][0]));
          }
          cubeVol /= 6;

          /* The orientation of cubeFaces above was
          (unfortunately) chosen to be one that results in
          negative volume for the cube. We could reorder
          cubeFaces, but instead, we will be lazy and simply
          invert the sign of cubeVol */
          cubeVol = -cubeVol;

          /* Negative Jacobians are set to 0 */
          if (cubeVol < 0)
            cubeVol = 0;

          /* If the Jacobian is too large, send it to 0 gracefully */
          if (cubeVol > 10)
            cubeVol = 20 - cubeVol;

          /* The line above will make any Jacobians greater
                         than 20 negative, so set these to zero */
          if (cubeVol < 0)
            cubeVol = 0;

          /* write jacobian corrected value into mri_unwarped */
          for (FF = 0; FF < mri->nframes; FF++)
          {
            (MRIFseq_vox(mri_unwarped, CC, RR, SS, FF)) *= cubeVol;
          } /* FF */

        } /* CC */

      } /* RR */

    } /* SS */

    free(voxelCornerCC_slab1);
    free(voxelCornerCC_slab2);
    free(voxelCornerRR_slab1);
    free(voxelCornerRR_slab2);
    free(voxelCornerSS_slab1);
    free(voxelCornerSS_slab2);

  } /* Jacobian correction */

  free(bead_dX);
  free(bead_dY);
  free(bead_dZ);

  MatrixFree(&M_tmp);
  MatrixFree(&M_CRS_2_XYZ);
  MatrixFree(&M_XYZ_2_CRS);
  MatrixFree(&M_XYZ_2_beadIJK);
  MatrixFree(&M_CRS_2_beadIJK);
  MatrixFree(&voxelCRS);
  MatrixFree(&tmpVec1);

  return(mri_unwarped);

}

int uGN_loadGradientData(char *unwarp_gradientType,
                         MATRIX *M_XYZ_2_beadIJK,
                         float **p_bead_dX,
                         float **p_bead_dY,
                         float **p_bead_dZ,
                         int *p_maxBeadI,
                         int *p_maxBeadJ,
                         int *p_maxBeadK)
{
  FILE *fp;
  float tmpMat[16];
  int numEls;
  char *FREESURFER_HOME;
  char dataFile[STRLEN];
  float tmpF;
  int tmpI;
  int BYTESWAP = 0;
  char fileDescriptor[STRLEN];

  FREESURFER_HOME = getenv("FREESURFER_HOME");
  if (!FREESURFER_HOME)
  {
    printf("Error: FREESURFER_HOME not defined\n");
    exit(1);
  }

  strcpy(dataFile, FREESURFER_HOME);

  /* Determine gradient type: sonata or allegra */
  if (strcmp(unwarp_gradientType, "sonata")  == 0)
  {
    strcat(dataFile, "/data/gradientNonlinearitySonata.mff");
    fp = fopen(dataFile, "rb");
  }
  else if (strcmp(unwarp_gradientType, "allegra") == 0)
  {
    strcat(dataFile, "/data/gradientNonlinearityAllegra.mff");
    fp = fopen(dataFile, "rb");
  }
  else if (strcmp(unwarp_gradientType, "GE") == 0)
  {
    strcat(dataFile, "/data/gradientNonlinearityGE.mff");
    fp = fopen(dataFile, "rb");
  }
  else
  {
    printf("Error: unknown scanner type\n");
    exit(1);
  }

  if (fp == NULL)
  {
    printf("Error: couldn't open %s\n", dataFile);
    exit(1);
  }

  /* Read in the file description (a string) */
  fscanf(fp, "%s\n", fileDescriptor);

  /* Check that distortion data will be able to be read in OK:
     i.e. need to verify endianness and 4-byte floats and ints */
  fread(&tmpI, sizeof(int), 1, fp);
  fread(&tmpF, sizeof(float), 1, fp);
  if ( (tmpI != 1) || (tmpF != 2) )
  {
    BYTESWAP = 1;

    if ( (sizeof(int) != 4) || (sizeof(float) != 4) )
    {
      printf("Error: because floats or ints are not 4-bytes\n");
      printf("on this architecture, there was a problem reading\n");
      printf("in the distortion offset data in %s/data/\n", FREESURFER_HOME);
      exit(1);
    }
  }

  /*
    printf("%s\n", fileDescriptor);
    printf("%d %f\n", tmpI, tmpF);
    byteswapbuffloat(&tmpI, 4);
    byteswapbuffloat(&tmpF, 4);
    printf("%d %f\n", tmpI, tmpF);
  */

  /* Read in the matrix that takes you from scanner XYZ coordinates to
     bead IJK indices */
  fread(tmpMat, sizeof(float), 16, fp);
  if (BYTESWAP)
    byteswapbuffloat(tmpMat, 4*16);

  *(MATRIX_RELT(M_XYZ_2_beadIJK, 1, 1)) = tmpMat[0];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 1, 2)) = tmpMat[1];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 1, 3)) = tmpMat[2];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 1, 4)) = tmpMat[3];

  *(MATRIX_RELT(M_XYZ_2_beadIJK, 2, 1)) = tmpMat[4];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 2, 2)) = tmpMat[5];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 2, 3)) = tmpMat[6];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 2, 4)) = tmpMat[7];

  *(MATRIX_RELT(M_XYZ_2_beadIJK, 3, 1)) = tmpMat[8];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 3, 2)) = tmpMat[9];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 3, 3)) = tmpMat[10];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 3, 4)) = tmpMat[11];

  *(MATRIX_RELT(M_XYZ_2_beadIJK, 4, 1)) = tmpMat[12];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 4, 2)) = tmpMat[13];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 4, 3)) = tmpMat[14];
  *(MATRIX_RELT(M_XYZ_2_beadIJK, 4, 4)) = tmpMat[15];

  /* Read in the number of bead displacement vectors */
  fread(p_maxBeadI, sizeof(int), 1, fp);
  if (BYTESWAP)
    byteswapbuffloat(p_maxBeadI, 4);
  fread(p_maxBeadJ, sizeof(int), 1, fp);
  if (BYTESWAP)
    byteswapbuffloat(p_maxBeadJ, 4);
  fread(p_maxBeadK, sizeof(int), 1, fp);
  if (BYTESWAP)
    byteswapbuffloat(p_maxBeadK, 4);

  numEls = (*p_maxBeadI)*(*p_maxBeadJ)*(*p_maxBeadK);

  /* Read in the bead displacement information */
  *p_bead_dX = malloc(numEls*sizeof(float));
  fread(*p_bead_dX, sizeof(float), numEls, fp);
  if (BYTESWAP)
    byteswapbuffloat(*p_bead_dX, 4*numEls);

  *p_bead_dY = malloc(numEls*sizeof(float));
  fread(*p_bead_dY, sizeof(float), numEls, fp);
  if (BYTESWAP)
    byteswapbuffloat(*p_bead_dY, 4*numEls);

  *p_bead_dZ = malloc(numEls*sizeof(float));
  fread(*p_bead_dZ, sizeof(float), numEls, fp);
  if (BYTESWAP)
    byteswapbuffloat(*p_bead_dZ, 4*numEls);

  fclose(fp);

  return(1);

}

int uGN_linInterp(float *bead_dX, float *bead_dY, float *bead_dZ,
                  float beadI, float beadJ, float beadK,
                  int maxBeadI, int maxBeadJ, int maxBeadK,
                  float *p_voxel_dX, float *p_voxel_dY, float *p_voxel_dZ)
{

  int index;
  int iii, jjj, kkk;
  int iiip1, jjjp1, kkkp1;
  float XXX[2][2][2], YYY[2][2][2], ZZZ[2][2][2];
  float aaa, bbb, ccc;


  /* Get the indices corresponding to the 8 integer points (i.e. cube
     corners) surrounding our real-valued point (beadI, beadJ, beadK) */
  iii = beadI;
  jjj = beadJ;
  kkk = beadK;
  iiip1 = iii + 1;
  jjjp1 = jjj + 1;
  kkkp1 = kkk + 1;

  /* Check that indices stay in bounds */
  if (iii < 1)
  {
    iii = 1;
    iiip1 = 1;
  }

  if (iii > maxBeadI)
  {
    iii = maxBeadI;
    iiip1 = maxBeadI;
  }

  if (jjj < 1)
  {
    jjj = 1;
    jjjp1 = 1;
  }

  if (jjj > maxBeadJ)
  {
    jjj = maxBeadJ;
    jjjp1 = maxBeadJ;
  }

  if (kkk < 1)
  {
    kkk = 1;
    kkkp1 = 1;
  }

  if (kkk > maxBeadK)
  {
    kkk = maxBeadK;
    kkkp1 = maxBeadK;
  }

  /* Read in the values at the cube corners */
  index = (iii-1) + (jjj-1)*maxBeadI + (kkk-1)*maxBeadI*maxBeadJ;
  XXX[0][0][0] = bead_dX[index];
  YYY[0][0][0] = bead_dY[index];
  ZZZ[0][0][0] = bead_dZ[index];

  index = (iiip1-1) + (jjj-1)*maxBeadI + (kkk-1)*maxBeadI*maxBeadJ;
  XXX[1][0][0] = bead_dX[index];
  YYY[1][0][0] = bead_dY[index];
  ZZZ[1][0][0] = bead_dZ[index];

  index = (iii-1) + (jjjp1-1)*maxBeadI + (kkk-1)*maxBeadI*maxBeadJ;
  XXX[0][1][0] = bead_dX[index];
  YYY[0][1][0] = bead_dY[index];
  ZZZ[0][1][0] = bead_dZ[index];

  index = (iii-1) + (jjj-1)*maxBeadI + (kkkp1-1)*maxBeadI*maxBeadJ;
  XXX[0][0][1] = bead_dX[index];
  YYY[0][0][1] = bead_dY[index];
  ZZZ[0][0][1] = bead_dZ[index];

  index = (iiip1-1) + (jjjp1-1)*maxBeadI + (kkk-1)*maxBeadI*maxBeadJ;
  XXX[1][1][0] = bead_dX[index];
  YYY[1][1][0] = bead_dY[index];
  ZZZ[1][1][0] = bead_dZ[index];

  index = (iiip1-1) + (jjj-1)*maxBeadI + (kkkp1-1)*maxBeadI*maxBeadJ;
  XXX[1][0][1] = bead_dX[index];
  YYY[1][0][1] = bead_dY[index];
  ZZZ[1][0][1] = bead_dZ[index];

  index = (iii-1) + (jjjp1-1)*maxBeadI + (kkkp1-1)*maxBeadI*maxBeadJ;
  XXX[0][1][1] = bead_dX[index];
  YYY[0][1][1] = bead_dY[index];
  ZZZ[0][1][1] = bead_dZ[index];

  index = (iiip1-1) + (jjjp1-1)*maxBeadI + (kkkp1-1)*maxBeadI*maxBeadJ;
  XXX[1][1][1] = bead_dX[index];
  YYY[1][1][1] = bead_dY[index];
  ZZZ[1][1][1] = bead_dZ[index];

  aaa = beadI - ((int) beadI);
  bbb = beadJ - ((int) beadJ);
  ccc = beadK - ((int) beadK);

  *p_voxel_dX = (XXX[0][0][0]*(1-aaa)*(1-bbb)*(1-ccc) +
                 XXX[1][0][0]*aaa*(1-bbb)*(1-ccc) +
                 XXX[0][1][0]*(1-aaa)*bbb*(1-ccc) +
                 XXX[0][0][1]*(1-aaa)*(1-bbb)*ccc +
                 XXX[1][1][0]*aaa*bbb*(1-ccc) +
                 XXX[1][0][1]*aaa*(1-bbb)*ccc +
                 XXX[0][1][1]*(1-aaa)*bbb*ccc +
                 XXX[1][1][1]*aaa*bbb*ccc);

  *p_voxel_dY = (YYY[0][0][0]*(1-aaa)*(1-bbb)*(1-ccc) +
                 YYY[1][0][0]*aaa*(1-bbb)*(1-ccc) +
                 YYY[0][1][0]*(1-aaa)*bbb*(1-ccc) +
                 YYY[0][0][1]*(1-aaa)*(1-bbb)*ccc +
                 YYY[1][1][0]*aaa*bbb*(1-ccc) +
                 YYY[1][0][1]*aaa*(1-bbb)*ccc +
                 YYY[0][1][1]*(1-aaa)*bbb*ccc +
                 YYY[1][1][1]*aaa*bbb*ccc);

  *p_voxel_dZ = (ZZZ[0][0][0]*(1-aaa)*(1-bbb)*(1-ccc) +
                 ZZZ[1][0][0]*aaa*(1-bbb)*(1-ccc) +
                 ZZZ[0][1][0]*(1-aaa)*bbb*(1-ccc) +
                 ZZZ[0][0][1]*(1-aaa)*(1-bbb)*ccc +
                 ZZZ[1][1][0]*aaa*bbb*(1-ccc) +
                 ZZZ[1][0][1]*aaa*(1-bbb)*ccc +
                 ZZZ[0][1][1]*(1-aaa)*bbb*ccc +
                 ZZZ[1][1][1]*aaa*bbb*ccc);



  return(1);

}

/* printf("!!! %f !!!\n", ((double) MRISvox(mri, 114, 114, 114))); */
/* MRIsampleVolume(mri, 114, 114, 114, &tmpmuk); */

