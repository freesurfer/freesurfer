/**
 * @brief DICOM 3.0 reading functions
 *
 * Copyright © 2011-2013 The General Hospital Corporation (Boston, MA) "MGH"
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

#include <ctype.h>
#include <dirent.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/file.h>
#include <sys/time.h>
#include <sys/timeb.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#ifndef Darwin
#include <malloc.h>
#else
void *malloc(size_t size);
#endif

#include <math.h>

#include "mri.h"

#include "diag.h"
#include "dti.h"
#include "fio.h"
#include "fsenv.h"
#include "macros.h"  // DEGREES
#include "mosaic.h"
#include "mri_identify.h"

// #include "affine.h"

#define _DICOMRead_SRC
#include "DICOMRead.h"
#undef _DICOMRead_SRC

static int DCMPrintCond(CONDITION cond);
void *ReadDICOMImage2(int nfiles, DICOMInfo **aDicomInfo, int startIndex);

static BOOL IsTagPresent[NUMBEROFTAGS];
static int sliceDirCosPresent;
static const char *jpegCompressed_UID = "1.2.840.10008.1.2.4";
static const char *rllEncoded_UID = "1.2.840.10008.1.2.5";

//#define _DEBUG

/*-----------------------------------------------------------------
  sdcmLoadVolume() - this loads a volume stored in Siemens DICOM
  format. It scans in the directory of dcmfile searching for all
  Siemens DICOM files. It then determines which of these files belong
  to the same run as dcmfile. If LoadVolume=1, then the pixel data is
  loaded, otherwise, only the header info is assigned.

  Determining which files belong in the same with dcmfile is a
  tricky process. It requires that all files be queried. It also
  requires that all runs be complete. A missing file will cause
  the process to crash even if the file does not belong to the
  same run as dcmfile.

  Notes:
  1. This currently works only for unsigned short data.
  2. It handles mosaics but not supermosaics.
  3. It assumes that slices within a mosaic are sorted in
  anatomical order.
  4. It handles multiple frames for mosaics and non-mosaics.

  SDCMListFile is a file with a list of Siemens DICOM files to
  be unpacked as one run. If using mri_convert, it can be passed
  with the --sdcmlist.
  -----------------------------------------------------------------*/
MRI *sdcmLoadVolume(const char *dcmfile, int LoadVolume, int nthonly)
{
  SDCMFILEINFO *sdfi;
  DCM_ELEMENT *element;
  SDCMFILEINFO **sdfi_list;
  int nthfile;
  int nlist;
  int ncols, nrows, nslices, nframes;
  int nmoscols, nmosrows;
  int mosrow, moscol, mosindex;
  int err, OutOfBounds, IsMosaic;
  int row, col, slice, frame, fid;
  unsigned short *pixeldata;
  MRI *vol, *voltmp;
  char **SeriesList;
  char *tmpstring;
  int Maj, Min, MinMin;
  double xs, ys, zs, xe, ye, ze, d, MinSliceScaleFactor, val;
  int nnlist;
  int vol_datatype;
  FSENV *env;
  std::string tmpfilestdout, cmd, tmpfile, FileNameUse;
  int IsCompressed, IsDWI;
  extern int sliceDirCosPresent;  // set when no ascii header

  xs = ys = zs = xe = ye = ze = d = 0.; /* to avoid compiler warnings */
  slice = 0;
  frame = 0; /* to avoid compiler warnings */

  sliceDirCosPresent = 0;  // assume not present

  /* split progress to 3 parts */
  int nstart = global_progress_range[0];
  int nend = global_progress_range[1];
  global_progress_range[1] = nstart + (nend - nstart) / 3;
  if (SDCMListFile != NULL)
    SeriesList = ReadSiemensSeries(SDCMListFile, &nlist, dcmfile);
  else
    SeriesList = ScanSiemensSeries(dcmfile, &nlist);

  if (SeriesList == NULL) {
    fprintf(stderr, "ERROR: could not find any files (SeriesList==NULL)\n");
    return (NULL);
  }

  if (nlist == 0) {
    fprintf(stderr, "ERROR: could not find any files (nlist==0)\n");
    return (NULL);
  }

  // for(nnlist=0; nnlist<nlist; nnlist++) fprintf(stdout,"%3d  %s\n",
  // nnlist,SeriesList[nnlist]);
  // fflush(stdout);

  global_progress_range[0] = global_progress_range[1];
  global_progress_range[1] += (nend - nstart) / 3;
  printf("INFO: loading series header info.\n");
  sdfi_list = LoadSiemensSeriesInfo(SeriesList, nlist);

  // free memory
  nnlist = nlist;
  while (nnlist--) free(SeriesList[nnlist]);
  free(SeriesList);

  printf("INFO: sorting.\n");
  SortSDCMFileInfo(sdfi_list, nlist);

  sdfiAssignRunNo2(sdfi_list, nlist);

  /* First File in the Run */
  if (nthonly < 0)
    sdfi = sdfi_list[0];
  else
    sdfi = sdfi_list[nthonly];

  /* There are some Siemens files don't have the slice dircos
     because they dont have ascii header (anonymization?).
     If this happens on a mosiac, then you cannot unpack it.
     For non-mosaics, recompute the slice dircos based on
     the image position of the first and last file. Note:
     slice dircos may be used to sort the files in the first
     place, but this should work either way.
  */
  if (sliceDirCosPresent == 0) {
    printf("\n");
    printf(
        "WARNING: file %s does not contain a Siemens ASCII header\n"
        "has this file been anonymized?\n",
        sdfi_list[0]->FileName);
    if (sdfi_list[0]->IsMosaic) {
      printf("ERROR: cannot unpack mosiacs without ASCII header\n");
      exit(1);
    }
    printf("Proceeding as best as I can ... \n");
    printf("\n");
    xs = sdfi_list[0]->ImgPos[0];
    ys = sdfi_list[0]->ImgPos[1];
    zs = sdfi_list[0]->ImgPos[2];
    xe = sdfi_list[nlist - 1]->ImgPos[0];
    ye = sdfi_list[nlist - 1]->ImgPos[1];
    ze = sdfi_list[nlist - 1]->ImgPos[2];
    d = sqrt((xe - xs) * (xe - xs) + (ye - ys) * (ye - ys) + (ze - zs) * (ze - zs));
    sdfi->Vs[0] = (xe - xs) / d;
    sdfi->Vs[1] = (ye - ys) / d;
    sdfi->Vs[2] = (ze - zs) / d;
    // Note: no need to change to RAS because ImagePos is already RAS
  }
  sdfiFixImagePosition(sdfi);
  sdfiVolCenter(sdfi);

  /* for easy access */
  ncols = sdfi->VolDim[0];
  nrows = sdfi->VolDim[1];
  nslices = sdfi->VolDim[2];
  nframes = sdfi->NFrames;
  IsMosaic = sdfi->IsMosaic;

  // verify nthframe is within the range
  if (nthonly >= 0) {
    if (nthonly > nframes - 1) {
      printf("ERROR: only has %d frames (%d - %d) but called for %d\n", nframes, 0, nframes - 1, nthonly);
      fflush(stdout);
      return (NULL);
    }
  }

  printf("INFO: (%3d %3d %3d), nframes = %d, ismosaic=%d\n", ncols, nrows, nslices, nframes, IsMosaic);
  fflush(stdout);

  /*
    6/19/18 Notes on rescaling: there are two possible types of
    rescaling. (1) The autoscale functor has been used.  This is
    indicated by 0x20, 0x4000. (2) The slice scale and intercept
    indicated by 0x28, 0x1052 and 0x1053. The implementation here can
    handle both simultaneously, though hopefully this will never
    occur. I think the autoscale functor was somewhat
    experimental. However, the slope and intercept can be common,
    especially in PET. I have found it in some ASL scans too (probably
    quantitative).
   */

  /* PW 2012/09/06: If the Autoscale functor has been used, save the MRI
     structure as floats, otherwise use shorts  */
  if (sdfi->UseSliceScaleFactor)
    vol_datatype = MRI_FLOAT;
  else
    vol_datatype = MRI_SHORT;
  printf("sdfi->UseSliceScaleFactor %d\n", sdfi->UseSliceScaleFactor);

  int DoRescale = 1;
  int RescaleNeeded = 0;
  // Go through all files and determine whether any need to rescale
  for (nthfile = 0; nthfile < nlist; nthfile++){
    if(sdfi_list[nthfile]->RescaleSlope != 1 || sdfi_list[nthfile]->RescaleIntercept != 0){
      RescaleNeeded = 1;
      break;
    }
  }
  if(RescaleNeeded){
    // must explicitly set FS_RESCALE_DICOM=0 to turn off rescaling
    if(getenv("FS_RESCALE_DICOM") != NULL && strcmp(getenv("FS_RESCALE_DICOM"),"0")==0){
      printf("INFO: nontrivial rescale factors are present but will not be applied because\n");
      printf("the FS_RESCALE_DICOM environment variable is set to 0.\n");
      printf("If you want to apply rescaling (intercept and slope), then unset that \n");
      printf("environment variable (or set it to non-zero) and re-run\n");
      printf("\n");
      DoRescale = 0;
    }
  }
  else {
    DoRescale = 0;
    printf("INFO: rescale not needed\n");
  }
  if(DoRescale){
    printf("INFO: applying rescale intercept and slope based on (0028,1052) (0028,1053).\n");
    printf("  If you do not want this, then set FS_RESCALE_DICOM to 0 and rerun.\n");
    vol_datatype = MRI_FLOAT;
  }
  // Use float if largest value is greater than short max
  for (nthfile = 0; nthfile < nlist; nthfile++)
    if (sdfi_list[nthfile]->LargestValue >= pow(2.0, 15)) vol_datatype = MRI_FLOAT;
  printf("datatype = %d, short=%d, float=%d\n", vol_datatype, MRI_SHORT, MRI_FLOAT);

  /** Allocate an MRI structure **/
  if (LoadVolume) {
    if (nthonly < 0)
      vol = MRIallocSequence(ncols, nrows, nslices, vol_datatype, nframes);
    else
      vol = MRIallocSequence(ncols, nrows, nslices, vol_datatype, 1);
    if (vol == NULL) {
      fprintf(stderr, "ERROR: could not alloc MRI volume\n");
      fflush(stderr);
      return (NULL);
    }
  }
  else {
    vol = MRIallocHeader(ncols, nrows, nslices, vol_datatype, nframes);
    if (vol == NULL) {
      fprintf(stderr, "ERROR: could not alloc MRI header \n");
      fflush(stderr);
      return (NULL);
    }
  }

  /* set the various paramters for the mri structure */
  vol->xsize = sdfi->VolRes[0]; /* x = col */
  vol->ysize = sdfi->VolRes[1]; /* y = row */
  vol->zsize = sdfi->VolRes[2]; /* z = slice */
  vol->x_r = sdfi->Vc[0];
  vol->x_a = sdfi->Vc[1];
  vol->x_s = sdfi->Vc[2];
  vol->y_r = sdfi->Vr[0];
  vol->y_a = sdfi->Vr[1];
  vol->y_s = sdfi->Vr[2];
  vol->z_r = sdfi->Vs[0];
  vol->z_a = sdfi->Vs[1];
  vol->z_s = sdfi->Vs[2];
  vol->c_r = sdfi->VolCenter[0];
  vol->c_a = sdfi->VolCenter[1];
  vol->c_s = sdfi->VolCenter[2];
  vol->ras_good_flag = 1;
  vol->te = sdfi->EchoTime;
  vol->ti = sdfi->InversionTime;
  vol->flip_angle = sdfi->FlipAngle;
  vol->FieldStrength = sdfi->FieldStrength;
  if (!sdfi->IsMosaic)
    vol->tr = sdfi->RepetitionTime;
  else {
    /* The TR definition will depend upon the software version */
    tmpstring = sdcmExtractNumarisVer(sdfi->NumarisVer, &Maj, &Min, &MinMin);
    if (tmpstring != NULL) {
      printf("Numaris Version: %s Maj = %d, Min=%d, MinMin = %d \n", sdfi->NumarisVer, Maj, Min, MinMin);
    }
    if (tmpstring != NULL && (Min == 1 && MinMin <= 6) && Maj < 4) {
      // This should only be run for pretty old data. I've lost
      // track as to which versions should do this. With Maj<4,
      // I'm pretty sure that this section of code will never
      // be run.  It might need to be run with version 4VA16
      // and earlier.
      printf("Computing TR with number of slices\n");
      vol->tr = sdfi->RepetitionTime * (sdfi->VolDim[2]);
    }
    else
      vol->tr = sdfi->RepetitionTime;
    /* Need to add any gap (eg, as in a hammer sequence */
    printf("Repetition Time = %g, TR = %g ms\n", sdfi->RepetitionTime, vol->tr);
    if (tmpstring != NULL) free(tmpstring);
  }

  // Phase Enc Direction
  if (sdfi->PhEncDir == NULL)
    vol->pedir = strcpyalloc("UNKNOWN");
  else {
    vol->pedir = strcpyalloc(sdfi->PhEncDir);
    str_toupper(vol->pedir);
  }
  printf("PE Dir %s %s\n", sdfi->PhEncDir, vol->pedir);
  fflush(stdout);

  // Load the AutoAlign Matrix, if one is there
  vol->AutoAlign = sdcmAutoAlignMatrix(dcmfile);

  /* Return now if we're not loading pixel data */
  if (!LoadVolume) {
    /* restore progress range */
    global_progress_range[0] = nstart;
    global_progress_range[1] = nend;
    return (vol);
  }

  /* Dump info about the first file to stdout */
  DumpSDCMFileInfo(stdout, sdfi);
  // verification of number of files vs. slices
  if (nlist > nslices * nframes) {
    fprintf(stderr, "ERROR: nlist (%d) > nslices (%d) x frames (%d)\n", nlist, nslices, nframes);
    fprintf(stderr, "ERROR: dump file list into fileinfo.txt.\n");
    fprintf(stderr, "ERROR: check for consistency\n");
    {
      FILE *fp;
      int i;
      fp = fopen("./fileinfo.txt", "w");
      for (i = 0; i < nlist; ++i) {
        DumpSDCMFileInfo(fp, sdfi_list[i]);
      }
      fclose(fp);
    }
    fprintf(stderr, "ERROR: set nlist = nslices*nframes.\n");
    nlist = nslices * nframes;
  }
  /* ------- Go through each file in the Run ---------*/
  global_progress_range[0] = global_progress_range[1];
  global_progress_range[1] += (nend - nstart) / 3;

  printf("UseSliceScaleFactor %d (slice 0: %g)\n", sdfi_list[0]->UseSliceScaleFactor, sdfi_list[0]->SliceScaleFactor);
  MinSliceScaleFactor = sdfi_list[0]->SliceScaleFactor;
  for (nthfile = 0; nthfile < nlist; nthfile++) {
    if (MinSliceScaleFactor > sdfi_list[nthfile]->SliceScaleFactor)
      MinSliceScaleFactor = sdfi_list[nthfile]->SliceScaleFactor;
  }

  IsDWI = 0;
  for (nthfile = 0; nthfile < nlist; nthfile++)
    if (sdfi_list[nthfile]->bval > 0) IsDWI = 1;
  printf("IsDWI = %d\n", IsDWI);
  if (IsDWI) {
    vol->bvals = MatrixAlloc(nframes, 1, MATRIX_REAL);
    vol->bvecs = MatrixAlloc(nframes, 3, MATRIX_REAL);
    vol->bvec_space = BVEC_SPACE_VOXEL;
  }

  env = FSENVgetenv();
  for (nthfile = 0; nthfile < nlist; nthfile++) {
    sdfi = sdfi_list[nthfile];

    /* Handle compression */
    // If changing this code, make sure to change similar code below
    if (strncmp(sdfi->TransferSyntaxUID, jpegCompressed_UID, 19) == 0 ||
        strncmp(sdfi->TransferSyntaxUID, rllEncoded_UID, 19) == 0) {
      // setenv DCMDICTPATH /usr/pubsw/packages/dcmtk/current/share/dcmtk/dicom.dic???
      IsCompressed = 1;
      tmpfile = std::string(env->tmpdir) + "/" + std::string(env->user) + ".tmp.decompressed.dcm.XXXXXX";
      std::vector<char> tmpfilechars(tmpfile.begin(), tmpfile.end());
      tmpfilechars.push_back(0); // Null terminate
      fid = mkstemp(tmpfilechars.data()); // mkstemp updates the "XXXXXX" at the end
      tmpfilechars.pop_back(); // Don't need the null terminator
      tmpfile = std::string(tmpfilechars.begin(), tmpfilechars.end()); // Copy back the modified string
      if (fid == -1) {
        printf("ERROR: could not create temp file for decompression %d\n", fid);
        exit(1);
      }
      close(fid);
      if (strncmp(sdfi->TransferSyntaxUID, jpegCompressed_UID, 19) == 0) {
        printf("JPEG compressed, decompressing\n");
	tmpfilestdout = tmpfile+".dcmdjpeg.out";
	cmd = std::string("fsdcmdecompress --i ")
	  + std::string(sdfi->FileName) + " --o " + tmpfile + " --jpeg >& " + tmpfilestdout;
      }
      if (strncmp(sdfi->TransferSyntaxUID, rllEncoded_UID, 19) == 0) {
        printf("RLE compressed, decompressing\n");
	tmpfilestdout = tmpfile + ".dcmdlrf.out";
	cmd = std::string("fsdcmdecompress --i ")
	  + std::string(sdfi->FileName) + " --o " + tmpfile + " --rle >& " + tmpfilestdout;
      }
      printf("cd %s\n", env->cwd);
      printf("%s\n", cmd.c_str());
      err = system(cmd.c_str()); // Gulp
      if (err == -1) {
        printf("ERROR: %d, see %s for more details\n", err, tmpfilestdout.c_str());
        exit(1);
      }
      FileNameUse = tmpfile; // Hmmmm
    }
    else {
      IsCompressed = 0;
      FileNameUse = std::string(sdfi->FileName);
    }

    /* Get the pixel data */
    element = GetElementFromFile(FileNameUse.c_str(), 0x7FE0, 0x10);
    if (element == NULL) {
      printf("ERROR: reading pixel data from %s\n", FileNameUse.c_str());
      MRIfree(&vol);
      exit(1);
    }
    if (IsCompressed) {
      unlink(tmpfile.c_str());
      unlink(tmpfilestdout.c_str());
    }

    pixeldata = (unsigned short *)(element->d.string);

    if (!IsMosaic) {
      /*---------------------------------------------*/
      /* It's not a mosaic -- load rows and cols from pixel data */
      if (nthfile == 0) {
        frame = 0;
        slice = 0;
      }
      if (Gdiag_no > 0) {
        printf("%3d %3d %3d    %s   %6.1f %6.1f %6.1f\n",
               nthfile,
               slice,
               frame,
               sdfi->FileName,
               sdfi->ImgPos[0],
               sdfi->ImgPos[1],
               sdfi->ImgPos[2]);
        fflush(stdout);
      }
      if (nthonly < 0) {  // do all frames
        for (row = 0; row < nrows; row++) {
          for (col = 0; col < ncols; col++) {
            if (sdfi->UseSliceScaleFactor) {
              /* PW 2012/09/06: Old way to scale data to the dynamic range of
                                a SHORT, but now we're saving as FLOATs.. */
              // val = 8*MinSliceScaleFactor*(*(pixeldata++))/sdfi->SliceScaleFactor;
              val = ((float)*(pixeldata++)) / sdfi->SliceScaleFactor;
            }
            else
              val = *(pixeldata++);
	    if(DoRescale)
	      val = val*sdfi->RescaleSlope + sdfi->RescaleIntercept;
            MRIsetVoxVal(vol, col, row, slice, frame, val);
          }
        }
      }
      // only copy a particular frame
      else if (frame == nthonly) {
        for (row = 0; row < nrows; row++) {
          for (col = 0; col < ncols; col++) {
            if (sdfi->UseSliceScaleFactor)
              val = ((float)*(pixeldata++)) / sdfi->SliceScaleFactor;
            else
              val = *(pixeldata++);
	    if(DoRescale)
	      val = val*sdfi->RescaleSlope + sdfi->RescaleIntercept;
            MRIsetVoxVal(vol, col, row, slice, 0, val);
          }
        }
      }
      if (IsDWI) {
        vol->bvals->rptr[frame + 1][1] = sdfi->bval;
        vol->bvecs->rptr[frame + 1][1] = sdfi->bvecx;
        vol->bvecs->rptr[frame + 1][2] = sdfi->bvecy;
        vol->bvecs->rptr[frame + 1][3] = sdfi->bvecz;
      }
      frame++;
      if (frame >= nframes) {
        frame = 0;
        slice++;
      }
    }
    else {  // is a mosaic
      /*---------------------------------------------*/
      /* It is a mosaic -- load entire volume for this frame from pixel data */
      frame = nthfile;
      nmoscols = sdfi->NImageCols;
      nmosrows = sdfi->NImageRows;
      for (row = 0; row < nrows; row++) {
        for (col = 0; col < ncols; col++) {
          for (slice = 0; slice < nslices; slice++) {
            /* compute the mosaic col and row from the volume
               col, row , and slice */
            err = VolSS2MosSS(col, row, slice, ncols, nrows, nmoscols, nmosrows, &moscol, &mosrow, &OutOfBounds);
            if (err || OutOfBounds) {
              FreeElementData(element);
              free(element);
              MRIfree(&vol);
              exit(1);
            }
            /* Compute the linear index into the block of pixel data */
            mosindex = moscol + mosrow * nmoscols;
            if (sdfi->UseSliceScaleFactor)
              val = ((float)*(pixeldata + mosindex)) / sdfi->SliceScaleFactor;
            else
              val = *(pixeldata + mosindex);
	    if(DoRescale)
	      val = val*sdfi->RescaleSlope + sdfi->RescaleIntercept;
            MRIsetVoxVal(vol, col, row, slice, frame, val);
          }
        }
      }
      if (IsDWI) {
        vol->bvals->rptr[frame + 1][1] = sdfi->bval;
        vol->bvecs->rptr[frame + 1][1] = sdfi->bvecx;
        vol->bvecs->rptr[frame + 1][2] = sdfi->bvecy;
        vol->bvecs->rptr[frame + 1][3] = sdfi->bvecz;
      }
    }

    FreeElementData(element);
    free(element);
    exec_progress_callback(nthfile, nlist, 0, 1);
  } /* for nthfile */

  /* Determine whether Siemens has reversed the slice order prior to
     packing into mosaic. This makes it inconsistent with the geometry
     in the dicom. If so, reverse it back without changing the geometry.
     Note that this does not have anything to do with the order in which
     the slices were acquired. However, if you assume a certain acquistion
     slice order model, it will be wrong unless this reversal is done.  */
  if (sdfiIsSliceOrderReversed(sdfi_list[0])) {
    printf("INFO: detected a Siemens slice order reversal in mosaic, so\n");
    printf(
        "      I'm going to reverse them back "
        "to their 'original' order.\n");
    voltmp = MRIreverseSliceOrder(vol, NULL);
    MRIfree(&vol);
    vol = voltmp;
  }
  else
    printf("INFO: no Siemens slice order reversal detected (good!). \n");

  if (IsDWI) {
    if (env->desired_bvec_space == BVEC_SPACE_SCANNER) {
      printf("Converting bvec to scanner space\n");
      DTIbvecChangeSpace(vol, BVEC_SPACE_SCANNER);
    }
  }

  while (nlist--) {
    // free strings
    if (sdfi_list[nlist]->FileName != NULL) free(sdfi_list[nlist]->FileName);
    if (sdfi_list[nlist]->StudyDate != NULL) free(sdfi_list[nlist]->StudyDate);
    if (sdfi_list[nlist]->StudyTime != NULL) free(sdfi_list[nlist]->StudyTime);
    if (sdfi_list[nlist]->PatientName != NULL) free(sdfi_list[nlist]->PatientName);
    if (sdfi_list[nlist]->SeriesTime != NULL) free(sdfi_list[nlist]->SeriesTime);
    if (sdfi_list[nlist]->AcquisitionTime != NULL) free(sdfi_list[nlist]->AcquisitionTime);
    if (sdfi_list[nlist]->ScannerModel != NULL) free(sdfi_list[nlist]->ScannerModel);
    if (sdfi_list[nlist]->NumarisVer != NULL) free(sdfi_list[nlist]->NumarisVer);
    if (sdfi_list[nlist]->PulseSequence != NULL) free(sdfi_list[nlist]->PulseSequence);
    if (sdfi_list[nlist]->ProtocolName != NULL) free(sdfi_list[nlist]->ProtocolName);
    if (sdfi_list[nlist]->PhEncDir != NULL) free(sdfi_list[nlist]->PhEncDir);
    if (sdfi_list[nlist]->TransferSyntaxUID != NULL) free(sdfi_list[nlist]->TransferSyntaxUID);
    // free struct
    free(sdfi_list[nlist]);
  }
  free(sdfi_list);

  /* restore progress range */
  global_progress_range[0] = nstart;
  global_progress_range[1] = nend;

  FSENVfree(&env);

  return (vol);
}

/*
  Note: This is a complete copy of sdcmLoadVolume,
  but applying both the autoscale to each
  frame (even for MOSAIC) and also resulting in a
  float volume, to avoid reduction
  of dynamic range by the cast to short.
  This is what is needed for autoscaled diffusion
  data. TW 03/22/2012
*/
MRI *sdcmLoadVolumeAutoScale(const char *dcmfile, int LoadVolume, int nthonly)
{
  SDCMFILEINFO *sdfi;
  DCM_ELEMENT *element;
  SDCMFILEINFO **sdfi_list;
  int nthfile;
  int nlist;
  int ncols, nrows, nslices, nframes;
  int nmoscols, nmosrows;
  int mosrow, moscol, mosindex;
  int err, OutOfBounds, IsMosaic;
  int row, col, slice, frame;
  unsigned short *pixeldata;
  MRI *vol, *voltmp;
  char **SeriesList;
  char *tmpstring, *pc = NULL, *pc2 = NULL;
  int Maj, Min, MinMin;
  double xs, ys, zs, xe, ye, ze, d, val;
  int nnlist;
  // int nthdir;
  DTI *dti;
  int TryDTI = 1, DoDTI = 1;
  extern int sliceDirCosPresent;  // set when no ascii header

  xs = ys = zs = xe = ye = ze = d = 0.; /* to avoid compiler warnings */
  slice = 0;
  frame = 0; /* to avoid compiler warnings */

  sliceDirCosPresent = 0;  // assume not present

  /* split progress to 3 parts */
  int nstart = global_progress_range[0];
  int nend = global_progress_range[1];
  global_progress_range[1] = nstart + (nend - nstart) / 3;
  if (SDCMListFile != NULL) {
    SeriesList = ReadSiemensSeries(SDCMListFile, &nlist, dcmfile);
  }
  else {
    SeriesList = ScanSiemensSeries(dcmfile, &nlist);
  }

  if (SeriesList == NULL) {
    fprintf(stderr, "ERROR: could not find any files (SeriesList==NULL)\n");
    return (NULL);
  }

  if (nlist == 0) {
    fprintf(stderr, "ERROR: could not find any files (nlist==0)\n");
    return (NULL);
  }

  fprintf(stderr, "WARNING: YOU ARE USING A BETA AUTOSCALE VERSION !\n");

  // for(nnlist=0; nnlist<nlist; nnlist++) fprintf(stdout,"%3d  %s\n",
  // nnlist,SeriesList[nnlist]);
  // fflush(stdout);

  global_progress_range[0] = global_progress_range[1];
  global_progress_range[1] += (nend - nstart) / 3;
  printf("INFO: loading series header info.\n");
  sdfi_list = LoadSiemensSeriesInfo(SeriesList, nlist);

  // free memory
  nnlist = nlist;
  while (nnlist--) {
    free(SeriesList[nnlist]);
  }
  free(SeriesList);

  printf("INFO: sorting.\n");
  SortSDCMFileInfo(sdfi_list, nlist);

  sdfiAssignRunNo2(sdfi_list, nlist);

  /* First File in the Run */
  if (nthonly < 0) {
    sdfi = sdfi_list[0];
  }
  else {
    sdfi = sdfi_list[nthonly];
  }

  /* There are some Siemens files don't have the slice dircos
     because they dont have ascii header (anonymization?).
     If this happens on a mosiac, then you cannot unpack it.
     For non-mosaics, recompute the slice dircos based on
     the image position of the first and last file. Note:
     slice dircos may be used to sort the files in the first
     place, but this should work either way.
  */
  if (sliceDirCosPresent == 0) {
    printf("\n");
    printf(
        "WARNING: file %s does not contain a Siemens ASCII header\n"
        "has this file been anonymized?\n",
        sdfi_list[0]->FileName);
    if (sdfi_list[0]->IsMosaic) {
      printf("ERROR: cannot unpack mosiacs without ASCII header\n");
      exit(1);
    }
    printf("Proceeding as best as I can ... \n");
    printf("\n");
    xs = sdfi_list[0]->ImgPos[0];
    ys = sdfi_list[0]->ImgPos[1];
    zs = sdfi_list[0]->ImgPos[2];
    xe = sdfi_list[nlist - 1]->ImgPos[0];
    ye = sdfi_list[nlist - 1]->ImgPos[1];
    ze = sdfi_list[nlist - 1]->ImgPos[2];
    d = sqrt((xe - xs) * (xe - xs) + (ye - ys) * (ye - ys) + (ze - zs) * (ze - zs));
    sdfi->Vs[0] = (xe - xs) / d;
    sdfi->Vs[1] = (ye - ys) / d;
    sdfi->Vs[2] = (ze - zs) / d;
    // Note: no need to change to RAS because ImagePos is already RAS
  }
  sdfiFixImagePosition(sdfi);
  sdfiVolCenter(sdfi);

  /* for easy access */
  ncols = sdfi->VolDim[0];
  nrows = sdfi->VolDim[1];
  nslices = sdfi->VolDim[2];
  nframes = sdfi->NFrames;
  IsMosaic = sdfi->IsMosaic;

  // verify nthframe is within the range
  if (nthonly >= 0) {
    if (nthonly > nframes - 1) {
      printf("ERROR: only has %d frames (%d - %d) but called for %d\n", nframes, 0, nframes - 1, nthonly);
      fflush(stdout);
      return (NULL);
    }
  }

  printf("INFO: (%3d %3d %3d), nframes = %d, ismosaic=%d\n", ncols, nrows, nslices, nframes, IsMosaic);
  fflush(stdout);

  /** Allocate an MRI structure **/
  if (LoadVolume) {
    if (nthonly < 0) {
      vol = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, nframes);
    }
    else {
      vol = MRIallocSequence(ncols, nrows, nslices, MRI_FLOAT, 1);
    }
    if (vol == NULL) {
      fprintf(stderr, "ERROR: could not alloc MRI volume\n");
      fflush(stderr);
      return (NULL);
    }
  }
  else {
    vol = MRIallocHeader(ncols, nrows, nslices, MRI_FLOAT, nframes);
    if (vol == NULL) {
      fprintf(stderr, "ERROR: could not alloc MRI header \n");
      fflush(stderr);
      return (NULL);
    }
  }

  /* set the various paramters for the mri structure */
  vol->xsize = sdfi->VolRes[0]; /* x = col */
  vol->ysize = sdfi->VolRes[1]; /* y = row */
  vol->zsize = sdfi->VolRes[2]; /* z = slice */
  vol->x_r = sdfi->Vc[0];
  vol->x_a = sdfi->Vc[1];
  vol->x_s = sdfi->Vc[2];
  vol->y_r = sdfi->Vr[0];
  vol->y_a = sdfi->Vr[1];
  vol->y_s = sdfi->Vr[2];
  vol->z_r = sdfi->Vs[0];
  vol->z_a = sdfi->Vs[1];
  vol->z_s = sdfi->Vs[2];
  vol->c_r = sdfi->VolCenter[0];
  vol->c_a = sdfi->VolCenter[1];
  vol->c_s = sdfi->VolCenter[2];
  vol->ras_good_flag = 1;
  vol->te = sdfi->EchoTime;
  vol->ti = sdfi->InversionTime;
  vol->flip_angle = sdfi->FlipAngle;
  vol->FieldStrength = sdfi->FieldStrength;
  if (!sdfi->IsMosaic) {
    vol->tr = sdfi->RepetitionTime;
  }
  else {
    /* The TR definition will depend upon the software version */
    tmpstring = sdcmExtractNumarisVer(sdfi->NumarisVer, &Maj, &Min, &MinMin);
    if (tmpstring != NULL) {
      printf("Numaris Version: %s Maj = %d, Min=%d, MinMin = %d \n", sdfi->NumarisVer, Maj, Min, MinMin);
    }
    if (tmpstring != NULL && (Min == 1 && MinMin <= 6) && Maj < 4) {
      // This should only be run for pretty old data. I've lost
      // track as to which versions should do this. With Maj<4,
      // I'm pretty sure that this section of code will never
      // be run.  It might need to be run with version 4VA16
      // and earlier.
      printf("Computing TR with number of slices\n");
      vol->tr = sdfi->RepetitionTime * (sdfi->VolDim[2]);
    }
    else {
      vol->tr = sdfi->RepetitionTime;
    }
    /* Need to add any gap (eg, as in a hammer sequence */
    printf("Repetition Time = %g, TR = %g ms\n", sdfi->RepetitionTime, vol->tr);
    if (tmpstring != NULL) {
      free(tmpstring);
    }
  }

  // Phase Enc Direction
  if (sdfi->PhEncDir == NULL) {
    vol->pedir = strcpyalloc("UNKNOWN");
  }
  else {
    vol->pedir = strcpyalloc(sdfi->PhEncDir);
    str_toupper(vol->pedir);
  }
  printf("PE Dir %s %s\n", sdfi->PhEncDir, vol->pedir);
  fflush(stdout);

  // Load the AutoAlign Matrix, if one is there
  vol->AutoAlign = sdcmAutoAlignMatrix(dcmfile);

  // Load DTI bvecs/bvals, if you can. If the procedure fails for some reason
  // you can setenv UNPACK_MGH_DTI 0.
  pc = SiemensAsciiTag(dcmfile, "sDiffusion.lDiffDirections", 0);
  pc2 = SiemensAsciiTag(dcmfile, "sWiPMemBlock.alFree[8]", 0);
  if (pc != NULL && pc2 != NULL) {
    printf("This looks like an MGH DTI volume\n");
    if (getenv("UNPACK_MGH_DTI") != NULL) {
      sscanf(getenv("UNPACK_MGH_DTI"), "%d", &TryDTI);
    }
    else {
      TryDTI = 1;
    }
    if (TryDTI) {
      DoDTI = 1;
    }
    else {
      DoDTI = 0;
    }
    if (!DoDTI) {
      printf("  but not getting bvec info because UNPACK_MGH_DTI is 0\n");
    }
  }
  else {
    DoDTI = 0;
  }
  if (DoDTI && !sdfi->IsMosaic) {
    printf("DTI is in non-moasic form, so cannot extract bvals/bvects\n");
    DoDTI = 0;
  }
  if (DoDTI) {
    printf("MGH DTI SeqPack Info\n");
    // Get b Values from header, based on sequence name.
    // Problem: nthfile = nthvolume when mosaics are used,
    // but not for non-mosaics.
    vol->bvals = MatrixAlloc(nframes, 1, MATRIX_REAL);
    for (nthfile = 0; nthfile < nlist; nthfile++) {
      // Go thru all the files in order to get all the directions
      sdfi = sdfi_list[nthfile];
      DTIparsePulseSeqName(sdfi->PulseSequence, &sdfi->bValue, &sdfi->nthDirection);
      // nthdir = sdfi->nthDirection;
      vol->bvals->rptr[nthfile + 1][1] = sdfi->bValue;
      printf("%d %s %lf %d\n", nthfile, sdfi->PulseSequence, sdfi->bValue, sdfi->nthDirection);
    }
    // Have to get vectors from archive
    dti = DTIstructFromSiemensAscii(dcmfile);
    if (dti == NULL) {
      printf("There was an error when tyring to load the gradient directions\n");
      printf("If you 'setenv UNPACK_MGH_DTI 0', it will not attempt to load\n");
      printf("the gradients\n");
      exit(1);
    }
    // vol->bvals = MatrixCopy(dti->bValue,NULL);
    vol->bvecs = MatrixCopy(dti->GradDir, NULL);
    DTIfree(&dti);
  }
  if (pc) {
    free(pc);
  }
  if (pc2) {
    free(pc2);
  }

  /* Return now if we're not loading pixel data */
  if (!LoadVolume) {
    /* restore progress range */
    global_progress_range[0] = nstart;
    global_progress_range[1] = nend;
    return (vol);
  }

  if (strcmp(sdfi->TransferSyntaxUID, "1.2.840.10008.1.2.4.70") == 0) {
    printf("ERROR: the pixel data cannot be loaded as it is JPEG compressed.\n");
    printf("       (Transfer Syntax UID: %s)\n", sdfi->TransferSyntaxUID);
    exit(1);
  }

  /* Dump info about the first file to stdout */
  DumpSDCMFileInfo(stdout, sdfi);
  // verification of number of files vs. slices
  if (nlist > nslices * nframes) {
    fprintf(stderr, "ERROR: nlist (%d) > nslices (%d) x frames (%d)\n", nlist, nslices, nframes);
    fprintf(stderr, "ERROR: dump file list into fileinfo.txt.\n");
    fprintf(stderr, "ERROR: check for consistency\n");
    {
      FILE *fp;
      int i;
      fp = fopen("./fileinfo.txt", "w");
      for (i = 0; i < nlist; ++i) {
        DumpSDCMFileInfo(fp, sdfi_list[i]);
      }
      fclose(fp);
    }
    fprintf(stderr, "ERROR: set nlist = nslices*nframes.\n");
    nlist = nslices * nframes;
  }
  /* ------- Go through each file in the Run ---------*/
  global_progress_range[0] = global_progress_range[1];
  global_progress_range[1] += (nend - nstart) / 3;

  /*
  printf("UseSliceScaleFactor %d (slice 0: %g)\n",sdfi_list[0]->UseSliceScaleFactor,
   sdfi_list[0]->SliceScaleFactor);
  MinSliceScaleFactor = sdfi_list[0]->SliceScaleFactor;
  for(nthfile = 0; nthfile < nlist; nthfile ++)
    if(MinSliceScaleFactor > sdfi_list[nthfile]->SliceScaleFactor)
      MinSliceScaleFactor = sdfi_list[nthfile]->SliceScaleFactor;
  */

  for (nthfile = 0; nthfile < nlist; nthfile++) {
    float ascale_factor = 1.0;

    sdfi = sdfi_list[nthfile];

    /* TW get the autoscale parameter here */
    element = GetElementFromFile(sdfi->FileName, 0x0020, 0x4000);
    if (strncmp(element->d.string, "Scale Factor:", strlen("Scale Factor:")) == 0) {
      ascale_factor = (float)strtod(&(element->d.string[strlen("Scale Factor:") + 1]), NULL);
    }

    /* Get the pixel data */
    element = GetElementFromFile(sdfi->FileName, 0x7FE0, 0x10);
    if (element == NULL) {
      fprintf(stderr, "ERROR: reading pixel data from %s\n", sdfi->FileName);
      MRIfree(&vol);
    }
    pixeldata = (unsigned short *)(element->d.string);

    if (!IsMosaic) {
      /*---------------------------------------------*/
      /* It's not a mosaic -- load rows and cols from pixel data */
      if (nthfile == 0) {
        frame = 0;
        slice = 0;
      }
      if (Gdiag_no > 0) {
        printf("%3d %3d %3d    %s   %6.1f %6.1f %6.1f\n",
               nthfile,
               slice,
               frame,
               sdfi->FileName,
               sdfi->ImgPos[0],
               sdfi->ImgPos[1],
               sdfi->ImgPos[2]);
        fflush(stdout);
      }
      if (nthonly < 0) {
        for (row = 0; row < nrows; row++) {
          for (col = 0; col < ncols; col++) {
            val = *(pixeldata++);
            MRIFseq_vox(vol, col, row, slice, frame) = (float)(val) / ascale_factor;
          }
        }
      }
      // only copy a particular frame
      else if (frame == nthonly) {
        for (row = 0; row < nrows; row++) {
          for (col = 0; col < ncols; col++) {
            MRIFvox(vol, col, row, slice) = (float)(*(pixeldata++)) / ascale_factor;
          }
        }
      }
      frame++;
      if (frame >= nframes) {
        frame = 0;
        slice++;
      }
    }
    else {
      /*---------------------------------------------*/
      /* It is a mosaic -- load entire volume for this frame from pixel data */
      frame = nthfile;
      nmoscols = sdfi->NImageCols;
      nmosrows = sdfi->NImageRows;
      for (row = 0; row < nrows; row++) {
        for (col = 0; col < ncols; col++) {
          for (slice = 0; slice < nslices; slice++) {
            /* compute the mosaic col and row from the volume
               col, row , and slice */
            err = VolSS2MosSS(col, row, slice, ncols, nrows, nmoscols, nmosrows, &moscol, &mosrow, &OutOfBounds);
            if (err || OutOfBounds) {
              FreeElementData(element);
              free(element);
              MRIfree(&vol);
              exit(1);
            }
            /* Compute the linear index into the block of pixel data */
            mosindex = moscol + mosrow * nmoscols;
            MRIFseq_vox(vol, col, row, slice, frame) = (float)(*(pixeldata + mosindex)) / ascale_factor;
          }
        }
      }
    }

    FreeElementData(element);
    free(element);
    exec_progress_callback(nthfile, nlist, 0, 1);
  } /* for nthfile */

  /* Determine whether Siemens has reversed the slice order prior to
     packing into mosaic. This makes it inconsistent with the geometry
     in the dicom. If so, reverse it back without changing the geometry.
     Note that this does not have anything to do with the order in which
     the slices were acquired. However, if you assume a certain acquistion
     slice order model, it will be wrong unless this reversal is done.  */
  if (sdfiIsSliceOrderReversed(sdfi_list[0])) {
    printf("INFO: detected a Siemens slice order reversal in mosaic, so\n");
    printf(
        "      I'm going to reverse them back "
        "to their 'original' order.\n");
    voltmp = MRIreverseSliceOrder(vol, NULL);
    MRIfree(&vol);
    vol = voltmp;
  }
  else {
    printf("INFO: no Siemens slice order reversal detected (good!). \n");
  }

  while (nlist--) {
    // free strings
    if (sdfi_list[nlist]->FileName != NULL) {
      free(sdfi_list[nlist]->FileName);
    }
    if (sdfi_list[nlist]->StudyDate != NULL) {
      free(sdfi_list[nlist]->StudyDate);
    }
    if (sdfi_list[nlist]->StudyTime != NULL) {
      free(sdfi_list[nlist]->StudyTime);
    }
    if (sdfi_list[nlist]->PatientName != NULL) {
      free(sdfi_list[nlist]->PatientName);
    }
    if (sdfi_list[nlist]->SeriesTime != NULL) {
      free(sdfi_list[nlist]->SeriesTime);
    }
    if (sdfi_list[nlist]->AcquisitionTime != NULL) {
      free(sdfi_list[nlist]->AcquisitionTime);
    }
    if (sdfi_list[nlist]->ScannerModel != NULL) {
      free(sdfi_list[nlist]->ScannerModel);
    }
    if (sdfi_list[nlist]->NumarisVer != NULL) {
      free(sdfi_list[nlist]->NumarisVer);
    }
    if (sdfi_list[nlist]->PulseSequence != NULL) {
      free(sdfi_list[nlist]->PulseSequence);
    }
    if (sdfi_list[nlist]->ProtocolName != NULL) {
      free(sdfi_list[nlist]->ProtocolName);
    }
    if (sdfi_list[nlist]->PhEncDir != NULL) {
      free(sdfi_list[nlist]->PhEncDir);
    }
    if (sdfi_list[nlist]->TransferSyntaxUID != NULL) {
      free(sdfi_list[nlist]->TransferSyntaxUID);
    }
    // free struct
    free(sdfi_list[nlist]);
  }
  free(sdfi_list);

  /* restore progress range */
  global_progress_range[0] = nstart;
  global_progress_range[1] = nend;

  return (vol);
}

/*!/
  \fn MATRIX *sdcmAutoAlignMatrix(const char *dcmfile)
  \brief Extracts the Auto Align Matrix from the Siemens ascii header.
*/
MATRIX *sdcmAutoAlignMatrix(const char *dcmfile)
{
  char *tmpstr;
  char sdcmtag[1000];
  int n, row, col;
  MATRIX *M;
  double v;

  /* Check for the existence of the 15th matrix element. If the
     auto align matrix exists, then this item must be 1. */
  tmpstr = SiemensAsciiTagEx(dcmfile, "sAutoAlign.dAAMatrix[15]", 0);
  if (tmpstr == NULL) {
    return (NULL);
  }

  printf("AutoAlign matrix detected \n");
  free(tmpstr);
  M = MatrixAlloc(4, 4, MATRIX_REAL);
  n = 0;
  for (row = 1; row <= 4; row++) {
    for (col = 1; col <= 4; col++) {
      sprintf(sdcmtag, "sAutoAlign.dAAMatrix[%d]", n);
      n++;
      tmpstr = SiemensAsciiTagEx(dcmfile, sdcmtag, 0);
      if (tmpstr == NULL) {
        v = 0;
      }
      else {
        sscanf(tmpstr, "%lf", &v);
        free(tmpstr);
      }
      M->rptr[row][col] = v;
      if (Gdiag_no > 0) {
        printf("n=%2d c=%d r=%d %lf ----------------------\n", n, col, row, v);
        printf("%s\n", sdcmtag);
        if (tmpstr) {
          printf("%s", tmpstr);
        }
      }
    }
  }

  printf("AutoAlign Matrix --------------------- \n");
  MatrixPrint(stdout, M);
  printf("\n");

  return (M);
}

/*---------------------------------------------------------------
  GetElementFromFile() - gets an element from a DICOM file. Returns
  a pointer to the object (or NULL upon failure).
  Author: Douglas Greve 9/6/2001
  ---------------------------------------------------------------*/
DCM_ELEMENT *GetElementFromFile(const char *dicomfile, long grpid, long elid)
{
  DCM_OBJECT *object = 0;
  CONDITION cond;
  DCM_ELEMENT *element;
  DCM_TAG tag;
  unsigned int rtnLength;
  void *Ctx = NULL;

  element = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));

  object = GetObjectFromFile(dicomfile, 0);
  if (object == NULL) {
    exit(1);
  }

  tag = DCM_MAKETAG(grpid, elid);
  cond = DCM_GetElement(&object, tag, element);
  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    free(element);
    return (NULL);
  }
  AllocElementData(element);
  cond = DCM_GetElementValue(&object, element, &rtnLength, &Ctx);
  /* Does Ctx have to be freed? */
  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    FreeElementData(element);
    free(element);
    return (NULL);
  }
  DCM_CloseObject(&object);

  COND_PopCondition(1); /********************************/

  return (element);
}
/*---------------------------------------------------------------
  GetObjectFromFile() - gets an object from a DICOM file. Returns
  a pointer to the object (or NULL upon failure).
  Author: Douglas Greve
  ---------------------------------------------------------------*/
DCM_OBJECT *GetObjectFromFile(const char *fname, unsigned long options)
{
  CONDITION cond;
  DCM_OBJECT *object = 0;
  int ok;

  // printf("     GetObjectFromFile(): %s %ld\n",fname,options);
  fflush(stdout);
  fflush(stderr);

  ok = IsDICOM(fname);
  if (!ok) {
    fprintf(stderr, "ERROR: %s is not a dicom file\n", fname);
    COND_DumpConditions();
    return (NULL);
  }
  options = options | DCM_ACCEPTVRMISMATCH;

  cond = DCM_OpenFile(fname, DCM_PART10FILE | options, &object);
  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    cond = DCM_OpenFile(fname, DCM_ORDERLITTLEENDIAN | options, &object);
  }
  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    cond = DCM_OpenFile(fname, DCM_ORDERBIGENDIAN | options, &object);
  }
  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    cond = DCM_OpenFile(fname, DCM_FORMATCONVERSION | options, &object);
  }
  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    COND_DumpConditions();
    return (NULL);
  }

  return (object);
}
/*-------------------------------------------------------------------
  AllocElementData() - allocates memory for the data portion of
  the element structure. Returns 1 if there's an error (otherwise 0).
  Use FreeElementData() to free the data portion.
  Author: Douglas N. Greve, 9/6/2001
  -------------------------------------------------------------------*/
int AllocElementData(DCM_ELEMENT *e)
{
  switch (e->representation) {
    case DCM_AE:
    case DCM_AS:
    case DCM_CS:
    case DCM_DA:
    case DCM_DS:
    case DCM_DT:
    case DCM_IS:
    case DCM_LO:
    case DCM_LT:
    case DCM_OB:
    case DCM_OW:
    case DCM_PN:
    case DCM_SH:
    case DCM_UN:  // unknown
    case DCM_ST:
    case DCM_TM:
    case DCM_UI:
      e->d.string = (char *)calloc(e->length + 17, sizeof(char));
      e->d.string[e->length] = '\0'; /* add null terminator */
      break;
    case DCM_SS:
      e->d.ss = (short *)calloc(1, sizeof(short));
      break;
    case DCM_SL:
      e->d.sl = (int *)calloc(1, sizeof(int));
      break;
    case DCM_SQ:
      fprintf(stderr, "Element is of type dcm_sq, not supported\n");
      return (1);
      break;
    case DCM_UL:
      e->d.ul = (unsigned int *)calloc(1, sizeof(unsigned int));
      break;
    case DCM_US:
      e->d.us = (unsigned short *)calloc(1, sizeof(unsigned short));
      // e->length is the byte count
      break;
    case DCM_AT:
      e->d.at = (DCM_TAG *)calloc(1, sizeof(DCM_TAG));
      break;
    case DCM_FD:
      e->d.fd = (double *)calloc(e->multiplicity, sizeof(double));
      break;
    case DCM_FL:
      e->d.fl = (float *)calloc(1, sizeof(float));
      break;
    default:
      fprintf(stderr, "AllocElementData: %d unrecognized\n", e->representation);
      return (1);
  }

  return (0);
}
/*---------------------------------------------------------------
  ElementValueString() - returns the value of the element as a
  null terminated string. Does not parse multiple valued elements.
  For string elements, it just copies the string and adds a
  terminator. For others, it just uses sprrintf.
  ---------------------------------------------------------------*/
char *ElementValueString(DCM_ELEMENT *e, int DoBackslash)
{
  char *evstring=NULL;
  unsigned int n, len;
  std::string tmpstr;
  tmpstr.resize(2048);

  memset(&tmpstr[0], 0, 2000);

  switch (e->representation) {
    case DCM_AE:
    case DCM_AS:
    case DCM_CS:
    case DCM_DA:
    case DCM_DS:
    case DCM_DT:
    case DCM_IS:
    case DCM_LO:
    case DCM_LT:
    case DCM_OW:
    case DCM_PN:
    case DCM_SH:
    case DCM_UN:  // unknown
    case DCM_ST:
    case DCM_TM:
    case DCM_UI:
      tmpstr = e->d.string;
      break;
    case DCM_OB:
      evstring = (char *) calloc(sizeof(char),e->length+1);
      len = e->length;
      for (n = 0; n < e->length; n++) {
        if (isprint(e->d.string[n]))
          evstring[n] = e->d.string[n];
        else if (e->d.string[n] == '\r' || e->d.string[n] == '\n' || e->d.string[n] == '\v')
          evstring[n] = '\n';
        else
          evstring[n] = ' ';
      }
      // printf("%s\n",tmpstr);
      break;
    case DCM_SS:
      tmpstr = std::to_string(static_cast<int>(*(e->d.ss)));
      break;
    case DCM_SL:
      tmpstr = std::to_string(static_cast<long>(*(e->d.sl)));
      break;
    case DCM_UL:
      tmpstr = std::to_string(static_cast<long>(*(e->d.ul)));
      break;
    case DCM_US:
      tmpstr = std::to_string(static_cast<int>(*(e->d.us)));
      break;
    case DCM_AT:
      tmpstr = std::to_string(static_cast<long>(*(e->d.at)));
      break;
    case DCM_FL:
      tmpstr = std::to_string(static_cast<float>(e->d.fd[0]));
      for (n = 1; n < e->multiplicity; n++) {
	tmpstr += " " + std::to_string(static_cast<float>(e->d.fd[n]));
      }
      break;
    case DCM_FD:
      tmpstr = std::to_string(static_cast<double>(e->d.fd[0]));
      for (n = 1; n < e->multiplicity; n++) {
	tmpstr += " " + std::to_string(static_cast<double>(e->d.fd[n]));
      }
      break;
    default:
      fprintf(stderr, "ElementValueString: %d unrecognized", e->representation);
      return (NULL);
  }

  if(evstring==NULL){
    len = tmpstr.size();
    evstring = (char *)calloc(len + 1, sizeof(char));
    memmove(evstring, tmpstr.data(), len + 1);
  }

  if (DoBackslash) {
    // replace backslashes with spaces
    for (n = 0; n < len; n++)
      if (evstring[n] == '\\') evstring[n] = ' ';
  }

  return (evstring);
}

/*-------------------------------------------------------------------
  FreeElementData() - frees memory allocated for the data portion of
  the element structure. Returns 1 if there's an error (otherwise 0).
  See also AllocElementData().
  Author: Douglas N. Greve, 9/6/2001
  -------------------------------------------------------------------*/
int FreeElementData(DCM_ELEMENT *e)
{
  switch (e->representation) {
    case DCM_AE:
    case DCM_AS:
    case DCM_CS:
    case DCM_DA:
    case DCM_DS:
    case DCM_DT:
    case DCM_IS:
    case DCM_LO:
    case DCM_LT:
    case DCM_OB:
    case DCM_OW:
    case DCM_PN:
    case DCM_SH:
    case DCM_ST:
    case DCM_TM:
    case DCM_UI:
      free(e->d.string);
      e->d.string = NULL;
      break;
    case DCM_SS:
      // free(&e->d.ss);
      free(e->d.ss);
      e->d.ss = NULL;
      break;
    case DCM_SL:
      // free(&e->d.sl);
      free(e->d.sl);
      e->d.sl = NULL;
      break;
    case DCM_SQ:
      // free(&e->d.sq);
      free(e->d.sq);
      e->d.sq = NULL;
      break;
    case DCM_UL:
      // free(&e->d.ul);
      free(e->d.ul);
      e->d.ul = NULL;
      break;
    case DCM_US:
      // free(&e->d.us);
      free(e->d.us);
      e->d.us = NULL;
      break;
    case DCM_AT:
      // free(&e->d.at);
      free(e->d.at);
      e->d.at = NULL;
      break;
    case DCM_FD:
      free(e->d.fd);
      break;
    case DCM_FL:
      free(e->d.fl);
      break;
    default:
      fprintf(stderr, "FreeElementData: %d unrecognized", e->representation);
      return (1);
  }

  return (0);
}
/*-------------------------------------------------------------------
  IsSiemensDICOM() - returns 1 if the file is a Siemens DICOM File.
  Returns 0 is its not a dicom file or not a Siemsn DICOM File.
  Checks DICOM tag (8,70).
  Author: Douglas N. Greve, 9/6/2001
  -------------------------------------------------------------------*/
int IsSiemensDICOM(const char *dcmfile)
{
  DCM_ELEMENT *e;

  // printf("Entering IsSiemensDICOM (%s)\n",dcmfile);
  fflush(stdout);
  fflush(stderr);

  if (!IsDICOM(dcmfile)) {
    fflush(stdout);
    fflush(stderr);
    // printf("Leaving IsSiemensDICOM (%s)\n",dcmfile);
    return (0);
  }
  fflush(stdout);
  fflush(stderr);

  e = GetElementFromFile(dcmfile, 0x8, 0x70);
  if (e == NULL) {
    printf(
        "WARNING: searching dicom file %s for "
        "Manufacturer tag 0x8, 0x70\n",
        dcmfile);
    printf("WARNING: the result could be a mess.\n");
    return (0);
  }
  fflush(stdout);
  fflush(stderr);

  /* Siemens appears to add a space onto the end of their
     Manufacturer string*/
  if (strcmp(e->d.string, "SIEMENS") != 0 && strcmp(e->d.string, "SIEMENS ") != 0) {
    FreeElementData(e);
    free(e);
    return (0);
  }
  FreeElementData(e);
  free(e);

  // printf("Leaving IsSiemensDICOM (%s)\n",dcmfile);
  return (1);
}

/* The original SiemensQsciiTag() is too slow         */
/* make sure that returned value be freed if non-null */
char *SiemensAsciiTagEx(const char *dcmfile, const char *TagString, int cleanup)
{
  static char filename[1024] = "";
  static char **lists = 0;
  static int count = 0;
  static int startOfAscii = 0;
  static int MAX_ASCIILIST = 512;
  static int INCREMENT = 64;

  int k;
  FILE *fp;
  char buf[1024];
  char command[1024 + 32];
  char *plist = 0;
  char VariableName[512];
  char *VariableValue = 0;
  char tmpstr2[512];
  int newSize;
  char **newlists = 0;

// Use this when debugging (gdb) because it will not fork
// return(SiemensAsciiTag(dcmfile, TagString, cleanup));

#ifdef Darwin
  // This routine calls the unix "strings" command which operates
  // differently on mac vs linux, so use the old, slow version of this
  // routine if it is a mac.
  return (SiemensAsciiTag(dcmfile, TagString, cleanup));
#endif

  if (getenv("USE_SIEMENSASCIITAG")) return (SiemensAsciiTag(dcmfile, TagString, cleanup));

  // cleanup section.  Make sure to set cleanup =1 at the final call
  // don't rely on TagString but the last flag only
  if (cleanup == 1) {
    for (int i = 0; i < count; ++i) {
      if (lists[i]) {
        free(lists[i]);
        lists[i] = 0;
      }
    }
    free(lists);
    lists = 0;
    count = 0;
    startOfAscii = 0;
    strcpy(filename, "");
    return ((char *)0);
  }

  // if the filename changed, then cache the ascii strings
  if (strcmp(dcmfile, filename) != 0) {
    if (lists == 0) {
      lists = (char **)calloc(sizeof(char *), MAX_ASCIILIST);
    }
    // initialized to be zero

    // Copy dcmfile to filename. If dcmfilename has parentheses then
    // the unix strings command below will fail. The code below puts
    // backslashes in front of any parens
    memset(&filename[0], 0, 1024);
    k = 0;
    for (unsigned int i = 0; i < strlen(dcmfile); i++) {
      if (dcmfile[i] == '(' || dcmfile[i] == ')' || dcmfile[i] == '[' || dcmfile[i] == ']') {
        filename[k] = 92;  // 92 is ascii dec for backslash
        k++;
      }
      filename[k] = dcmfile[i];
      k++;
    }

    // free allocated list of strings
    for (int i = 0; i < count; ++i) {
      if (lists[i]) {
        free(lists[i]);
        lists[i] = 0;
      }
    }
    // now build up string lists ///////////////////////////////////
    startOfAscii = 0;
    strcpy(command, "strings ");
    strcat(command, filename);
    // Note: popen creates a fork, which can be a problem in gdb
    if ((fp = popen(command, "r")) == NULL) {
      fprintf(stderr, "could not open pipe for %s\n", filename);
      return 0;
    }
    count = 0;
    while (fgets(buf, 1024, fp)) {
      // replace CR with null
      char *p = strrchr(buf, '\n');
      if (p) {
        *p = '\0';
      }

      // check the region
      //      if (strncmp(buf, "### ASCCONV BEGIN ###", 21)==0)
      // it seems the connectome scanner has begin fields that
      // do not include the chars ' ###'
      if (strncmp(buf, "### ASCCONV BEGIN", 17) == 0) {
        startOfAscii = 1;
      }
      else if (strncmp(buf, "### ASCCONV END ###", 19) == 0) {
        startOfAscii = 0;
      }

      if (startOfAscii == 1) {
        // printf("%d:%d, %s\n", count, strlen(buf), buf);
        plist = (char *)malloc(strlen(buf) + 1);
        strcpy(plist, buf);
        lists[count] = plist;
        ++count;
        // if we exhaused the list pointer
        // we have to realloc
        if (count == MAX_ASCIILIST) {
          newSize = MAX_ASCIILIST + INCREMENT;
          newlists = (char **)realloc(lists, newSize * sizeof(char *));
          if (newlists != 0)  // if not failed
          {
            lists = newlists;  // update the pointer
            // initialize uninitialized list
            for (int i = 0; i < INCREMENT; ++i) {
              lists[MAX_ASCIILIST + i] = (char *)0;
            }
            MAX_ASCIILIST = newSize;
          }
          else {
            // this should not happen, but hey just in case.
            // Ascii tag is not essential and thus allow it to pass.
            fprintf(stderr,
                    "cannot store any more ASCII tags. "
                    "Tagname is %s\n",
                    TagString);
          }
        }
      }
    }
    pclose(fp);
  }
  // build up string lists available
  // search the tag
  for (int i = 0; i < count; ++i) {
    // get the variable name (the first string)
    VariableName[0] = 0;
    sscanf(lists[i], "%s %*s %*s", VariableName);
    if (VariableName[0] && (strcmp(VariableName, TagString) == 0)) {
      /* match found. get the value (the third string) */
      sscanf(lists[i], "%*s %*s %s", tmpstr2);
      VariableValue = (char *)calloc(strlen(tmpstr2) + 17, sizeof(char));
      memmove(VariableValue, tmpstr2, strlen(tmpstr2));
    }
    else {
      continue;
    }
  }
  // show any errors present
  fflush(stdout);
  fflush(stderr);

  return VariableValue;
}

/*-----------------------------------------------------------------
  SiemensAsciiTag() - siemens dicom files have some data stored as a
  block of ASCII text. Each line in the block has the form:
  VariableName = VariableValue
  The begining of the block is delineated by the line
  ### ASCCONV BEGIN
  The end of the block is delineated by the line
  ### ASCCONV END ###
  This function searches this block for a variable named TagString
  and returns the Value as a string.

  It returns NULL if
  1. It's not a Siemens DICOM File
  2. The begining of the ASCII block cannot be found
  3. There is no match with the TagString

  Note: flag does not do anything. Just there to make compatible
  with SiemensAsciiTagEx(). If flag=1, returns NULL;

  Note: Prior to Feb 25, 2011, the beginning of the ascii header
  was coded by the string "### ASCCONV BEGIN ###" -- same as
  now but with three hash symbols. The Skyra DICOMS changed this
  convention, but it was easily fixed by removing the three hash
  marks.

  Author: Douglas N. Greve, 9/6/2001
  -----------------------------------------------------------------*/
char *SiemensAsciiTag(const char *dcmfile, const char *TagString, int flag)
{
  char linestr[4000];
  char tmpstr2[500];
  FILE *fp;
  int dumpline, nthchar;
  char *rt;
  const char *BeginStr;
  int LenBeginStr;
  char *TestStr;
  int nTest;
  char VariableName[500];
  char *VariableValue;

  if (flag) {
    return (NULL);
  }

  // printf("Entering SiemensAsciiTag() \n");fflush(stdout);fflush(stderr);
  // printf("dcmfile = %s, tagstr = %s\n",dcmfile,TagString);
  fflush(stdout);
  fflush(stderr);

  BeginStr = "### ASCCONV BEGIN";
  LenBeginStr = strlen(BeginStr);
  TestStr = (char *)calloc(LenBeginStr + 1, sizeof(char));

  // Don't require it to be a dicom file -- could be an info dump output
  // from mri_probedicom.
  // if(!IsSiemensDICOM(dcmfile)) return(NULL);

  fp = fopen(dcmfile, "r");
  if (fp == NULL) {
    printf("ERROR: could not open dicom file %s\n", dcmfile);
    exit(1);
  }

  /* This section steps through the file char-by-char until
     the BeginStr is matched */
  dumpline = 0;
  nthchar = 0;
  while (1) {
    fseek(fp, nthchar, SEEK_SET);
    nTest = fread(TestStr, sizeof(char), LenBeginStr, fp);
    if (nTest != LenBeginStr) {
      break;
    }
    if (strcmp(TestStr, BeginStr) == 0) {
      fseek(fp, nthchar, SEEK_SET);
      dumpline = 1;
      break;
    }
    nthchar++;
  }
  free(TestStr);

  if (!dumpline) {
    fclose(fp);  // must close
    fflush(stdout);
    fflush(stderr);
    return (NULL); /* Could not match Begin String */
  }

  /* Once the Begin String has been matched, this section
     searches each line until the TagString is matched
     or until the End String is matched */
  VariableValue = NULL;
  while (1) {
    rt = fgets(linestr, 4000, fp);
    if (rt == NULL) {
      break;
    }

    if (strncmp(linestr, "### ASCCONV END ###", 19) == 0) {
      break;
    }

    VariableName[0] = 0;
    linestr[3999] = 0;
    sscanf(linestr, "%s %*s %*s", VariableName);

    if (strlen(VariableName) != strlen(TagString)) {
      continue;
    }

    if (strcmp(VariableName, TagString) == 0) {
      /* match found */
      sscanf(linestr, "%*s %*s %s", tmpstr2);
      VariableValue = (char *)calloc(strlen(tmpstr2) + 1, sizeof(char));
      memmove(VariableValue, tmpstr2, strlen(tmpstr2));
      break;
    }
  }
  fclose(fp);

  // printf("Leaving SiemensAsciiTag() \n");fflush(stdout);fflush(stderr);
  fflush(stdout);
  fflush(stderr);

  return (VariableValue);
}
/*-----------------------------------------------------------------------
  dcmGetVolRes - Gets the volume resolution (mm) from a DICOM File. The
  column and row resolution is obtained from tag (28,30). This tag is stored
  as a string of the form "ColRes\RowRes". The slice resolution is obtained
  from tag (18,50). The slice thickness may not be correct for mosaics.
  See sdcmMosaicSliceRes().
  Author: Douglas N. Greve, 9/6/2001
  -----------------------------------------------------------------------*/
int dcmGetVolRes(const char *dcmfile, float *ColRes, float *RowRes, float *SliceRes)
{
  DCM_ELEMENT *e;
  char *s;
  int ns, n;
  int slash_not_found;
  int tag_not_found = 0;

  /* Load the Pixel Spacing - this is a string of the form:
     ColRes\RowRes   */
  e = GetElementFromFile(dcmfile, 0x28, 0x30);
  if (e == NULL) {
    return (1);
  }

  /* Put it in a temporary sting */
  s = e->d.string;

  /* Go through each character looking for the backslash */
  slash_not_found = 1;
  ns = strlen(s);
  for (n = 0; n < ns; n++) {
    if (s[n] == '\\') {
      s[n] = ' ';
      slash_not_found = 0;
      break;
    }
  }
  if (slash_not_found) {
    return (1);
  }

  sscanf(s, "%f %f", ColRes, RowRes);

  FreeElementData(e);
  free(e);

  /*E* // Load the Spacing Between Slices
    e = GetElementFromFile(dcmfile, 0x18, 0x88);
    if(e == NULL){
    // If 18,88 does not exist, load the slice thickness
    e = GetElementFromFile(dcmfile, 0x18, 0x50);
    if(e == NULL) return(1);
    }
    sscanf(e->d.string,"%f",SliceRes);
    FreeElementData(e); free(e);
  */

  if (AutoSliceResElTag) {
    printf("Automatically determining SliceResElTag\n");
    e = GetElementFromFile(dcmfile, 0x18, 0x23);
    if (e != NULL) {
      if (strcmp(e->d.string, "3D") == 0)
        SliceResElTag1 = 0x50;
      else
        SliceResElTag1 = 0x88;
    }
    else
      printf("Tag 18,23 is null, cannot automatically determine SliceResElTag\n");
    printf("SliceResElTag order is %lx then %lx\n", SliceResElTag1, SliceResElTag2);
  }
  /* By default, the slice resolution is determined from 18,88. If
     that does not exist, then 18,50 is used. For siemens mag res
     angiogram (MRAs), 18,50 must be used first */
  e = GetElementFromFile(dcmfile, 0x18, SliceResElTag1);
  if (e == NULL)
    tag_not_found = 1;
  else {
    sscanf(e->d.string, "%f", SliceRes);
    if (*SliceRes == 0) {
      tag_not_found = 1;  // tag found but was zero
      FreeElementData(e);
      free(e);
    }
  }
  if (tag_not_found) {  // so either no tag or tag was zero
    e = GetElementFromFile(dcmfile, 0x18, SliceResElTag2);
    if (e == NULL) return (1);  // no tag
    sscanf(e->d.string, "%f", SliceRes);
    FreeElementData(e);
    free(e);
    if (*SliceRes == 0) return (1);  // tag exists but zero
  }
  // fprintf(stderr, "SliceRes=%f\n",*SliceRes);
  // FreeElementData(e); free(e);

  return (0);
}
/*-----------------------------------------------------------------------
  dcmGetSeriesNo - Gets the series number from tag (20,11).
  Returns -1 if error.
  Author: Douglas N. Greve, 9/25/2001
  -----------------------------------------------------------------------*/
int dcmGetSeriesNo(const char *dcmfile)
{
  DCM_ELEMENT *e;
  int SeriesNo;

  e = GetElementFromFile(dcmfile, 0x20, 0x11);
  if (e == NULL) {
    return (-1);
  }

  sscanf(e->d.string, "%d", &SeriesNo);

  FreeElementData(e);
  free(e);

  return (SeriesNo);
}
/*-----------------------------------------------------------------------
  dcmGetNRows - Gets the number of rows in the image from tag (28,10).
  Note that this is the number of rows in the image regardless of
  whether its a mosaic.
  Returns -1 if error.
  Author: Douglas N. Greve, 9/6/2001
  -----------------------------------------------------------------------*/
int dcmGetNRows(const char *dcmfile)
{
  DCM_ELEMENT *e;
  int NRows;

  e = GetElementFromFile(dcmfile, 0x28, 0x10);
  if (e == NULL) {
    return (-1);
  }

  NRows = *(e->d.us);

  if (e->representation != DCM_US) {
    printf("bad element for %s\n", dcmfile);
  }

  FreeElementData(e);
  free(e);

  return (NRows);
}
/*-----------------------------------------------------------------------
  dcmGetNCols - Gets the number of columns in the image from tag (28,11).
  Note that this is the number of columns in the image regardless of
  whether its a mosaic.
  Returns -1 if error.
  Author: Douglas N. Greve, 9/6/2001
  -----------------------------------------------------------------------*/
int dcmGetNCols(const char *dcmfile)
{
  DCM_ELEMENT *e;
  int NCols;

  e = GetElementFromFile(dcmfile, 0x28, 0x11);
  if (e == NULL) {
    return (-1);
  }

  NCols = *(e->d.us);

  FreeElementData(e);
  free(e);

  return (NCols);
}
/*-----------------------------------------------------------------------
  dcmImageDirCos - Gets the RAS direction cosines for the col and row of
  the image based on DICOM tag (20,37). Vcx is the x-component of the
  unit vector that points from the center of one voxel to the center of
  an adjacent voxel in the next higher column within the same row and
  slice.
  Returns 1 if error.
  Author: Douglas N. Greve, 9/10/2001
  -----------------------------------------------------------------------*/
int dcmImageDirCos(const char *dcmfile, float *Vcx, float *Vcy, float *Vcz, float *Vrx, float *Vry, float *Vrz)
{
  DCM_ELEMENT *e;
  char *s;
  unsigned int n;
  int nbs;
  float rms;

  /* Load the direction cosines - this is a string of the form:
     Vcx\Vcy\Vcz\Vrx\Vry\Vrz */
  e = GetElementFromFile(dcmfile, 0x20, 0x37);
  if (e == NULL) {
    return (1);
  }
  s = e->d.string;

  /* replace back slashes with spaces */
  nbs = 0;
  for (n = 0; n < strlen(s); n++) {
    if (s[n] == '\\') {
      s[n] = ' ';
      nbs++;
    }
  }

  if (nbs != 5) {
    return (1);
  }

  sscanf(s, "%f %f %f %f %f %f ", Vcx, Vcy, Vcz, Vrx, Vry, Vrz);

  /* Convert Vc from LPS to RAS and Normalize */
  rms = sqrt((*Vcx) * (*Vcx) + (*Vcy) * (*Vcy) + (*Vcz) * (*Vcz));
  (*Vcx) /= -rms;
  (*Vcy) /= -rms;
  (*Vcz) /= +rms;

  /* Convert Vr from LPS to RAS and Normalize */
  rms = sqrt((*Vrx) * (*Vrx) + (*Vry) * (*Vry) + (*Vrz) * (*Vrz));
  (*Vrx) /= -rms;
  (*Vry) /= -rms;
  (*Vrz) /= +rms;

  FreeElementData(e);
  free(e);

  // Not sure why I had this here to begin with.
  // This tag should indicate HFS (for example).
  // 0020 0020 CS REL Patient Orientation
  // e = GetElementFromFile(dcmfile, 0x20, 0x20);
  // if (e==NULL) return (1);
  // s = e->d.string;
  // FreeElementData(e);
  // free(e);

  return (0);
}
/*-----------------------------------------------------------------------
  dcmImagePosition - Gets the RAS position of the center of the CRS=0
  voxel based on DICOM tag (20,32). Note that only the z component
  is valid for mosaics. See also sdfiFixImagePosition().
  Returns 1 if error.
  Author: Douglas N. Greve, 9/10/2001
  -----------------------------------------------------------------------*/
int dcmImagePosition(const char *dcmfile, float *x, float *y, float *z)
{
  DCM_ELEMENT *e;
  char *s;
  unsigned int n;
  int nbs;

  /* Load the Image Position: this is a string of the form:
     x\y\z  */
  e = GetElementFromFile(dcmfile, 0x20, 0x32);
  if (e == NULL) {
    return (1);
  }
  s = e->d.string;

  /* replace back slashes with spaces */
  nbs = 0;
  for (n = 0; n < strlen(s); n++) {
    if (s[n] == '\\') {
      s[n] = ' ';
      nbs++;
    }
  }

  if (nbs != 2) {
    return (1);
  }

  sscanf(s, "%f %f %f ", x, y, z);

  /* Convert from LPS to RAS */
  (*x) *= -1.0;
  (*y) *= -1.0;

  FreeElementData(e);
  free(e);

  return (0);
}
/*-----------------------------------------------------------------------
  sdcmSliceDirCos - Gets the RAS direction cosines for the slice base on
  the Siemens ASCII header. In the ASCII header, there are components
  of the form sSliceArray.asSlice[0].sNormal.dXXX, where XXX is Sag, Cor,
  and/or Tra. These form the vector that is perpendicular to the
  slice plane in the direction of increasing slice number. If absent,
  the component is set to zero.

  Notes:
  1. Converts from Siemens/DICOM LIS to RAS.
  2. Normalizes the vector.

  Returns 1 if error.
  Author: Douglas N. Greve, 9/10/2001
  -----------------------------------------------------------------------*/
int sdcmSliceDirCos(const char *dcmfile, float *Vsx, float *Vsy, float *Vsz)
{
  char *tmpstr;
  float rms;

  if (!IsSiemensDICOM(dcmfile)) {
    return (1);
  }

  tmpstr = SiemensAsciiTagEx(dcmfile, "sSliceArray.asSlice[0].sNormal.dSag", 0);
  if (tmpstr != NULL) {
    sscanf(tmpstr, "%f", Vsx);
    free(tmpstr);
  }

  tmpstr = SiemensAsciiTagEx(dcmfile, "sSliceArray.asSlice[0].sNormal.dCor", 0);
  if (tmpstr != NULL) {
    sscanf(tmpstr, "%f", Vsy);
    free(tmpstr);
  }

  tmpstr = SiemensAsciiTagEx(dcmfile, "sSliceArray.asSlice[0].sNormal.dTra", 0);
  if (tmpstr != NULL) {
    sscanf(tmpstr, "%f", Vsz);
    free(tmpstr);
  }

  if (*Vsx == 0 && *Vsy == 0 && *Vsz == 0) {
    sliceDirCosPresent = 0;
    return (1);
  }

  /* Convert from LPS to RAS */
  (*Vsx) *= -1.0;
  (*Vsy) *= -1.0;

  /* Normalize */
  rms = sqrt((*Vsx) * (*Vsx) + (*Vsy) * (*Vsy) + (*Vsz) * (*Vsz));

  (*Vsx) /= rms;
  (*Vsy) /= rms;
  (*Vsz) /= rms;

  sliceDirCosPresent = 1;

  return (0);
}

/*-----------------------------------------------------------------------
  sdcmIsMosaic() - tests whether a siemens dicom file has a mosaic image.
  If it is a mosaic it returns 1 and puts the dimension of the volume
  into pNcols, pNrows, pNslices, pNframes. These pointer arguments can
  be NULL, in which case they are ignored.

  This function works by computing the expected number of rows and columns
  assuming that the image is not a mosaic.  This is done by dividing the
  Phase Encode FOV by the image resolution in the phase encode direction
  (and same with Readout FOV).

  Author: Douglas N. Greve, 9/6/2001
  -----------------------------------------------------------------------*/
int sdcmIsMosaic(const char *dcmfile, int *pNcols, int *pNrows, int *pNslices, int *pNframes)
{
  DCM_ELEMENT *e;
  char *PhEncDir;
  int Nrows, Ncols;
  float ColRes, RowRes, SliceRes;
  int NrowsExp, NcolsExp;
  int NimagesMosaic, NmosaicSideLen;
  float PhEncFOV, ReadOutFOV;
  int err, IsMosaic;
  char *tmpstr;

  if (!IsSiemensDICOM(dcmfile)) {
    return (0);
  }

  tmpstr = getenv("SDCM_ISMOSAIC_OVERRIDE");
  if (tmpstr != NULL) {
    sscanf(tmpstr, "%d", &IsMosaic);
    printf("Env Override: IsMosaic = %d\n", IsMosaic);
    return (IsMosaic);
  }

  IsMosaic = 0;

  /* Get the phase encode direction: should be COL or ROW */
  /* COL means that each row is a different phase encode (??)*/
  e = GetElementFromFile(dcmfile, 0x18, 0x1312);
  if (e == NULL) {
    return (0);
  }
  PhEncDir = deblank((char *)e->d.string);
  FreeElementData(e);
  free(e);

  Nrows = dcmGetNRows(dcmfile);
  if (Nrows == -1) {
    return (0);
  }

  Ncols = dcmGetNCols(dcmfile);
  if (Ncols == -1) {
    return (0);
  }

  /* 2019-10-05, mu40: try to derive dimensions from Siemens' private
   * NumberOfImagesInMosaic field first, which represents the number of slices
   * in the run. Note that mosaics are always square, i.e. filled with empty
   * slices at the end. */
  e = GetElementFromFile(dcmfile, 0x19, 0x100a);
  NimagesMosaic = 0;
  if (e != NULL) {
    IsMosaic = 1;
    NimagesMosaic = (int)*(e->d.us);
    NmosaicSideLen = ceil(sqrt(NimagesMosaic));
    NrowsExp = Nrows / NmosaicSideLen;
    NcolsExp = Ncols / NmosaicSideLen;
  }
  else {
    tmpstr = SiemensAsciiTagEx(dcmfile, "sSliceArray.asSlice[0].dPhaseFOV", 0);
    if (tmpstr == NULL) {
      return (0);
    }
    sscanf(tmpstr, "%f", &PhEncFOV);
    free(tmpstr);

    tmpstr = SiemensAsciiTagEx(dcmfile, "sSliceArray.asSlice[0].dReadoutFOV", 0);
    if (tmpstr == NULL) {
      return (0);
    }
    sscanf(tmpstr, "%f", &ReadOutFOV);
    free(tmpstr);

    err = dcmGetVolRes(dcmfile, &ColRes, &RowRes, &SliceRes);
    if (err) {
      return (-1);
    }

    if (strncmp(PhEncDir, "COL", 3) == 0) {
      /* Each row is a different phase encode */
      NrowsExp = (int)(rint(PhEncFOV / RowRes));
      NcolsExp = (int)(rint(ReadOutFOV / ColRes));
    }
    else {
      /* Each column is a different phase encode */
      NrowsExp = (int)(rint(ReadOutFOV / RowRes));
      NcolsExp = (int)(rint(PhEncFOV / ColRes));
    }
  }

  if (NrowsExp != Nrows || NcolsExp != Ncols) {
    IsMosaic = 1;
    if (pNrows != NULL) {
      tmpstr = getenv("NROWS_OVERRIDE");
      if (tmpstr == NULL) {
        *pNrows = NrowsExp;
      }
      else {
        sscanf(tmpstr, "%d", pNrows);
        printf("Overriding number of rows with %d\n", *pNrows);
      }
    }
    if (pNcols != NULL) {
      tmpstr = getenv("NCOLS_OVERRIDE");
      if (tmpstr == NULL) {
        *pNcols = NcolsExp;
      }
      else {
        sscanf(tmpstr, "%d", pNcols);
        printf("Overriding number of columns with %d\n", *pNcols);
      }
    }
    if (pNslices != NULL) {
      tmpstr = getenv("NSLICES_OVERRIDE"); // was NSLICES_OVERRIDE_BCHWAUNIE
      if (tmpstr == NULL && NimagesMosaic > 0) {
        *pNslices = NimagesMosaic;
      }
      else if (tmpstr == NULL) {
        tmpstr = SiemensAsciiTagEx(dcmfile, "sSliceArray.lSize", 0);
        if (tmpstr == NULL) {
          return (0);
        }
        sscanf(tmpstr, "%d", pNslices);
        free(tmpstr);
      }
      else {
        sscanf(tmpstr, "%d", pNslices);
        printf("Overriding number of slices with %d\n", *pNslices);
      }
    }
    if (pNframes != NULL) {
      tmpstr = SiemensAsciiTagEx(dcmfile, "lRepetitions", 0);
      if (tmpstr == NULL) {
        return (0);
      }
      sscanf(tmpstr, "%d", pNframes);
      (*pNframes)++;
      free(tmpstr);
    }
  }
  free(PhEncDir);

  return (IsMosaic);
}
/*----------------------------------------------------------------
  GetSDCMFileInfo() - this fills a SDCMFILEINFO structure for a
  single Siemens DICOM file. Some of the data are filled from
  the DICOM header and some from the Siemens ASCII header. The
  pixel data are not loaded.
  ----------------------------------------------------------------*/
SDCMFILEINFO *GetSDCMFileInfo(const char *dcmfile)
{
  DCM_OBJECT *object = 0;
  SDCMFILEINFO *sdcmfi;
  CONDITION cond;
  DCM_TAG tag;
  DCM_ELEMENT *e;
  int l;
  unsigned short ustmp = 0;
  double dtmp = 0;
  char *strtmp, *strtmp2, *pc;
  int retval, nDiffDirections, nB0;
  double xr, xa, xs, yr, ya, ys, zr, za, zs;
  int DoDWI;

  if (!IsSiemensDICOM(dcmfile)) {
    return (NULL);
  }

  sdcmfi = (SDCMFILEINFO *)calloc(1, sizeof(SDCMFILEINFO));

  fflush(stdout);
  fflush(stderr);
  object = GetObjectFromFile(dcmfile, 0);
  fflush(stdout);
  fflush(stderr);
  if (object == NULL) {
    exit(1);
  }

  l = strlen(dcmfile);
  sdcmfi->FileName = (char *)calloc(l + 1, sizeof(char));
  memmove(sdcmfi->FileName, dcmfile, l);

  // This stores the 'Transfer Syntax Unique Identification',
  // which reports the structure of the image data, revealing
  // whether the data has been compressed. See:
  // http://www.psychology.nottingham.ac.uk/staff/cr1/dicom.html
  tag = DCM_MAKETAG(0x2, 0x10);
  cond = GetString(&object, tag, &sdcmfi->TransferSyntaxUID);

  tag = DCM_MAKETAG(0x10, 0x10);
  cond = GetString(&object, tag, &sdcmfi->PatientName);

  tag = DCM_MAKETAG(0x8, 0x20);
  cond = GetString(&object, tag, &sdcmfi->StudyDate);

  tag = DCM_MAKETAG(0x8, 0x30);
  cond = GetString(&object, tag, &sdcmfi->StudyTime);

  tag = DCM_MAKETAG(0x8, 0x31);
  cond = GetString(&object, tag, &sdcmfi->SeriesTime);

  tag = DCM_MAKETAG(0x8, 0x32);
  cond = GetString(&object, tag, &sdcmfi->AcquisitionTime);

  tag = DCM_MAKETAG(0x8, 0x1090);
  cond = GetString(&object, tag, &sdcmfi->ScannerModel);

  tag = DCM_MAKETAG(0x18, 0x1020);
  cond = GetString(&object, tag, &sdcmfi->NumarisVer);

  tag = DCM_MAKETAG(0x18, 0x24);
  cond = GetString(&object, tag, &sdcmfi->PulseSequence);

  tag = DCM_MAKETAG(0x18, 0x1030);
  cond = GetString(&object, tag, &sdcmfi->ProtocolName);
  if (strlen(sdcmfi->ProtocolName) == 0) {
    sdcmfi->ProtocolName = strcpyalloc("PROTOCOL_UNKOWN");
  }

  tag = DCM_MAKETAG(0x20, 0x11);
  cond = GetUSFromString(&object, tag, &ustmp);
  if (cond != DCM_NORMAL) {
    printf("WARNING: No Series Number (20,11) found in %s\n", sdcmfi->FileName);
    sdcmfi->ErrorFlag = 1;
  }
  sdcmfi->SeriesNo = (int)ustmp;
  sdcmfi->RunNo = sdcmfi->SeriesNo - 1;

  tag = DCM_MAKETAG(0x20, 0x13);
  cond = GetUSFromString(&object, tag, &ustmp);
  if (cond != DCM_NORMAL) {
    printf("WARNING: No Image Number (20,13) found in %s\n", sdcmfi->FileName);
    sdcmfi->ErrorFlag = 1;
  }
  sdcmfi->ImageNo = (int)ustmp;

  sdcmfi->UseSliceScaleFactor = 0;
  sdcmfi->SliceScaleFactor = 1;
  if (getenv("FS_NO_SLICE_SCALE_FACTOR") == NULL) {
    tag = DCM_MAKETAG(0x20, 0x4000);
    cond = GetString(&object, tag, &strtmp);
    if (cond == DCM_NORMAL) {
      if (strncmp(strtmp, "Scale Factor", 12) == 0) {
        sscanf(strtmp, "%*s %*s %lf", &sdcmfi->SliceScaleFactor);
        // printf("Slice Scale Factor %lf\n",sdcmfi->SliceScaleFactor);
        sdcmfi->UseSliceScaleFactor = 1;
      }
    }
  }

  sdcmfi->RescaleIntercept = 0;
  tag = DCM_MAKETAG(0x28, 0x1052);
  cond = GetString(&object, tag, &strtmp);
  if(cond == DCM_NORMAL) {
    sscanf(strtmp, "%lf", &sdcmfi->RescaleIntercept);
    free(strtmp);
    if(sdcmfi->RescaleIntercept != 0.0 && Gdiag_no > 0)
      printf("Info: %d RescaleIntercept = %lf \n",DCM_ImageNumber,sdcmfi->RescaleIntercept);
  }
  sdcmfi->RescaleSlope = 1.0;
  tag = DCM_MAKETAG(0x28, 0x1053);
  cond = GetString(&object, tag, &strtmp);
  if(cond == DCM_NORMAL) {
    sscanf(strtmp, "%lf", &sdcmfi->RescaleSlope);
    free(strtmp);
    if(sdcmfi->RescaleSlope != 1.0 && Gdiag_no > 0)
      printf("Info: %d RescaleSlope = %lf \n",DCM_ImageNumber,sdcmfi->RescaleSlope);
  }
  //DCMcheckInterceptSlope(object);

  tag = DCM_MAKETAG(0x18, 0x86);
  cond = GetUSFromString(&object, tag, &ustmp);
  sdcmfi->EchoNo = (int)ustmp;

  tag = DCM_MAKETAG(0x18, 0x1020);
  cond = GetDoubleFromString(&object, tag, &dtmp);

  tag = DCM_MAKETAG(0x18, 0x1314);
  cond = GetDoubleFromString(&object, tag, &dtmp);
  if (cond == DCM_NORMAL) {
    sdcmfi->FlipAngle = (float)M_PI * dtmp / 180.0;
  }
  else {
    sdcmfi->FlipAngle = 0;
  }

  tag = DCM_MAKETAG(0x18, 0x81);
  cond = GetDoubleFromString(&object, tag, &dtmp);
  sdcmfi->EchoTime = (float)dtmp;

  tag = DCM_MAKETAG(0x18, 0x87);
  cond = GetDoubleFromString(&object, tag, &dtmp);
  sdcmfi->FieldStrength = (float)dtmp;

  tag = DCM_MAKETAG(0x18, 0x82);
  cond = GetDoubleFromString(&object, tag, &dtmp);
  if (cond == DCM_NORMAL)
    sdcmfi->InversionTime = (float)dtmp;
  else
    sdcmfi->InversionTime = -1;

  e = GetElementFromFile(dcmfile, 0x28, 0x107);
  if (e)
    sdcmfi->LargestValue = (float)*(e->d.us);
  else
    sdcmfi->LargestValue = 0;

  /* Get the phase encode direction: should be COL or ROW */
  /* COL means that each row is a different phase encode (??)*/
  // Phase Enc Dir (DNG)
  tag = DCM_MAKETAG(0x18, 0x1312);
  cond = GetString(&object, tag, &strtmp);
  sdcmfi->PhEncDir = deblank(strtmp);
  free(strtmp);

  tag = DCM_MAKETAG(0x18, 0x80);
  cond = GetDoubleFromString(&object, tag, &dtmp);
  sdcmfi->RepetitionTime = (float)dtmp;

  strtmp = SiemensAsciiTagEx(dcmfile, "lRepetitions", 0);
  if (strtmp != NULL) {
    // This can cause problems with DTI scans if lRepetitions is actually set
    sscanf(strtmp, "%d", &(sdcmfi->lRepetitions));
    free(strtmp);
  }
  else {
    strtmp = SiemensAsciiTag(dcmfile, "sDiffusion.lDiffDirections", 0);
    strtmp2 = SiemensAsciiTag(dcmfile, "sWiPMemBlock.alFree[8]", 0);
    if (strtmp != NULL && strtmp2 != NULL) {
      sscanf(strtmp, "%d", &nDiffDirections);
      sscanf(strtmp, "%d", &nB0);
      sdcmfi->lRepetitions = nB0 + nDiffDirections - 1;
      DTIparsePulseSeqName(sdcmfi->PulseSequence, &sdcmfi->bValue, &sdcmfi->nthDirection);
      // printf("nDiffDirections = %d, nB0 = %d\n",nDiffDirections,nB0);
      // printf("%s %g %d\n",sdcmfi->PulseSequence,
      //     sdcmfi->bValue, sdcmfi->nthDirection);
    }
    else {
      sdcmfi->lRepetitions = 0;
    }
    if (strtmp) {
      free(strtmp);
    }
    if (strtmp2) {
      free(strtmp2);
    }
  }
  sdcmfi->NFrames = sdcmfi->lRepetitions + 1;
  /* This is not the last word on NFrames. See sdfiAssignRunNo().*/

  strtmp = SiemensAsciiTagEx(dcmfile, "sSliceArray.lSize", 0);
  if (strtmp != NULL) {
    sscanf(strtmp, "%d", &(sdcmfi->SliceArraylSize));
    free(strtmp);
  }
  else {
    sdcmfi->SliceArraylSize = 0;
  }

  strtmp = SiemensAsciiTagEx(dcmfile, "sSliceArray.asSlice[0].dPhaseFOV", 0);
  if (strtmp != NULL) {
    sscanf(strtmp, "%f", &(sdcmfi->PhEncFOV));
    free(strtmp);
  }
  else {
    sdcmfi->PhEncFOV = 0;
  }

  strtmp = SiemensAsciiTagEx(dcmfile, "sSliceArray.asSlice[0].dReadoutFOV", 0);
  if (strtmp != NULL) {
    sscanf(strtmp, "%f", &(sdcmfi->ReadoutFOV));
    free(strtmp);
  }
  else {
    sdcmfi->ReadoutFOV = 0;
  }

  sdcmfi->NImageRows = dcmGetNRows(dcmfile);
  if (sdcmfi->NImageRows < 0) {
    printf("WARNING: Could not determine number of image rows in %s\n", sdcmfi->FileName);
    sdcmfi->ErrorFlag = 1;
  }
  sdcmfi->NImageCols = dcmGetNCols(dcmfile);
  if (sdcmfi->NImageCols < 0) {
    printf("WARNING: Could not determine number of image cols in %s\n", sdcmfi->FileName);
    sdcmfi->ErrorFlag = 1;
  }

  dcmImagePosition(dcmfile, &(sdcmfi->ImgPos[0]), &(sdcmfi->ImgPos[1]), &(sdcmfi->ImgPos[2]));

  dcmImageDirCos(dcmfile,
                 &(sdcmfi->Vc[0]),
                 &(sdcmfi->Vc[1]),
                 &(sdcmfi->Vc[2]),
                 &(sdcmfi->Vr[0]),
                 &(sdcmfi->Vr[1]),
                 &(sdcmfi->Vr[2]));

  /* The following may return 1 (Vs[i] = 0 for all i) when there is no
     ASCII header (anonymization?). This is a show-stopper for mosaics.
     For non-mosaics, it is recoverable because we can sort the files
     and compute the slice dir cos from the image position.*/
  retval = sdcmSliceDirCos(dcmfile, &(sdcmfi->Vs[0]), &(sdcmfi->Vs[1]), &(sdcmfi->Vs[2]));

  sdcmfi->IsMosaic = sdcmIsMosaic(dcmfile, NULL, NULL, NULL, NULL);

  /* If could not get sliceDirCos, then we calculate an initial value.
     This might not be used at all. If it is used, then it is only
     used to sort the files.  If this initial sliceDirCos is wrong,
     then it can only be wrong by a sign. This would cause the slices
     may be reversed relative to how they were acquired. The final
     geometry will be correct because the sliceDirCos is recomputed
     based on the final sorting. The only time this would create a
     problem is if someone thought that the first slice in the volume
     was the first slice acquired.
  */
  if (retval == 1 && sdcmfi->IsMosaic == 0) {
    /* we have x_(r,a,s) and y_(r,a,s).  z_(r,a,s) must be
       orthogonal to these */
    /* get the cross product of x_(r,a,s) and y_(r,a,s) is
       in proportion to z_(r,a,s) */
    /* also x_(r,a,s) and y_(r,a,s) are normalized and
       thus cross product is also normalized */
    xr = sdcmfi->Vc[0];
    xa = sdcmfi->Vc[1];
    xs = sdcmfi->Vc[2];
    yr = sdcmfi->Vr[0];
    ya = sdcmfi->Vr[1];
    ys = sdcmfi->Vr[2];
    zr = xa * ys - xs * ya;
    za = xs * yr - xr * ys;
    zs = xr * ya - xa * yr;
    // Assume left-handed
    sdcmfi->Vs[0] = -zr;
    sdcmfi->Vs[1] = -za;
    sdcmfi->Vs[2] = -zs;
    /* Confirm sign by two files later  */
  }

  dcmGetVolRes(dcmfile, &(sdcmfi->VolRes[0]), &(sdcmfi->VolRes[1]), &(sdcmfi->VolRes[2]));

  if (sdcmfi->IsMosaic) {
    sdcmIsMosaic(dcmfile, &(sdcmfi->VolDim[0]), &(sdcmfi->VolDim[1]), &(sdcmfi->VolDim[2]), &(sdcmfi->NFrames));
  }
  else {
    sdcmfi->VolDim[0] = sdcmfi->NImageCols;
    sdcmfi->VolDim[1] = sdcmfi->NImageRows;
  }

  DoDWI = 1;
  pc = getenv("FS_LOAD_DWI");
  if (pc == NULL)
    DoDWI = 1;
  else if (strcmp(pc, "0") == 0)
    DoDWI = 0;

  if (DoDWI) {
    double bval, xbvec, ybvec, zbvec;
    int err;
    err = dcmGetDWIParams(object, &bval, &xbvec, &ybvec, &zbvec);
    if (err) {
      printf("ERROR: GetSDCMFileInfo(): dcmGetDWIParams() %d\n", err);
      printf("DICOM File: %s\n", dcmfile);
      printf("break %s:%d\n", __FILE__, __LINE__);
      return (NULL);
    }
    if (Gdiag_no > 0)
      printf("GetSDCMFileInfo(): DWI: %s %d %lf %lf %lf %lf\n", dcmfile, err, bval, xbvec, ybvec, zbvec);
    sdcmfi->bval = bval;
    sdcmfi->bvecx = xbvec;
    sdcmfi->bvecy = ybvec;
    sdcmfi->bvecz = zbvec;
  }
  else {
    sdcmfi->bval = 0;
    sdcmfi->bvecx = 0;
    sdcmfi->bvecy = 0;
    sdcmfi->bvecz = 0;
  }

  // cleanup Ascii storage
  SiemensAsciiTagEx(dcmfile, (char *)0, 1);

  cond = DCM_CloseObject(&object);

  /* Clear the condition stack to prevent overflow */
  COND_PopCondition(1);

  return (sdcmfi);
}
/*----------------------------------------------------------*/
int DumpSDCMFileInfo(FILE *fp, SDCMFILEINFO *sdcmfi)
{
  fprintf(fp, "FileName \t\t%s\n", sdcmfi->FileName);
  fprintf(fp, "Identification\n");
  fprintf(fp, "\tNumarisVer        %s\n", sdcmfi->NumarisVer);
  fprintf(fp, "\tScannerModel      %s\n", sdcmfi->ScannerModel);
  fprintf(fp, "\tPatientName       %s\n", sdcmfi->PatientName);
  fprintf(fp, "Date and time\n");
  fprintf(fp, "\tStudyDate         %s\n", sdcmfi->StudyDate);
  fprintf(fp, "\tStudyTime         %s\n", sdcmfi->StudyTime);
  fprintf(fp, "\tSeriesTime        %s\n", sdcmfi->SeriesTime);
  fprintf(fp, "\tAcqTime           %s\n", sdcmfi->AcquisitionTime);
  fprintf(fp, "Acquisition parameters\n");
  fprintf(fp, "\tPulseSeq          %s\n", sdcmfi->PulseSequence);
  fprintf(fp, "\tProtocol          %s\n", sdcmfi->ProtocolName);
  fprintf(fp, "\tPhEncDir          %s\n", sdcmfi->PhEncDir);
  fprintf(fp, "\tEchoNo            %d\n", sdcmfi->EchoNo);
  fprintf(fp, "\tFlipAngle         %g\n", DEGREES(sdcmfi->FlipAngle));
  fprintf(fp, "\tEchoTime          %g\n", sdcmfi->EchoTime);
  fprintf(fp, "\tInversionTime     %g\n", sdcmfi->InversionTime);
  fprintf(fp, "\tRepetitionTime    %g\n", sdcmfi->RepetitionTime);
  fprintf(fp, "\tPhEncFOV          %g\n", sdcmfi->PhEncFOV);
  fprintf(fp, "\tReadoutFOV        %g\n", sdcmfi->ReadoutFOV);
  fprintf(fp, "Image information\n");
  fprintf(fp, "\tRunNo             %d\n", sdcmfi->RunNo);
  fprintf(fp, "\tSeriesNo          %d\n", sdcmfi->SeriesNo);
  fprintf(fp, "\tImageNo           %d\n", sdcmfi->ImageNo);
  fprintf(fp, "\tNImageRows        %d\n", sdcmfi->NImageRows);
  fprintf(fp, "\tNImageCols        %d\n", sdcmfi->NImageCols);
  fprintf(fp, "\tNFrames           %d\n", sdcmfi->NFrames);
  fprintf(fp, "\tSliceArraylSize   %d\n", sdcmfi->SliceArraylSize);
  fprintf(fp, "\tIsMosaic          %d\n", sdcmfi->IsMosaic);

  fprintf(fp, "\tImgPos            %8.4f %8.4f %8.4f \n", sdcmfi->ImgPos[0], sdcmfi->ImgPos[1], sdcmfi->ImgPos[2]);

  fprintf(fp, "\tVolRes            %8.4f %8.4f %8.4f \n", sdcmfi->VolRes[0], sdcmfi->VolRes[1], sdcmfi->VolRes[2]);
  fprintf(fp, "\tVolDim            %3d      %3d      %3d \n", sdcmfi->VolDim[0], sdcmfi->VolDim[1], sdcmfi->VolDim[2]);
  fprintf(fp, "\tVc                %8.4f %8.4f %8.4f \n", sdcmfi->Vc[0], sdcmfi->Vc[1], sdcmfi->Vc[2]);
  fprintf(fp, "\tVr                %8.4f %8.4f %8.4f \n", sdcmfi->Vr[0], sdcmfi->Vr[1], sdcmfi->Vr[2]);
  fprintf(fp, "\tVs                %8.4f %8.4f %8.4f \n", sdcmfi->Vs[0], sdcmfi->Vs[1], sdcmfi->Vs[2]);

  fprintf(
      fp, "\tVolCenter         %8.4f %8.4f %8.4f \n", sdcmfi->VolCenter[0], sdcmfi->VolCenter[1], sdcmfi->VolCenter[2]);

  if (sdcmfi->TransferSyntaxUID != NULL) {
    fprintf(fp, "\tTransferSyntaxUID %s\n", sdcmfi->TransferSyntaxUID);
  }

  return (0);
}

/*-----------------------------------------------------------------*/
int FreeSDCMFileInfo(SDCMFILEINFO **ppsdcmfi)
{
  SDCMFILEINFO *p;

  p = *ppsdcmfi;

  if (p->FileName != NULL) {
    free(p->FileName);
  }
  if (p->TransferSyntaxUID != NULL) {
    free(p->TransferSyntaxUID);
  }
  if (p->PatientName != NULL) {
    free(p->PatientName);
  }
  if (p->StudyDate != NULL) {
    free(p->StudyDate);
  }
  if (p->StudyTime != NULL) {
    free(p->StudyTime);
  }
  if (p->SeriesTime != NULL) {
    free(p->SeriesTime);
  }
  if (p->AcquisitionTime != NULL) {
    free(p->AcquisitionTime);
  }
  if (p->PulseSequence != NULL) {
    free(p->PulseSequence);
  }
  if (p->ProtocolName != NULL) {
    free(p->ProtocolName);
  }
  if (p->PhEncDir != NULL) {
    free(p->PhEncDir);
  }

  free(*ppsdcmfi);
  return (0);
}
/*-----------------------------------------------------------------
  sdcmExtractNumarisVer() - extacts the NUMARIS version string. The
  string can be something like "syngo MR 2002B 4VA21A", in which case
  only the last string is the actual version string. Or it can look
  like "4VA12B". In either case, the last string is the version
  string. The input is the string obtained from DICOM element
  (18,1020). The Major number is the first number in the string. The
  Minor number is the number after VA.  The MinorMinor number is the
  number the Minor.
  -------------------------------------------------------------------*/
char *sdcmExtractNumarisVer(const char *e_18_1020, int *Maj, int *Min, int *MinMin)
{
  int l, n, m;
  char *ver;

  l = strlen(e_18_1020);
  if (l < 6) {
    printf(
        "Cannot parse NUMARIS version string %s\n"
        "found in dicom tag 18,1020 (string len < 6)\n",
        e_18_1020);
    return (NULL);
  }

  /* dont copy blanks at the end */
  n = l - 1;
  while (n >= 0 && e_18_1020[n] == ' ') {
    n--;
  }
  if (n < 0) {
    printf(
        "Could not parse NUMARIS version string %s\n"
        "found in dicom tag 18,1020 (all blanks)\n",
        e_18_1020);
    return (NULL);
  }

  /* count length of the first non-blank string */
  m = 0;
  while (n >= 0 && e_18_1020[n] != ' ') {
    m++;
    n--;
  }
  n++;

  if (m != 6) {
    printf(
        "Could not parse NUMARIS version string %s\n"
        "found in dicom tag 18,1020 (len = %d != 6)\n",
        e_18_1020,
        m);
    return (NULL);
  }

  ver = (char *)calloc(sizeof(char), m + 1);
  strncpy(ver, &e_18_1020[n], m);

  /* Now determine major and minor release numbers */
  if (ver[1] != 'V' || !(ver[2] == 'A' || ver[2] == 'B')) {
    printf(
        "Could not parse NUMARIS version string %s\n"
        "found in dicom tag 18,1020 (VA or VB not in 2nd and 3rd)\n",
        e_18_1020);
    return (NULL);
  }

  *Maj = ver[0] - 48;
  *Min = ver[3] - 48;
  *MinMin = ver[4] - 48;
  // printf("Maj = %d, Min = %d, MinMin = %d\n",*Maj,*Min,*MinMin);

  return (ver);
}
/*--------------------------------------------------------------------
  ScanSiemensDCMDir() - similar to ScanDir but returns only files that
  are Siemens DICOM Files. It also returns a pointer to an array of
  SDCMFILEINFO structures.

  Author: Douglas Greve.
  Date: 09/10/2001
  *------------------------------------------------------------------*/
SDCMFILEINFO **ScanSiemensDCMDir(const char *PathName, int *NSDCMFiles)
{
  struct dirent **NameList;
  int i, pathlength;
  int NFiles;
  char tmpstr[1000];
  SDCMFILEINFO **sdcmfi_list;
  int pct, sumpct;
  FILE *fp;

  char *pname = (char *)calloc(strlen(PathName) + 1, sizeof(char));
  strcpy(pname, PathName);

  /* Remove all trailing forward slashes from pname */
  pathlength = strlen(pname);
  if (pname[pathlength - 1] == '/') {
    for (i = pathlength - 1; i >= 0; i--) {
      if (pname[i] == '/') {
        pname[i] = ' ';
        break;
      }
    }
  }
  pathlength = strlen(pname);

  /* select all directory entries, and sort them by name */
  NFiles = scandir(pname, &NameList, 0, alphasort);

  if (NFiles < 0) {
    fprintf(stderr, "WARNING: No files found in %s\n", pname);
    free(pname);
    return (NULL);
  }
  fprintf(stderr, "INFO: Found %d files in %s\n", NFiles, pname);

  /* Count the number of Siemens DICOM Files */
  fprintf(stderr, "INFO: counting Siemens Files\n");
  (*NSDCMFiles) = 0;
  for (i = 0; i < NFiles; i++) {
    sprintf(tmpstr, "%s/%s", pname, NameList[i]->d_name);
    if (IsSiemensDICOM(tmpstr)) {
      (*NSDCMFiles)++;
    }
  }

  fprintf(stderr, "INFO: found %d Siemens Files\n", *NSDCMFiles);

  if (*NSDCMFiles == 0) {
    free(pname);
    return (NULL);
  }

  sdcmfi_list = (SDCMFILEINFO **)calloc(*NSDCMFiles, sizeof(SDCMFILEINFO *));

  fprintf(stderr, "INFO: scanning info from Siemens Files\n");

  if (SDCMStatusFile != NULL) {
    fprintf(stderr, "INFO: status file is %s\n", SDCMStatusFile);
  }

  fprintf(stderr, "%2d ", 0);
  sumpct = 0;
  (*NSDCMFiles) = 0;
  for (i = 0; i < NFiles; i++) {
    // fprintf(stderr,"%4d ",i);
    pct = rint(100 * (i + 1) / NFiles) - sumpct;
    if (pct >= 2) {
      sumpct += pct;
      fprintf(stderr, "%3d ", sumpct);
      fflush(stderr);
      if (SDCMStatusFile != NULL) {
        fp = fopen(SDCMStatusFile, "w");
        if (fp != NULL) {
          fprintf(fp, "%3d\n", sumpct);
          fclose(fp);
        }
      }
    }

    sprintf(tmpstr, "%s/%s", pname, NameList[i]->d_name);
    if (IsSiemensDICOM(tmpstr)) {
      sdcmfi_list[*NSDCMFiles] = GetSDCMFileInfo(tmpstr);
      if (sdcmfi_list[*NSDCMFiles] == NULL) {
        return (NULL);
      }
      (*NSDCMFiles)++;
    }
  }
  fprintf(stderr, "\n");

  // free memory
  while (NFiles--) {
    free(NameList[NFiles]);
  }
  free(NameList);

  free(pname);

  return (sdcmfi_list);
}
/*--------------------------------------------------------------------
  LoadSiemensSeriesInfo() - loads header info from each of the nList
  files listed in SeriesList. This list is obtained from either
  ReadSiemensSeries() or ScanSiemensSeries().
  Author: Douglas Greve.
  Date: 09/10/2001
  *------------------------------------------------------------------*/
SDCMFILEINFO **LoadSiemensSeriesInfo(char **SeriesList, int nList)
{
  SDCMFILEINFO **sdfi_list;
  int n;

  // printf("LoadSiemensSeriesInfo()\n");

  sdfi_list = (SDCMFILEINFO **)calloc(nList, sizeof(SDCMFILEINFO *));

  for (n = 0; n < nList; n++) {
    fflush(stdout);
    fflush(stderr);

    // printf("%3d %s ---------------\n",n,SeriesList[n]);
    if (!IsSiemensDICOM(SeriesList[n])) {
      fprintf(stderr, "ERROR: %s is not a Siemens DICOM File\n", SeriesList[n]);
      fflush(stderr);
      free(sdfi_list);
      return (NULL);
    }

    // printf("Getting file info %s ---------------\n",SeriesList[n]);
    fflush(stdout);
    fflush(stderr);
    sdfi_list[n] = GetSDCMFileInfo(SeriesList[n]);
    if (sdfi_list[n] == NULL) {
      fprintf(stderr, "ERROR: reading %s \n", SeriesList[n]);
      fflush(stderr);
      free(sdfi_list);
      return (NULL);
    }
    exec_progress_callback(n, nList, 0, 1);
  }
  fprintf(stderr, "\n");
  fflush(stdout);
  fflush(stderr);

  return (sdfi_list);
}
/*--------------------------------------------------------------------
  ReadSiemensSeries() - returns a list of file names as found in the
  given ListFile. This file should contain a list of file names
  separated by white space. Any directory names are stripped off and
  replaced with the dirname of dcmfile. If dcmfile is not part of the
  list, it is added to the list.  No attempt is made to assure that
  the files are Siemens DICOM files or that they even exist. The
  resulting list is of the same form produced by ScanSiemensSeries();

  Author: Douglas Greve.
  Date: 09/25/2001
  *------------------------------------------------------------------*/
char **ReadSiemensSeries(const char *ListFile, int *nList, const char *dcmfile)
{
  FILE *fp;
  char **SeriesList;
  char tmpstr[1000];
  int n;
  char *dcmbase, *dcmdir, *fbase;
  int AddDCMFile;

  if (!IsSiemensDICOM(dcmfile)) {
    fprintf(stderr,
            "ERROR (ReadSiemensSeries): %s is not a "
            "Siemens DICOM file\n",
            dcmfile);
    return (NULL);
  }

  dcmdir = fio_dirname(dcmfile);
  dcmbase = fio_basename(dcmfile, NULL);

  fp = fopen(ListFile, "r");
  if (fp == NULL) {
    fprintf(stderr, "ERROR: could not open %s for reading\n", ListFile);
    return (NULL);
  }

  /* Count the numer of files in the list. Look for dcmfile and
     turn off the flag if it is found.*/
  AddDCMFile = 1;
  (*nList) = 0;
  while (fscanf(fp, "%s", tmpstr) != EOF) {
    fbase = fio_basename(tmpstr, NULL);
    if (strcmp(dcmbase, fbase) == 0) {
      AddDCMFile = 0;
    }
    (*nList)++;
    free(fbase);
  }

  if ((*nList) == 0) {
    fprintf(stderr, "ERROR: no files found in %s\n", ListFile);
    fflush(stderr);
    return (NULL);
  }
  fseek(fp, 0, SEEK_SET); /* go back to the begining */

  if (AddDCMFile) {
    fprintf(stderr, "INFO: adding dcmfile to list.\n");
    fflush(stderr);
    (*nList)++;
  }

  fprintf(stderr, "INFO: found %d files in list file\n", (*nList));
  fflush(stderr);

  SeriesList = (char **)calloc((*nList), sizeof(char *));

  for (n = 0; n < ((*nList) - AddDCMFile); n++) {
    if (fscanf(fp, "%s", tmpstr) != 1) {
      fprintf(stderr, "ERROR: could not scan string\n");
    }
    fbase = fio_basename(tmpstr, NULL);
    sprintf(tmpstr, "%s/%s", dcmdir, fbase);
    SeriesList[n] = (char *)calloc(strlen(tmpstr) + 1, sizeof(char));
    memmove(SeriesList[n], tmpstr, strlen(tmpstr));
    free(fbase);
  }
  if (AddDCMFile) {
    SeriesList[n] = (char *)calloc(strlen(dcmfile) + 1, sizeof(char));
    memmove(SeriesList[n], dcmfile, strlen(dcmfile));
  }

  fclose(fp);
  free(dcmdir);
  free(dcmbase);
  return (SeriesList);
}
/*--------------------------------------------------------------------
  ScanSiemensSeries() - scans a directory for Siemens DICOM files with
  the same Series Number as the given Siemens DICOM file. Returns a
  list of file names (including dcmfile), including the path. The
  resulting list is of the same form produced by ReadSiemensSeries();

  Author: Douglas Greve.
  Date: 09/25/2001
  *------------------------------------------------------------------*/
char **ScanSiemensSeries(const char *dcmfile, int *nList)
{
  int SeriesNo, SeriesNoTest;
  char *PathName;
  int NFiles, i;
  struct dirent **NameList;
  char **SeriesList;
  char tmpstr[1000];

  if (!IsSiemensDICOM(dcmfile)) {
    fprintf(stderr,
            "ERROR (ScanSiemensSeries): "
            "%s is not a Siemens DICOM file\n",
            dcmfile);
    fflush(stderr);
    return (NULL);
  }

  printf("Getting Series No \n");
  SeriesNo = dcmGetSeriesNo(dcmfile);
  if (SeriesNo == -1) {
    fprintf(stderr, "ERROR: reading series number from %s\n", dcmfile);
    fflush(stderr);
    return (NULL);
  }

  printf("Scanning Directory \n");
  /* select all directory entries, and sort them by name */
  PathName = fio_dirname(dcmfile);
  NFiles = scandir(PathName, &NameList, 0, alphasort);

  if (NFiles < 0) {
    fprintf(stderr, "WARNING: No files found in %s\n", PathName);
    fflush(stderr);
    return (NULL);
  }
  fprintf(stderr, "INFO: Found %d files in %s\n", NFiles, PathName);
  fprintf(stderr, "INFO: Scanning for Series Number %d\n", SeriesNo);
  fflush(stderr);

  /* Alloc enough memory for everyone */
  SeriesList = (char **)calloc(NFiles, sizeof(char *));
  (*nList) = 0;
  for (i = 0; i < NFiles; i++) {
    sprintf(tmpstr, "%s/%s", PathName, NameList[i]->d_name);
    // printf("Testing %s ----------------------------------\n",tmpstr);
    IsSiemensDICOM(tmpstr);
    if (IsSiemensDICOM(tmpstr)) {
      SeriesNoTest = dcmGetSeriesNo(tmpstr);
      if (SeriesNoTest == SeriesNo) {
        SeriesList[*nList] = (char *)calloc(strlen(tmpstr) + 1 + 8, sizeof(char));
        memmove(SeriesList[*nList], tmpstr, strlen(tmpstr));
        // printf("%3d  %s\n",*nList,SeriesList[*nList]);
        (*nList)++;
      }
    }
    exec_progress_callback(i, NFiles, 0, 1);
  }
  fprintf(stderr, "INFO: found %d files in series\n", *nList);
  fflush(stderr);

  // free memory
  while (NFiles--) {
    free(NameList[NFiles]);
  }
  free(NameList);

  if (*nList == 0) {
    free(SeriesList);
    return (NULL);
  }
  free(PathName);

  return (SeriesList);
}

/*-----------------------------------------------------------*/
int SortSDCMFileInfo(SDCMFILEINFO **sdcmfi_list, int nlist)
{
  qsort(sdcmfi_list, nlist, sizeof(SDCMFILEINFO **), CompareSDCMFileInfo);
  return (0);
}

/*-----------------------------------------------------------
  CompareSDCMFileInfo() - compares two siemens dicom files for
  the purposes of sorting them. They are sorted with qsort
  (via SortSDCMFileInfo()) which sorts them in ascending order.
  If CompareSDCMFileInfo() returns a -1, is means that the
  first file will appear ahead of the second file in the list,
  and vice versa if a +1 is returned.

  If ErrorFlag is set, then funny things will happen.

  Overall, files are sorted so that the run number always increases,
  regardless of whether a run is a mosaic or non-mosaic. The assignment
  of run number to a file is accomplished with sdfiAssignRunNo().

  Within a non-mosaic ran the first file will sort to that of the
  first frame of the first slice of the first run. If there are
  multiple frames in the run, then the next file will be that of the
  second frame of the first slice of the first run. When all the
  frames of the first slice have been exhausted, the next file will
  contain the first frame of the second slice, and so on.

  Within mosaic runs, the files are sorted in order of temporal
  acquision.

  The first comparison is based on Series Number. If the first file
  has a lower series number, a -1 is returned (ie, those with lower
  series numbers will appear ahead of those with higher numbers in the
  sorted list).  If they both have the same series number, the next
  comparison is evaluated.

  The second comparison is based on relative Slice Position. This is
  determined from the slice direction cosine (SDC) and the XYZ
  coordinates of the first voxel in each file. The SDC points in the
  direction of increasing slice number and is compared to the vector
  pointing from the first voxel of the first file to the first voxel
  of the second file. If they are parallel, then the first file has
  a lower slice number than the second, and it will appear earlier
  in the sorted list (ie, a -1 is returned). The dot product between
  the vector and the SDC is used to determine whether it parallel
  (dot = 1) or anti-parallel (dot = -1). To avoid floating point
  errors, the dot must be greater than 0.5 (parallel) or less than
  -0.5 (anti-parallel), otherwise the next comparison is evaluated.

  The third comparison is the Image Number, which indicates the temporal
  sequence (except in mosaics). Those that occur earlier in time (ie,
  have a lower image number) will appear earlier in the sorted list.

  The fourth comparison is the Echo Number.

  If the two files cannot be discriminated from these comparisions,
  a warning is printed and a 0 is returned.

  Notes and Assumptions:

  1. The series number always changes from one run to another.
  2. For mosaics, the series number increments with each frame.
  3. For mosaics, the image number is meaningless. Because of
  note 2, the image number comparison should never be reached.
  4. For non-mosaics, image number increases for later acquisitions.

  2019-10-05, mu40: only consider slice locations for non-mosaics, as this
  condition may be mistakenly triggered e.g. if prospective motion correction
  was being used. Assumptions 2 and 3 were observed to be violated in HCP
  rs-fMRI data.
  -----------------------------------------------------------*/
int CompareSDCMFileInfo(const void *a, const void *b)
{
  SDCMFILEINFO *sdcmfi1;
  SDCMFILEINFO *sdcmfi2;
  int n;
  float dv[3], dvsum2, dot;
  // float actm1, actm2;

  sdcmfi1 = *((SDCMFILEINFO **)a);
  sdcmfi2 = *((SDCMFILEINFO **)b);

  if (sdcmfi1->ErrorFlag) {
    return (-1);
  }
  if (sdcmfi2->ErrorFlag) {
    return (+1);
  }

  /* Sort by Series Number */
  if (sdcmfi1->SeriesNo < sdcmfi2->SeriesNo) {
    return (-1);
  }
  if (sdcmfi1->SeriesNo > sdcmfi2->SeriesNo) {
    return (+1);
  }

  /* ------ Sort by Slice Position -------- */
  if (!sdcmfi1->IsMosaic || !sdcmfi2->IsMosaic) {
    /* Compute vector from the first to the second */
    dvsum2 = 0;
    for (n = 0; n < 3; n++) {
      dv[n] = sdcmfi2->ImgPos[n] - sdcmfi1->ImgPos[n];
      dvsum2 += (dv[n] * dv[n]);
    }
    for (n = 0; n < 3; n++) {
      dv[n] /= sqrt(dvsum2); /* normalize */
    }
    /* Compute dot product with Slice Normal vector */
    dot = 0;
    for (n = 0; n < 3; n++) {
      dot += (dv[n] * sdcmfi1->Vs[n]);
    }
    // Sort by slice position.
    if (dot > +0.5) {
      return (-1);
    }
    if (dot < -0.5) {
      return (+1);
    }
  }

  /* Sort by Image Number (Temporal Sequence) */
  if (sdcmfi1->ImageNo < sdcmfi2->ImageNo) {
    return (-1);
  }
  if (sdcmfi1->ImageNo > sdcmfi2->ImageNo) {
    return (+1);
  }

  /* Sort by Echo Number  */
  if (sdcmfi1->EchoNo < sdcmfi2->EchoNo) {
    return (-1);
  }
  if (sdcmfi1->EchoNo > sdcmfi2->EchoNo) {
    return (+1);
  }

  /* Sort by Acquisition Time */
  /* This has been commented out because it should be
     redundant with ImageNo */
  // sscanf(sdcmfi1->AcquisitionTime,"%f",&actm1);
  // sscanf(sdcmfi2->AcquisitionTime,"%f",&actm2);
  // if(actm1 < actm2) return(-1);
  // if(actm1 > actm2) return(+1);

  printf(
      "WARNING: files are not found to be different and "
      "cannot be sorted\n");
  printf("File1: %s\n", sdcmfi1->FileName);
  printf("File2: %s\n", sdcmfi2->FileName);

  return (0);
}
/*-----------------------------------------------------------
  sdfiAssignRunNo2() - assigns run number based on series number
  -----------------------------------------------------------*/
int sdfiAssignRunNo2(SDCMFILEINFO **sdfi_list, int nlist)
{
  SDCMFILEINFO *sdfi, *sdfitmp, *sdfi0;
  int nthfile, NRuns, nthrun, nthslice, nthframe;
  int nfilesperrun = 0, firstpass, nframes;
  char *FirstFileName = 0;
  int *RunList = 0, *RunNoList = 0;

  nframes = 0; /* to stop compiler warnings */

  RunNoList = sdfiRunNoList(sdfi_list, nlist, &NRuns);
  if (NRuns == 0) {
    return (NRuns);
  }

  for (nthrun = 0; nthrun < NRuns; nthrun++) {
    FirstFileName = sdfiFirstFileInRun(RunNoList[nthrun], sdfi_list, nlist);
    RunList = sdfiRunFileList(FirstFileName, sdfi_list, nlist, &nfilesperrun);

    sdfi0 = sdfi_list[RunList[0]];

    if (sdfi0->IsMosaic) {
      /* 2019-10-05, mu40: do not set the error flag if the number of files in
       * the run exceeds the number of repetitions, as the run is most likely
       * not truncated. The lRepetition field may not have been set, e.g. in
       * vNavs. */
      sdfi0->NFrames = nfilesperrun;
      if (nfilesperrun != (sdfi0->lRepetitions + 1)) {
        fprintf(stderr, "WARNING: Run %d appears to be truncated\n", nthrun + 1);
        fprintf(stderr, "  Files Found: %d, Files Expected (lRep+1): %d\n", nfilesperrun, (sdfi0->lRepetitions + 1));
        DumpSDCMFileInfo(stderr, sdfi0);
        fflush(stderr);
        if (nfilesperrun < (sdfi0->lRepetitions + 1)) {
            sdfi0->ErrorFlag = 1;
        }
      }
    }

    else {
      /* It is NOT a mosaic */

      nthfile = 0;
      nthslice = 0;
      firstpass = 1;

      while (nthfile < nfilesperrun) {
        sdfi = sdfi_list[RunList[nthfile]];
        sdfitmp = sdfi_list[RunList[nthfile]];
        nthframe = 0;
        while (sdfiSameSlicePos(sdfi, sdfitmp)) {
          nthframe++;
          nthfile++;
          if (nthfile < nfilesperrun) {
            sdfitmp = sdfi_list[RunList[nthfile]];
          }
          else {
            break;
          }
        }
        if (firstpass) {
          firstpass = 0;
          nframes = nthframe;
        }
        if (nthframe != nframes) {
          fprintf(stderr, "WARNING: Run %d appears to be truncated\n", RunNoList[nthrun]);
          fprintf(stderr, "  Slice = %d, nthframe = %d, nframes = %d, %d\n", nthslice, nthframe, nframes, firstpass);
          fflush(stderr);
          sdfi0->ErrorFlag = 1;
          break;
        }
        nthslice++;
      } /* end loop over files in the run */

      sdfi0->VolDim[2] = nthslice;
      sdfi0->NFrames = nframes;
    } /* end if it is NOT a mosaic */

    /* Update the parameters for all files in the run */
    for (nthfile = 0; nthfile < nfilesperrun; nthfile++) {
      sdfi = sdfi_list[RunList[nthfile]];
      sdfi->VolDim[2] = sdfi0->VolDim[2];
      sdfi->NFrames = sdfi0->NFrames;
      sdfi->ErrorFlag = sdfi0->ErrorFlag;
    }

    if (RunList) {
      free(RunList);
      RunList = 0;
    }
    if (FirstFileName) {
      free(FirstFileName);
      FirstFileName = 0;
    }
  } /* end loop over runs */
  if (RunNoList) {
    free(RunNoList);
    RunNoList = 0;
  }
  return (NRuns);
}
/*-----------------------------------------------------------
  sdfiRunNoList() - returns a list of run numbers
  -----------------------------------------------------------*/
int *sdfiRunNoList(SDCMFILEINFO **sdfi_list, int nlist, int *NRuns)
{
  SDCMFILEINFO *sdfi;
  int nthfile, PrevRunNo;
  int *RunNoList;
  int nthrun;

  *NRuns = sdfiCountRuns(sdfi_list, nlist);
  if (*NRuns == 0) {
    return (NULL);
  }

  RunNoList = (int *)calloc(*NRuns, sizeof(int));

  nthrun = 0;
  PrevRunNo = -1;
  for (nthfile = 0; nthfile < nlist; nthfile++) {
    sdfi = sdfi_list[nthfile];
    if (PrevRunNo == sdfi->RunNo) {
      continue;
    }
    PrevRunNo = sdfi->RunNo;
    RunNoList[nthrun] = sdfi->RunNo;
    fprintf(stderr, "RunNo = %d\n", sdfi->RunNo);
    nthrun++;
  }
  return (RunNoList);
}
/*-----------------------------------------------------------
  sdfiCountRuns() - counts the number of runs in the list
  -----------------------------------------------------------*/
int sdfiCountRuns(SDCMFILEINFO **sdfi_list, int nlist)
{
  SDCMFILEINFO *sdfi;
  int nthfile, NRuns, PrevRunNo;

  NRuns = 0;
  PrevRunNo = -1;
  for (nthfile = 0; nthfile < nlist; nthfile++) {
    sdfi = sdfi_list[nthfile];
    if (PrevRunNo == sdfi->RunNo) {
      continue;
    }
    PrevRunNo = sdfi->RunNo;
    NRuns++;
  }
  return (NRuns);
}
/*-----------------------------------------------------------
  sdfiCountFilesInRun() - counts the number of files in a
  given run. This differs from sdfiNFilesInRun() in that
  this takes a run number instead of a file name.
  -----------------------------------------------------------*/
int sdfiCountFilesInRun(int RunNo, SDCMFILEINFO **sdfi_list, int nlist)
{
  SDCMFILEINFO *sdfi;
  int nthfile, NFilesInRun;

  NFilesInRun = 0;
  for (nthfile = 0; nthfile < nlist; nthfile++) {
    sdfi = sdfi_list[nthfile];
    if (sdfi->RunNo == RunNo) {
      NFilesInRun++;
    }
  }
  return (NFilesInRun);
}
/*-----------------------------------------------------------
  sdfiFirstFileInRun() - returns the name of the first file
  in the given run.
  -----------------------------------------------------------*/
char *sdfiFirstFileInRun(int RunNo, SDCMFILEINFO **sdfi_list, int nlist)
{
  SDCMFILEINFO *sdfi;
  int nthfile, len;
  char *FirstFileName;

  for (nthfile = 0; nthfile < nlist; nthfile++) {
    sdfi = sdfi_list[nthfile];
    if (sdfi->RunNo == RunNo) {
      len = strlen(sdfi->FileName);
      FirstFileName = (char *)calloc(len + 1, sizeof(char));
      memmove(FirstFileName, sdfi->FileName, len);
      return (FirstFileName);
    }
  }
  fprintf(stderr, "WARNING: could not find Run %d in list\n", RunNo);
  return (NULL);
}
/*-----------------------------------------------------------
  sdfiAssignRunNo() - assigns a run number to each file. It is assumed
  that the list has already been ordered by SortSDCMFileInfo(), which
  sorts (according to CompareSDCMFileInfo()) by Series Number, Slice
  Position, and Image Number.

  The run number that each file is associated with is determined
  incrementally, starting with the first file in the list.

  If the first file IS NOT A MOSAIC, then all the subsequent files
  in the list with the same series number are assigned to the
  same run as the first file. In this case, the total number of
  files in the run is equal to the number of slices times the number
  of frames. The number of frames is determined by the number of
  files with the same Image Position. The NFrames element of the
  structure for each file in the run is changed to reflect this.

  If the first file IS A MOSAIC, then the number of files in the run
  is assumed to equal the number of frames as indicated by
  NFrames = lRepeptitions + 1.

  The process is repeated with the first file of the second run, and
  so on, until all the files are accounted for.

  Returns:
  Number of Runs Found (no error)
  0 if there was an error

  Notes:
  1. List must have been sorted by SortSDCMFileInfo().
  2. For non-mosaics, the number of slices is determined here.
  2. For non-mosaics with multiple frames, the number of frames
  (ie, ->NFrames) is also determined here.
  -----------------------------------------------------------*/
int sdfiAssignRunNo(SDCMFILEINFO **sdcmfi_list, int nfiles)
{
  int nthfile, nthrun, nthslice, nthframe, serno, sernotest;
  int nfilesperrun;
  int ncols, nrows, nslices, nframes;
  SDCMFILEINFO *sdfi, *sdfitmp;
  char *FirstFile;
  int FirstFileNo, n, nthfileperrun;
  // char *tmpstr;

  // tmpstr = NULL;

  nthfile = 0;
  nthrun = 0;
  nthframe = 0;

#ifdef _DEBUG
  printf("    File    NthFl Ser Img  NFrs   Run NthFlRun\n");
#endif

  while (nthfile < nfiles) {
    nfilesperrun = 0;

    sdfi = sdcmfi_list[nthfile];
    FirstFile = sdfi->FileName;
    FirstFileNo = nthfile;

    if (!sdfi->IsMosaic) {
      /*--------------------------------------------------*/
      /* Its not a mosaic --- sort by series no and frame */
      ncols = sdfi->NImageCols;
      nrows = sdfi->NImageRows;
      serno = sdfi->SeriesNo;
      sernotest = serno;
      nthslice = 0;

      while (sernotest == serno) {
        sdfitmp = sdfi;
        nthframe = 0;

        while (sdfiSameSlicePos(sdfitmp, sdfi)) {
          sdfi->RunNo = nthrun;
#ifdef _DEBUG
          tmpstr = fio_basename(sdfi->FileName, NULL);
          printf("%10s %4d  %2d  %3d  %3d   %3d   %3d  (%g,%g,%g)\n",
                 tmpstr,
                 nthfile,
                 sdfi->SeriesNo,
                 sdfi->ImageNo,
                 sdfi->NFrames,
                 nthrun,
                 nfilesperrun + 1,
                 sdfi->ImgPos[0],
                 sdfi->ImgPos[1],
                 sdfi->ImgPos[2]);
          fflush(stdout);
          free(tmpstr);
#endif
          nthframe++;
          nthfile++;
          nfilesperrun++;
          if (nthfile < nfiles) {
            sdfi = sdcmfi_list[nthfile];
            sernotest = sdfi->SeriesNo;
          }
          else {
            sernotest = -1;
            break;
          }
        }
        nthslice++;
      }
      nslices = nthslice;
      nframes = nthframe;
      if (nfilesperrun != (nslices * nframes)) {
        printf("ERROR: not enough files for run %d\n", nthrun);
        printf("nslices = %d, nframes = %d\n", nslices, nframes);
        printf("nexpected = %d, nfound = %d\n", (nslices * nframes), nfilesperrun);
        printf("FirstFile: %s\n", FirstFile);
        return (0);
      }
    } /* if(mosaic) */
    else {
      /*--------------------------------------------------*/
      /* It is a mosaic -- count out the number of frames */
      nfilesperrun = sdfi->NFrames;
      nframes = sdfi->NFrames; /* this is for compat with non-mos */
      nthfileperrun = 0;
      sdcmIsMosaic(FirstFile, &ncols, &nrows, &nslices, NULL);
      for (nthframe = 0; nthframe < sdfi->NFrames; nthframe++) {
        if (nthfile >= nfiles) {
          fprintf(stdout, "ERROR: not enough files for %d frames of run %d\n", sdfi->NFrames, nthrun);
          fprintf(stdout, "%s\n", FirstFile);
          fflush(stdout);
          return (0);
        }
        sdfi = sdcmfi_list[nthfile];
        sdfi->RunNo = nthrun;

#ifdef _DEBUG
        tmpstr = fio_basename(sdfi->FileName, NULL);
        printf("%10s %4d  %2d  %3d  %3d   %3d   %3d \n",
               tmpstr,
               nthfile,
               sdfi->SeriesNo,
               sdfi->ImageNo,
               sdfi->NFrames,
               sdfi->RunNo,
               nthfileperrun + 1);
        fflush(stdout);
        free(tmpstr);
#endif

        nthfileperrun++;
        nthfile++;
      }
    } /* if(not mosaic) */

    /* Make sure each fileinfo in the run has a complete
       set of volume dimension information */
    for (n = FirstFileNo; n < FirstFileNo + nfilesperrun; n++) {
      sdfi = sdcmfi_list[n];
      sdfi->VolDim[0] = ncols;
      sdfi->VolDim[1] = nrows;
      sdfi->VolDim[2] = nslices;
      sdfi->NFrames = nframes; /* has no effect for mosaic */
    }

#ifdef _DEBUG
    tmpstr = fio_basename(FirstFile, NULL);
    printf("%2d %10s %3d   (%3d %3d %3d %3d)\n", nthrun, tmpstr, nfilesperrun, ncols, nrows, nslices, sdfi->NFrames);
    free(tmpstr);
#endif

    nthrun++;
  }

  return (nthrun);
}
/*------------------------------------------------------------------
  sdfiRunNo() - this returns the Run Number of the given file name.
  The list is searched until the FileName member matches that of
  dcmfile. The RunNo member is then returned. The sdfi_list must
  have been sorted into runs with sdfiAssignRunNo().
  ------------------------------------------------------------------*/
int sdfiRunNo(const char *dcmfile, SDCMFILEINFO **sdfi_list, int nlist)
{
  int nthfile;
  SDCMFILEINFO *sdfi;

  for (nthfile = 0; nthfile < nlist; nthfile++) {
    sdfi = sdfi_list[nthfile];
    if (strcmp(dcmfile, sdfi->FileName) == 0) {
      return (sdfi->RunNo);
    }
  }

  fprintf(stderr, "WARNING: file %s not found in list\n", dcmfile);

  return (-1);
}
/*------------------------------------------------------------------
  sdfiNFilesInRun() - this returns the number of files in the run
  with dcmfile. The Run Number for dcmfile is obtained using sdfiRunNo().
  The list is searched for all the files with the same Run Number.
  ------------------------------------------------------------------*/
int sdfiNFilesInRun(const char *dcmfile, SDCMFILEINFO **sdfi_list, int nlist)
{
  int nthfile;
  int RunNo;
  int NFilesInRun;
  SDCMFILEINFO *sdfi;

  RunNo = sdfiRunNo(dcmfile, sdfi_list, nlist);
  if (RunNo < 0) {
    return (1);
  }

  NFilesInRun = 0;
  for (nthfile = 0; nthfile < nlist; nthfile++) {
    sdfi = sdfi_list[nthfile];
    if (sdfi->RunNo == RunNo) {
      NFilesInRun++;
    }
  }

  return (NFilesInRun);
}
/*------------------------------------------------------------------
  sdfiRunFileList() - returns a list of indices from sdfi_list that
  share the same Run Number. The number in the list is passed as
  NRunList.
  ------------------------------------------------------------------*/
int *sdfiRunFileList(const char *dcmfile, SDCMFILEINFO **sdfi_list, int nlist, int *NRunList)
{
  int nthfile, nthfileinrun;
  SDCMFILEINFO *sdfi;
  int *RunList;
  int RunNo;

  /* get the run number for this dcmfile */
  RunNo = sdfiRunNo(dcmfile, sdfi_list, nlist);
  if (RunNo < 0) {
    return (NULL);
  }

  /* get the number of files in this run */
  *NRunList = sdfiNFilesInRun(dcmfile, sdfi_list, nlist);

  /* alloc the list of file numbers */
  RunList = (int *)calloc(*NRunList, sizeof(int));

  nthfileinrun = 0;
  for (nthfile = 0; nthfile < nlist; nthfile++) {
    sdfi = sdfi_list[nthfile];
    if (sdfi->RunNo == RunNo) {
      if (nthfileinrun >= *NRunList) {
        free(RunList);
        fprintf(stderr,
                "ERROR: found more than %d files "
                "in the run for DICOM file %s\n",
                *NRunList,
                sdfi->FileName);
        return (NULL);
      }
      RunList[nthfileinrun] = nthfile;
      // printf("%3d  %3d\n",nthfileinrun,RunList[nthfileinrun]);
      nthfileinrun++;
    }
  }
  return (RunList);
}
/*-----------------------------------------------------------------------
  sdfiFixImagePosition() - Fixes the (RAS) Image Position of a Siemens
  mosaic.  In Siemens DICOM files, the Image Position (20,32) is
  incorrect for mosaics. The Image Position is supposed to be the XYZ
  location of the first voxel. It is actually what would be the
  location of the first voxel if a single image the size of the mosaic
  had been collected centered at the first slice in the mosaic. The
  correction is done by using info from the Siemens ASCII header which
  indicates the center of each slice and the direction cosine info to
  compute the location of the corner of the first slice.  This routine
  used to do the correction by computing the center of the first slice
  (which coincides with the center of the mosaic), and then computing
  the actual XYZ of the first slice based on its center, in-plane
  resolution, in-plane dimension, and direction cosines.  But this
  resulted in image positions that were slightly off from those derived
  the from non-mosaic volumes with identical slice prescriptions. This
  change was made on 7/28/05.
  -----------------------------------------------------------------------*/
int sdfiFixImagePosition(SDCMFILEINFO *sdfi)
{
  char *strtmp, *dcmfile;
  MATRIX *ras_c, *R, *crs_c, *ras0, *shift;
  int r;

  if (!sdfi->IsMosaic) {
    return (0);
  }

  R = MatrixAlloc(3, 3, MATRIX_REAL);
  for (r = 0; r < 3; r++) {
    R->rptr[r + 1][1] = sdfi->Vc[r] * sdfi->VolRes[0];
    R->rptr[r + 1][2] = sdfi->Vr[r] * sdfi->VolRes[1];
    R->rptr[r + 1][3] = sdfi->Vs[r] * sdfi->VolRes[2];
  }

  /* 2019-10-05, mu40: add alternative method for fixing mosaic position without
   * using the ASCII header. This is based on the previous routine that Doug
   * replaced in 2005. */
  if (getenv("FS_MOSAIC_FIX_NOASCII")) {
    printf("INFO: fixing mosaic center without using ASCII header\n");
    shift = MatrixAlloc(3, 1, MATRIX_REAL);
    shift->rptr[1][1] = (sdfi->NImageCols - sdfi->VolDim[0]) / 2.0;
    shift->rptr[2][1] = (sdfi->NImageRows - sdfi->VolDim[1]) / 2.0;
    shift->rptr[3][1] = 0;
    MatrixMultiplyD(R, shift, shift);
    sdfi->ImgPos[0] += shift->rptr[1][1];
    sdfi->ImgPos[1] += shift->rptr[2][1];
    sdfi->ImgPos[2] += shift->rptr[3][1];
    MatrixFree(&shift);
  }
  else {
    // Center of first slice
    ras_c = MatrixAlloc(3, 1, MATRIX_REAL);
    crs_c = MatrixAlloc(3, 1, MATRIX_REAL);

    dcmfile = sdfi->FileName;
    strtmp = SiemensAsciiTag(dcmfile, "sSliceArray.asSlice[0].sPosition.dSag", 0);
    if (strtmp != NULL) {
      sscanf(strtmp, "%f", &(ras_c->rptr[1][1]));
      ras_c->rptr[1][1] *= -1.0;
      free(strtmp);
    }
    strtmp = SiemensAsciiTag(dcmfile, "sSliceArray.asSlice[0].sPosition.dCor", 0);
    if (strtmp != NULL) {
      sscanf(strtmp, "%f", &(ras_c->rptr[2][1]));
      ras_c->rptr[2][1] *= -1.0;
      free(strtmp);
    }
    strtmp = SiemensAsciiTag(dcmfile, "sSliceArray.asSlice[0].sPosition.dTra", 0);
    if (strtmp != NULL) {
      sscanf(strtmp, "%f", &(ras_c->rptr[3][1]));
      free(strtmp);
    }

    crs_c->rptr[1][1] = sdfi->VolDim[0] / 2.0;
    crs_c->rptr[2][1] = sdfi->VolDim[1] / 2.0;
    crs_c->rptr[3][1] = 0;  // first slice

    ras0 = MatrixMultiply(R, crs_c, NULL);
    ras0 = MatrixSubtract(ras_c, ras0, ras0);

    sdfi->ImgPos[0] = ras0->rptr[1][1];
    sdfi->ImgPos[1] = ras0->rptr[2][1];
    sdfi->ImgPos[2] = ras0->rptr[3][1];

    MatrixFree(&ras_c);
    MatrixFree(&crs_c);
    MatrixFree(&ras0);
  }

  MatrixFree(&R);
  return (0);
}
/*-----------------------------------------------------------------------
  sdfiVolCenter() - Computes the RAS XYZ "center" of the volume for
  the given the first Siemens DICOM file of a run. For mosaics, this
  assumes that the Image Position has been fixed (see
  sdfiFixImagePosition()). Center is in quotes because it is not
  really the center. Rather it is the RAS XYZ position at voxel CRS =
  [Nc Nr Ns]/2. The true center would be at CRS = ([Nc Nr Ns]-1)/2.
  This definition of the "center" is consistent with the c_r, c_a,
  c_s in the MRI struct.

  Author: Douglas N. Greve, 9/12/2001. Updated 12/19/02.
  -----------------------------------------------------------------------*/
int sdfiVolCenter(SDCMFILEINFO *sdfi)
{
  int r, c;
  float Mdc[3][3], FoV[3];

  for (r = 0; r < 3; r++) {
    Mdc[r][0] = sdfi->Vc[r];
    Mdc[r][1] = sdfi->Vr[r];
    Mdc[r][2] = sdfi->Vs[r];
  }

  /* Note: for true center use (sdfi->VolDim[r]-1) */
  for (r = 0; r < 3; r++) {
    FoV[r] = sdfi->VolRes[r] * sdfi->VolDim[r];
  }

  /* sdfi->Volcenter[] is actually the RAS value of the corner voxel */
  /* (0,0,z).  Thus, if you use the first slice, then the value is   */
  /* the same as the translation part                                */
  /* Thus, you just add the Mdc*fov + translation gives c_(r,a,s)    */
  for (r = 0; r < 3; r++) {
    sdfi->VolCenter[r] = sdfi->ImgPos[r];
    for (c = 0; c < 3; c++) {
      sdfi->VolCenter[r] += Mdc[r][c] * FoV[c] / 2.0;
    }
  }
  return (0);
}

/*-----------------------------------------------------------------------
  sdfiSameSlicePos() - Determines whether the images in the two files
  have the same slice position.
  Author: Douglas N. Greve, 9/12/2001
  -----------------------------------------------------------------------*/
int sdfiSameSlicePos(SDCMFILEINFO *sdfi1, SDCMFILEINFO *sdfi2)
{
  static float eps = 0;
  if (eps == 0) {
    char *pc;
    pc = getenv("FS_SAME_SLICE_THRESH");
    if (pc == NULL)
      eps = .000001;
    else
      sscanf(pc, "%f", &eps);
    printf("sdfiSameSlicePos() eps = %f\n", eps);
  }

  if (fabs(sdfi1->ImgPos[0] - sdfi2->ImgPos[0]) > eps) {
    return (0);
  }
  if (fabs(sdfi1->ImgPos[1] - sdfi2->ImgPos[1]) > eps) {
    return (0);
  }
  if (fabs(sdfi1->ImgPos[2] - sdfi2->ImgPos[2]) > eps) {
    return (0);
  }
  return (1);
}

/*-----------------------------------------------------------------------
  sdfiIsSliceOrderReversed() - determine whether the slice order was
  reversed when packed into a mosiac. Slice reversals in non-mosaics
  are handled properly in the sorting routines. The slice order may be
  reversed due to setting "ImageNumbering" on the scanner console.
  This is sometimes set so that the images are displayed in a certain
  order (eg, head-to-foot) instead of how they were collected (eg,
  ascending or descending). This makes the slice normal direction
  cosine information inconsistent with the actual slice direction.
  Note that this does not have anything to do with the order in which
  the slices were acquired. However, if you assume a certain
  acquistion slice order model, it will be wrong unless this reversal
  is done.

  Author: Douglas N. Greve, 7/28/05
  -----------------------------------------------------------------------*/
int sdfiIsSliceOrderReversed(SDCMFILEINFO *sdfi)
{
  int correv, sagrev, trarev;
  char *dcmfile, *strtmp;

  // Slice reversals in non-mosaics are handled properly
  // in the file sorting routines because each slice is
  // in a different file.
  if (!sdfi->IsMosaic) {
    return (0);
  }
  dcmfile = sdfi->FileName;

  sagrev = 0;
  correv = 0;
  trarev = 0;

  /* These tags are set to 0x1 when a reversal has occured in the
     relevant dimension (sort of), otherwise they will be empty.
     Note that more than one may be set, but only applies to
     the slices. Eg, if Sag and Cor are set but the slice direction
     is primarily Cor, then only the slice order is reversed
     (ie, it is not reversed in Sag, which would be a left-right
     flip).
  */
  strtmp = SiemensAsciiTag(dcmfile, "sSliceArray.ucImageNumbSag", 0);
  if (strtmp != NULL) {
    sagrev = 1;
    free(strtmp);
  }
  strtmp = SiemensAsciiTag(dcmfile, "sSliceArray.ucImageNumbCor", 0);
  if (strtmp != NULL) {
    correv = 1;
    free(strtmp);
  }
  strtmp = SiemensAsciiTag(dcmfile, "sSliceArray.ucImageNumbTra", 0);
  if (strtmp != NULL) {
    trarev = 1;
    free(strtmp);
  }
  printf("sagrev = %d, correv =%d, trarev = %d\n", sagrev, correv, trarev);
  printf("Vs = %g %g %g\n", sdfi->Vs[0], sdfi->Vs[1], sdfi->Vs[2]);

  // No reversal if none of the ImageNumb flags are set
  if (!sagrev && !correv && !trarev) {
    return (0);
  }

  // Slices are primarily Sag and there is a sag reversal
  if ((fabs(sdfi->Vs[0]) > fabs(sdfi->Vs[1])) && (fabs(sdfi->Vs[0]) > fabs(sdfi->Vs[2])) && sagrev) {
    return (1);
  }

  // Slices are primarily Cor and there is a cor reversal
  if ((fabs(sdfi->Vs[1]) > fabs(sdfi->Vs[0])) && (fabs(sdfi->Vs[1]) > fabs(sdfi->Vs[2])) && correv) {
    return (1);
  }

  // Slices are primarily Axial and there is an axial reversal
  if ((fabs(sdfi->Vs[2]) > fabs(sdfi->Vs[0])) && (fabs(sdfi->Vs[2]) > fabs(sdfi->Vs[1])) && trarev) {
    return (1);
  }

  // If it gets here, there is no reversal
  return (0);
}

/*--------------------------------------------------------------*/
/* Sebastien's Routines are below */
/*--------------------------------------------------------------*/

/*******************************************************
   PrintDICOMInfo
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: structure DICOMInfo
   output: prints specific DICOM fields
*******************************************************/
void PrintDICOMInfo(DICOMInfo *dcminfo)
{
  int i;
  char str[256];

  printf("-------------------------------------------------\n");
  printf("DICOM meta-header\n\n");

  printf("file name\t\t\t%s\n", dcminfo->FileName);

  printf("Date and time\n");
  if (IsTagPresent[DCM_StudyDate]) {
    sprintf(str, "%s", dcminfo->StudyDate);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tstudy date\t\t%s\n", str);

  if (IsTagPresent[DCM_StudyTime]) {
    sprintf(str, "%s", dcminfo->StudyTime);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tstudy time\t\t%s\n", str);

  if (IsTagPresent[DCM_SeriesTime]) {
    sprintf(str, "%s", dcminfo->SeriesTime);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tseries time\t\t%s\n", str);

  if (IsTagPresent[DCM_AcquisitionTime]) {
    sprintf(str, "%s", dcminfo->AcquisitionTime);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tacquisition time\t%s\n", str);

  printf("Identification\n");
  if (IsTagPresent[DCM_PatientName]) {
    sprintf(str, "%s", dcminfo->PatientName);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tpatient name\t\t%s\n", str);

  if (IsTagPresent[DCM_Manufacturer]) {
    sprintf(str, "%s", dcminfo->Manufacturer);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tmanufacturer\t\t%s\n", str);

  printf("Dimensions\n");
  if (IsTagPresent[DCM_Rows]) {
    sprintf(str, "%d", dcminfo->Rows);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tnumber of rows\t\t%s\n", str);

  if (IsTagPresent[DCM_Columns]) {
    sprintf(str, "%d", dcminfo->Columns);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tnumber of columns\t%s\n", str);

  printf("\tnumber of frames\t%d\n", dcminfo->NumberOfFrames);

  if (IsTagPresent[DCM_xsize]) {
    sprintf(str, "%g", dcminfo->xsize);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tpixel width\t\t%s\n", str);

  if (IsTagPresent[DCM_ysize]) {
    sprintf(str, "%g", dcminfo->ysize);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tpixel height\t\t%s\n", str);

  if (IsTagPresent[DCM_SliceThickness]) {
    sprintf(str, "%g", dcminfo->SliceThickness);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tslice thickness\t\t%s\n", str);

  printf("\tfield of view\t\t%g\n", dcminfo->FieldOfView);

  if (IsTagPresent[DCM_ImageNumber]) {
    sprintf(str, "%d", dcminfo->ImageNumber);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\timage number\t\t%s (might be not reliable)\n", str);

  if (IsTagPresent[DCM_TransferSyntaxUID]) {
    printf("\ttransfer syntax UID\t%s\n", dcminfo->TransferSyntaxUID);
  }

  printf("Acquisition parameters\n");
  if (IsTagPresent[DCM_EchoTime]) {
    sprintf(str, "%g", dcminfo->EchoTime);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\techo time\t\t%s\n", str);

  if (IsTagPresent[DCM_RepetitionTime]) {
    sprintf(str, "%g", dcminfo->RepetitionTime);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\trepetition time\t\t%s\n", str);

  if (IsTagPresent[DCM_InversionTime]) {
    sprintf(str, "%g", dcminfo->InversionTime);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tinversion time\t\t%s\n", str);

  if (IsTagPresent[DCM_EchoNumber]) {
    sprintf(str, "%d", dcminfo->EchoNumber);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\techo number\t\t%s\n", str);

  if (IsTagPresent[DCM_FlipAngle]) {
    sprintf(str, "%g", dcminfo->FlipAngle);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tflip angle\t\t%s\n", str);

  if (IsTagPresent[DCM_BitsAllocated]) {
    sprintf(str, "%d", dcminfo->BitsAllocated);
  }
  else {
    strcpy(str, "not found");
  }
  printf("\tbits allocated\t\t%s\n", str);

  printf("Spatial information\n");

  printf("\tfirst image position\t");
  if (IsTagPresent[DCM_ImagePosition]) {
    for (i = 0; i < 3; i++) {
      printf("%g ", dcminfo->FirstImagePosition[i]);
    }
    printf("\n");
  }
  else {
    printf("notfound\n");
  }

  printf("\tlast image position\t");
  if (IsTagPresent[DCM_ImagePosition]) {
    for (i = 0; i < 3; i++) {
      printf("%g ", dcminfo->LastImagePosition[i]);
    }
    printf("\n");
  }
  else {
    printf("notfound\n");
  }

  printf("\timage orientation\t");
  if (IsTagPresent[DCM_ImageOrientation]) {
    for (i = 0; i < 6; i++) {
      printf("%g ", dcminfo->ImageOrientation[i]);
    }
    printf("\n");
  }
  else {
    printf("notfound\n");
  }

  printf("-------------------------------------------------\n\n");
}

CONDITION GetString(DCM_OBJECT **object, DCM_TAG tag, char **st)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void *ctx;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    *st = (char *)calloc(10, sizeof(char));
    memmove(*st, "unknown", 8);
    return cond;
  }
  *st = (char *)calloc(attribute.length + 9, sizeof(char));
  attribute.d.string = *st;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  return cond;
}

CONDITION GetUSFromString(DCM_OBJECT **object, DCM_TAG tag, unsigned short *us)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char *s;
  void *ctx;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  s = (char *)calloc(attribute.length + 1, sizeof(char));
  attribute.d.string = s;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *us = (unsigned short)atoi(s);
  free(s);
  return cond;
}

CONDITION GetShortFromString(DCM_OBJECT **object, DCM_TAG tag, short *sh)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char *s;
  void *ctx;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  s = (char *)calloc(attribute.length + 1, sizeof(char));
  attribute.d.string = s;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *sh = (short)atoi(s);
  free(s);
  return cond;
}

CONDITION GetUSFromUS(DCM_OBJECT **object, DCM_TAG tag, unsigned short *us)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void *ctx, *ot;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  ot = (void *)malloc(attribute.length + 1);
  attribute.d.ot = ot;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *us = *attribute.d.us;
  free(ot);
  return cond;
}

CONDITION GetShortFromShort(DCM_OBJECT **object, DCM_TAG tag, short *ss)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void *ctx, *ot;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  ot = (void *)malloc(attribute.length + 1);
  attribute.d.ot = ot;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *ss = *attribute.d.ss;
  free(ot);
  return cond;
}

CONDITION GetPixelData_Save(DCM_OBJECT **object, DCM_TAG tag, unsigned short **ad)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void *ctx, *ot;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  ot = (void *)malloc(attribute.length + 1);
  attribute.d.ot = ot;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *ad = (unsigned short *)ot;
  return cond;
}

CONDITION GetPixelData(DCM_OBJECT **object, DCM_TAG tag, void **ad)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  void *ctx, *ot;
  unsigned int pixelDataLength;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }

  pixelDataLength = attribute.length;
  if (pixelDataLength == 0xFFFFFFFF) {
    printf("ERROR: Invalid DICOM pixel data length = %8.8Xh\n", pixelDataLength);
    printf("       DICOM group:element 7FE0:0010\n");
    exit(1);
  }

  ot = (void *)malloc(pixelDataLength + 1);
  attribute.d.ot = ot;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *ad = ot;
  return cond;
}

CONDITION GetDoubleFromString(DCM_OBJECT **object, DCM_TAG tag, double *d)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char *s;
  void *ctx;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  s = (char *)calloc(attribute.length + 1, sizeof(char));
  attribute.d.string = s;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);
  *d = atof(s);
  free(s);
  return cond;
}

CONDITION GetMultiDoubleFromString(DCM_OBJECT **object, DCM_TAG tag, double *d[], int multiplicity)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char *s, *ss;
  void *ctx;
  unsigned int i;
  int j, mult;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  s = (char *)calloc(attribute.length + 1, sizeof(char));
  ss = (char *)calloc(attribute.length + 1, sizeof(char));
  attribute.d.string = s;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);

  i = 0;
  j = 0;
  mult = 0;
  while (mult < multiplicity && i < attribute.length) {
    j = 0;
    while (s[i] != '\\' && i < attribute.length) {
      ss[j++] = s[i++];
    }
    i++;
    ss[j] = '\0';
    (*d)[mult++] = atof(ss);
  }
  free(s);
  free(ss);
  return cond;
}

CONDITION GetMultiShortFromString(DCM_OBJECT **object, DCM_TAG tag, short *us[], int multiplicity)
{
  DCM_ELEMENT attribute;
  CONDITION cond;
  char *s, *ss;
  void *ctx;
  unsigned int i;
  int j, mult;

  attribute.tag = tag;
  cond = DCM_GetElement(object, tag, &attribute);
  if (cond != DCM_NORMAL) {
    return cond;
  }
  s = (char *)calloc(attribute.length + 1, sizeof(char));
  ss = (char *)calloc(attribute.length + 1, sizeof(char));
  attribute.d.string = s;
  ctx = NULL;
  cond = DCM_GetElementValue(object, &attribute, &attribute.length, &ctx);

  i = 0;
  j = 0;
  mult = 0;
  while (mult < multiplicity && i < attribute.length) {
    j = 0;
    while (s[i] != '\\' && i < attribute.length) {
      ss[j++] = s[i++];
    }
    i++;
    ss[j] = '\0';
    (*us)[mult++] = atoi(ss);
  }
  free(s);
  free(ss);
  return cond;
}

/*******************************************************
   GetDICOMInfo
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: file name, structure DICOMInfo, boolean, image number
   output: fills in structure DICOMInfo with DICOM meta-header fields
*******************************************************/

CONDITION GetDICOMInfo(const char *fname, DICOMInfo *dcminfo, BOOL ReadImage, int ImageNumber)
{
  DCM_OBJECT **object = (DCM_OBJECT **)calloc(1, sizeof(DCM_OBJECT *));
  DCM_TAG tag;
  CONDITION cond, cond2 = DCM_NORMAL;
  double *tmp = (double *)calloc(10, sizeof(double));
  short *itmp = (short *)calloc(3, sizeof(short));
  int i;
  char *strtmp = NULL, *pc;
  int DoDWI;

  // Transfer Syntax UIDs
  // see http://www.psychology.nottingham.ac.uk/staff/cr1/dicom.html
  // for a discussion on this element
  char uid_buf[64];

  cond = DCM_OpenFile(fname, DCM_PART10FILE | DCM_ACCEPTVRMISMATCH, object);
  if (cond != DCM_NORMAL) cond = DCM_OpenFile(fname, DCM_ORDERLITTLEENDIAN | DCM_ACCEPTVRMISMATCH, object);
  if (cond != DCM_NORMAL) cond = DCM_OpenFile(fname, DCM_ORDERBIGENDIAN | DCM_ACCEPTVRMISMATCH, object);
  if (cond != DCM_NORMAL) cond = DCM_OpenFile(fname, DCM_FORMATCONVERSION | DCM_ACCEPTVRMISMATCH, object);
  if (cond != DCM_NORMAL) {
    printf("not DICOM images, sorry...\n");
    exit(1);
  }

  dcminfo->FileName = (char *)calloc(strlen(fname) + 1, sizeof(char));
  strcpy(dcminfo->FileName, fname);
  dcminfo->FileName[strlen(fname)] = '\0';

  // transfer syntax UID (data format identifier, ie, raw or JPEG compressed)
  // It can be converted to non-jpeg with
  // setenv DCMDICTPATH /usr/pubsw/packages/dcmtk/current/share/dcmtk/dicom.dic
  // dcmdjpeg +te jpgdicom newdicom
  tag = DCM_MAKETAG(0x2, 0x10);
  cond = GetString(object, tag, &dcminfo->TransferSyntaxUID);
  if (cond != DCM_NORMAL) {
    dcminfo->TransferSyntaxUID = (char *)calloc(256, sizeof(char));
    strcpy(dcminfo->TransferSyntaxUID, "UNKNOWN");
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag TransferSystaxUID not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_TransferSyntaxUID] = true;
  }

  // JPEG compressed data is *not* supported by freesurfer.
  // It can be converted to non-jpeg with
  // setenv DCMDICTPATH /usr/pubsw/packages/dcmtk/current/share/dcmtk/dicom.dic
  // dcmdjpeg +te jpgdicom newdicom
  strncpy(uid_buf, dcminfo->TransferSyntaxUID, sizeof(uid_buf)-1);
  uid_buf[strlen(jpegCompressed_UID)] = 0;
  if (strcmp(uid_buf, jpegCompressed_UID) == 0) {
    // Don't do anything about this until pixel data are loaded
    // printf("WARNING: JPEG compressed image data not supported!\n");
    // printf("         (Transfer Syntax UID: %s)\n",dcminfo->TransferSyntaxUID);
    // Do not exit here because we may want only to read the header and
    // not the pixel data.
    // exit(1);
  }
  // RLL encoded data is *not* supported by freesurfer
  strncpy(uid_buf, dcminfo->TransferSyntaxUID, sizeof(uid_buf)-1);
  uid_buf[strlen(rllEncoded_UID)] = 0;
  if (strcmp(uid_buf, rllEncoded_UID) == 0) {
    printf("ERROR: RLL-encoded image data not supported!\n");
    printf("       (Transfer Syntax UID: %s)\n", dcminfo->TransferSyntaxUID);
    exit(1);
  }

  // manufacturer
  tag = DCM_MAKETAG(0x8, 0x70);
  cond = GetString(object, tag, &dcminfo->Manufacturer);
  if (cond != DCM_NORMAL) {
    dcminfo->Manufacturer = (char *)calloc(256, sizeof(char));
    strcpy(dcminfo->Manufacturer, "UNKNOWN");
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag Manufacturer not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_Manufacturer] = true;
  }

  // patient name
  tag = DCM_MAKETAG(0x10, 0x10);
  cond = GetString(object, tag, &dcminfo->PatientName);
  if (cond != DCM_NORMAL) {
    dcminfo->PatientName = (char *)calloc(256, sizeof(char));
    strcpy(dcminfo->PatientName, "UNKNOWN");
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag PatientName not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_PatientName] = true;
  }

  // study date
  tag = DCM_MAKETAG(0x8, 0x20);
  cond = GetString(object, tag, &dcminfo->StudyDate);
  if (cond != DCM_NORMAL || strcmp(dcminfo->StudyDate, "00000000") == 0) {
    dcminfo->StudyDate = (char *)calloc(256, sizeof(char));
    strcpy(dcminfo->StudyDate, "UNKNOWN");
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag StudyDate not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_StudyDate] = true;
  }

  // study time
  tag = DCM_MAKETAG(0x8, 0x30);
  cond = GetString(object, tag, &dcminfo->StudyTime);
  if (cond != DCM_NORMAL || strcmp(dcminfo->StudyTime, "00000000") == 0) {
    dcminfo->StudyTime = (char *)calloc(256, sizeof(char));
    strcpy(dcminfo->StudyTime, "UNKNOWN");
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag StudyTime not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_StudyTime] = true;
  }

  // series time
  tag = DCM_MAKETAG(0x8, 0x31);
  cond = GetString(object, tag, &dcminfo->SeriesTime);
  if (cond != DCM_NORMAL || strcmp(dcminfo->SeriesTime, "00000000") == 0) {
    dcminfo->SeriesTime = (char *)calloc(256, sizeof(char));
    strcpy(dcminfo->SeriesTime, "UNKNOWN");
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag SeriesTime not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_SeriesTime] = true;
  }

  // acquisition time
  tag = DCM_MAKETAG(0x8, 0x31);
  cond = GetString(object, tag, &dcminfo->AcquisitionTime);
  if (cond != DCM_NORMAL || strcmp(dcminfo->AcquisitionTime, "00000000") == 0) {
    dcminfo->AcquisitionTime = (char *)calloc(256, sizeof(char));
    strcpy(dcminfo->AcquisitionTime, "UNKNOWN");
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag AcquisitionTime not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_AcquisitionTime] = true;
  }

  // slice thickness (use Slice Spacing)
  tag = DCM_MAKETAG(0x18, 0x88);
  cond = GetDoubleFromString(object, tag, &dcminfo->SliceThickness);
  if (cond != DCM_NORMAL || dcminfo->SliceThickness == 0.0) {
#ifdef _DEBUG
    printf("WARNING: tag Slice Spacing not found in %s\n", fname);
#endif
    // try Slice Thickness tag
    tag = DCM_MAKETAG(0x18, 0x50);
    cond = GetDoubleFromString(object, tag, &dcminfo->SliceThickness);
    if (cond != DCM_NORMAL || dcminfo->SliceThickness == 0.0) {
      dcminfo->SliceThickness = 0.0;
      cond2 = cond;
#ifdef _DEBUG
      printf("WARNING: tag Slice Thickness not found\n");
#endif
    }
    else {
      IsTagPresent[DCM_SliceThickness] = true;
    }
  }
  else {
    IsTagPresent[DCM_SliceThickness] = true;
  }

  // image number
  tag = DCM_MAKETAG(0x20, 0x13);
  cond = GetUSFromString(object, tag, &dcminfo->ImageNumber);
  if ((cond != DCM_NORMAL || dcminfo->ImageNumber == 0) && ImageNumber != -1) {
    dcminfo->ImageNumber = ImageNumber;
    cond2 = cond;
    printf("WARNING: tag ImageNumber not found in %s\n", fname);
  }
  else {
    IsTagPresent[DCM_ImageNumber] = true;
  }

  // series number
  tag = DCM_MAKETAG(0x20, 0x11);
  cond = GetUSFromString(object, tag, &dcminfo->SeriesNumber);
  if (cond != DCM_NORMAL || dcminfo->SeriesNumber == 0) {
    cond2 = cond;
    printf("WARNING: tag SeriesNumber not found in %s\n", fname);
  }
  else {
    IsTagPresent[DCM_SeriesNumber] = true;
  }

  // rows
  tag = DCM_MAKETAG(0x28, 0x10);
  cond = GetUSFromUS(object, tag, &dcminfo->Rows);
  if (cond != DCM_NORMAL) {
    dcminfo->Rows = 0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag Rows not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_Rows] = true;
  }

  // columns
  tag = DCM_MAKETAG(0x28, 0x11);
  cond = GetUSFromUS(object, tag, &dcminfo->Columns);
  if (cond != DCM_NORMAL) {
    dcminfo->Columns = 0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag Columns not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_Columns] = true;
  }

  // pixel spacing
  tag = DCM_MAKETAG(0x28, 0x30);
  cond = GetMultiDoubleFromString(object, tag, &tmp, 2);
  if (cond != DCM_NORMAL || tmp[0] == 0.0 || tmp[1] == 0.0) {
    dcminfo->xsize = 0.0;
    dcminfo->ysize = 0.0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag Pixel spacing not found in %s\n", fname);
#endif
  }
  else {
    dcminfo->xsize = tmp[0];
    dcminfo->ysize = tmp[1];
    IsTagPresent[DCM_xsize] = true;
    IsTagPresent[DCM_ysize] = true;
  }

  // bits allocated
  tag = DCM_MAKETAG(0x28, 0x100);
  cond = GetUSFromUS(object, tag, &dcminfo->BitsAllocated);
  if (cond != DCM_NORMAL) {
    dcminfo->BitsAllocated = 0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag BitsAllocated not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_BitsAllocated] = true;
  }

  // repetition time
  tag = DCM_MAKETAG(0x18, 0x80);
  cond = GetDoubleFromString(object, tag, &dcminfo->RepetitionTime);
  if (cond != DCM_NORMAL) {
    dcminfo->RepetitionTime = 0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag RepetitionTime not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_RepetitionTime] = true;
  }

  // echo time
  tag = DCM_MAKETAG(0x18, 0x81);
  cond = GetDoubleFromString(object, tag, &dcminfo->EchoTime);
  if (cond != DCM_NORMAL) {
    dcminfo->EchoTime = 0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag EchoTime not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_EchoTime] = true;
  }

  tag = DCM_MAKETAG(0x18, 0x87);
  cond = GetDoubleFromString(object, tag, &dcminfo->FieldStrength);
  if (cond != DCM_NORMAL) {
    dcminfo->FieldStrength = 0;
    cond2 = cond;
  }

  // flip angle
  tag = DCM_MAKETAG(0x18, 0x1314);
  cond = GetDoubleFromString(object, tag, &dcminfo->FlipAngle);
  if (cond != DCM_NORMAL) {
    dcminfo->FlipAngle = 0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag FlipAngle not found in %s\n", fname);
#endif
  }
  else {
    dcminfo->FlipAngle = M_PI * dcminfo->FlipAngle / 180.0;
    IsTagPresent[DCM_FlipAngle] = true;
  }

  // inversion time
  tag = DCM_MAKETAG(0x18, 0x82);
  cond = GetDoubleFromString(object, tag, &dcminfo->InversionTime);
  if (cond != DCM_NORMAL) {
    dcminfo->InversionTime = 0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag InversionTime not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_InversionTime] = true;
  }

  // echo number
  tag = DCM_MAKETAG(0x18, 0x86);
  cond = GetShortFromString(object, tag, &dcminfo->EchoNumber);
  if (cond != DCM_NORMAL) {
    dcminfo->EchoNumber = 0;
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag EchoNumber not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_EchoNumber] = true;
  }

  // image position
  // watch out: DICOM gives the image position in LPS coordinates
  //            We are interested in RAS coordinates (-1,-1,1)!
  tag = DCM_MAKETAG(0x20, 0x32);
  cond = GetMultiDoubleFromString(object, tag, &tmp, 3);
  if (cond != DCM_NORMAL) {
    for (i = 0; i < 3; i++) {
      dcminfo->ImagePosition[i] = 0;
    }
    cond2 = cond;
#ifdef _DEBUG
    printf("WARNING: tag image position not found in %s\n", fname);
#endif
  }
  else {
    IsTagPresent[DCM_ImagePosition] = true;
    for (i = 0; i < 3; i++) {
      dcminfo->ImagePosition[i] = tmp[i];
    }
  }

  // image orientation
  tag = DCM_MAKETAG(0x20, 0x37);
  cond = GetMultiDoubleFromString(object, tag, &tmp, 6);
  if (cond != DCM_NORMAL) {
    for (i = 0; i < 6; i++) {
      dcminfo->ImageOrientation[i] = 0;
    }
    cond2 = cond;
    printf("WARNING: tag image orientation not found in %s\n", fname);
  }
  else {
    IsTagPresent[DCM_ImageOrientation] = true;
    for (i = 0; i < 6; i++) {
      dcminfo->ImageOrientation[i] = tmp[i];
    }
    dcminfo->Vc[0] = dcminfo->ImageOrientation[0];
    dcminfo->Vc[1] = dcminfo->ImageOrientation[1];
    dcminfo->Vc[2] = dcminfo->ImageOrientation[2];
    dcminfo->Vr[0] = dcminfo->ImageOrientation[3];
    dcminfo->Vr[1] = dcminfo->ImageOrientation[4];
    dcminfo->Vr[2] = dcminfo->ImageOrientation[5];
    // Set slice direction to 0 for now. Do not try to set this based
    // on cross product of Vc and Vr as the sign is still ambiguous.
    // Vs will be computed latter (by DCMSliceDir) based on the
    // direction from the first slice/file in the series to the 2nd
    // slice.
    dcminfo->Vs[0] = 0;
    dcminfo->Vs[1] = 0;
    dcminfo->Vs[2] = 0;
  }

  // Phase Enc Dir.
  tag = DCM_MAKETAG(0x18, 0x1312);
  cond = GetString(object, tag, &strtmp);
  if (cond != DCM_NORMAL) {
    dcminfo->PhEncDir = NULL;
  }
  else {
    dcminfo->PhEncDir = deblank(strtmp);
    free(strtmp);
  }

  dcminfo->RescaleIntercept = 0;
  //DCMcheckInterceptSlope(*object);
  tag = DCM_MAKETAG(0x28, 0x1052);
  cond = GetString(object, tag, &strtmp);
  if(cond == DCM_NORMAL) {
    sscanf(strtmp, "%lf", &dcminfo->RescaleIntercept);
    free(strtmp);
    if(dcminfo->RescaleIntercept != 0.0 && Gdiag_no > 0)
      printf("Info: %d RescaleIntercept = %lf \n",DCM_ImageNumber,dcminfo->RescaleIntercept);
  }
  dcminfo->RescaleSlope = 1.0;
  tag = DCM_MAKETAG(0x28, 0x1053);
  cond = GetString(object, tag, &strtmp);
  if(cond == DCM_NORMAL) {
    sscanf(strtmp, "%lf", &dcminfo->RescaleSlope);
    free(strtmp);
    if(dcminfo->RescaleSlope != 1.0 && Gdiag_no > 0)
      printf("Info: %d RescaleSlope = %lf \n",DCM_ImageNumber,dcminfo->RescaleSlope);
  }

  DoDWI = 1;
  pc = getenv("FS_LOAD_DWI");
  if (pc == NULL)
    DoDWI = 1;
  else if (strcmp(pc, "0") == 0)
    DoDWI = 0;

  if (DoDWI) {
    double bval, xbvec, ybvec, zbvec;
    int err;
    err = dcmGetDWIParams(*object, &bval, &xbvec, &ybvec, &zbvec);
    if (err) {
      printf("ERROR: GetDICOMInfo(): dcmGetDWIParams() %d\n", err);
      printf("DICOM File: %s\n", fname);
      printf("break %s:%d\n", __FILE__, __LINE__);
      return ((CONDITION)1);
    }
    if (Gdiag_no > 0) printf("GetDICOMInfo(): DWI: %s %d %lf %lf %lf %lf\n", fname, err, bval, xbvec, ybvec, zbvec);
    dcminfo->bval = bval;
    dcminfo->bvecx = xbvec;
    dcminfo->bvecy = ybvec;
    dcminfo->bvecz = zbvec;
  }
  else {
    dcminfo->bval = 0;
    dcminfo->bvecx = 0;
    dcminfo->bvecy = 0;
    dcminfo->bvecz = 0;
  }

  // pixel data
  if (ReadImage) {
    tag = DCM_MAKETAG(0x7FE0, 0x10);
    cond = GetPixelData(object, tag, &dcminfo->PixelData);
    if (cond != DCM_NORMAL) {
      dcminfo->PixelData = NULL;
      cond2 = cond;
#ifdef _DEBUG
      printf("WARNING: tag pixel data not found in %s\n", fname);
#endif
    }
  }

  cond = DCM_CloseObject(object);
  if (cond != DCM_NORMAL) {
    cond2 = cond;
  }

  dcminfo->FieldOfView = dcminfo->xsize * dcminfo->Rows;

  /* Clear the condition stack to prevent overflow */
  COND_PopCondition(1);

  free(tmp);
  free(itmp);
  return cond2;
}

/*******************************************************
   ReadDICOMImage
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: array of file names, number of elements
   output: pixels image (8 or 16 bits), and DICOM
           informations stored in first array element
*******************************************************/

void *ReadDICOMImage(int nfiles, DICOMInfo **aDicomInfo)
{
  int n, i, j, bitsAllocated, numberOfFrames, offset, nvox;
  DICOMInfo dcminfo;
  unsigned char *PixelData8, *v8 = NULL;
  unsigned short *PixelData16, *v16 = NULL;
  int npix;
  double firstPosition[3], lastPosition[3];
  unsigned short int min16 = 65535, max16 = 0;
  unsigned short int min8 = 255, max8 = 0;

  memset(firstPosition, 0, sizeof(firstPosition));
  memset(lastPosition, 0, sizeof(lastPosition));

  npix = aDicomInfo[0]->Columns * aDicomInfo[0]->Rows;
  numberOfFrames = nfiles;
  nvox = npix * numberOfFrames;
  aDicomInfo[0]->NumberOfFrames = numberOfFrames;
  bitsAllocated = aDicomInfo[0]->BitsAllocated;

  printf("reading DICOM image...\n");

  switch (bitsAllocated) {
    case 8:
      v8 = (unsigned char *)calloc(nvox, sizeof(unsigned char));
      if (v8 == NULL) {
        printf("ReadDICOMImage: allocation problem\n");
        exit(1);
      }
      break;
    case 16:
      v16 = (unsigned short *)calloc(nvox, sizeof(unsigned short));
      if (v16 == NULL) {
        printf("ReadDICOMImage: allocation problem\n");
        exit(1);
      }
      break;
  }  // switch

  for (n = 0; n < nfiles; n++) {
    GetDICOMInfo(aDicomInfo[n]->FileName, &dcminfo, TRUE, n);

    if (n == 0)
      for (i = 0; i < 3; i++) {
        aDicomInfo[0]->FirstImagePosition[i] = dcminfo.ImagePosition[i];
      }
    if (n == nfiles - 1)
      for (i = 0; i < 3; i++) {
        aDicomInfo[0]->LastImagePosition[i] = dcminfo.ImagePosition[i];
      }

    offset = npix * n;

    switch (dcminfo.BitsAllocated) {
      case 8:
        PixelData8 = (unsigned char *)dcminfo.PixelData;
        for (j = 0; j < npix; j++) {
          v8[offset + j] = PixelData8[j];
          if (PixelData8[j] > max8) {
            max8 = PixelData8[j];
          }
          else if (PixelData8[j] < min8) {
            min8 = PixelData8[j];
          }
        }
        aDicomInfo[0]->min8 = min8;
        aDicomInfo[0]->max8 = max8;
        free(PixelData8);
        break;

      case 16:
        PixelData16 = (unsigned short *)dcminfo.PixelData;
        for (j = 0; j < npix; j++) {
          v16[offset + j] = PixelData16[j];
          if (PixelData16[j] > max16) {
            max16 = PixelData16[j];
          }
          else if (PixelData16[j] < min16) {
            min16 = PixelData16[j];
          }
        }
        aDicomInfo[0]->min16 = min16;
        aDicomInfo[0]->max16 = max16;
        free(PixelData16);
        break;
    }
    exec_progress_callback(n, nfiles, 0, 1);
  }

  for (i = 0; i < 3; i++) {
    aDicomInfo[0]->ImagePosition[i] = (firstPosition[i] + lastPosition[i]) / 2.;
  }

  switch (bitsAllocated) {
    case 8:
      return ((void *)v8);
      break;
    case 16:
      return ((void *)v16);
      break;
    default:
      return NULL;
  }
}

// ReadDICOMImage2
//
///////////////////////////////////////////////////////////////////
void *ReadDICOMImage2(int nfiles, DICOMInfo **aDicomInfo, int startIndex)
{
  int n, i, j, bitsAllocated, numberOfFrames, offset, nvox;
  DICOMInfo dcminfo;
  unsigned char *PixelData8, *v8 = NULL;
  unsigned short *PixelData16, *v16 = NULL;
  int npix;
  double firstPosition[3], lastPosition[3];
  unsigned short int min16 = 65535, max16 = 0;
  unsigned short int min8 = 255, max8 = 0;

  memset(firstPosition, 0, sizeof(firstPosition));
  memset(lastPosition, 0, sizeof(lastPosition));

  npix = aDicomInfo[startIndex]->Columns * aDicomInfo[startIndex]->Rows;
  numberOfFrames = nfiles;
  nvox = npix * numberOfFrames;
  // filling missing info
  aDicomInfo[startIndex]->NumberOfFrames = numberOfFrames;
  bitsAllocated = aDicomInfo[startIndex]->BitsAllocated;

  printf("reading DICOM image...\n");

  switch (bitsAllocated) {
    case 8:
      v8 = (unsigned char *)calloc(nvox, sizeof(unsigned char));
      if (v8 == NULL) {
        printf("ReadDICOMImage: allocation problem\n");
        exit(1);
      }
      break;
    case 16:
      v16 = (unsigned short *)calloc(nvox, sizeof(unsigned short));
      if (v16 == NULL) {
        printf("ReadDICOMImage: allocation problem\n");
        exit(1);
      }
      break;
  }  // switch

  for (n = 0; n < nfiles; n++) {
    GetDICOMInfo(aDicomInfo[n + startIndex]->FileName, &dcminfo, TRUE, n);

    if (n == 0)
      // fill missing info
      for (i = 0; i < 3; i++) aDicomInfo[startIndex]->FirstImagePosition[i] = dcminfo.ImagePosition[i];
    if (n == nfiles - 1)
      // fill missing info
      for (i = 0; i < 3; i++) aDicomInfo[startIndex]->LastImagePosition[i] = dcminfo.ImagePosition[i];

    offset = npix * n;

    switch (dcminfo.BitsAllocated) {
      case 8:
        PixelData8 = (unsigned char *)dcminfo.PixelData;
        for (j = 0; j < npix; j++) {
          v8[offset + j] = PixelData8[j];
          if (PixelData8[j] > max8) {
            max8 = PixelData8[j];
          }
          else if (PixelData8[j] < min8) {
            min8 = PixelData8[j];
          }
        }
        // fill missing info
        aDicomInfo[startIndex]->min8 = min8;
        aDicomInfo[startIndex]->max8 = max8;
        free(PixelData8);
        break;

      case 16:
        PixelData16 = (unsigned short *)dcminfo.PixelData;
        for (j = 0; j < npix; j++) {
          v16[offset + j] = PixelData16[j];
          if (PixelData16[j] > max16) {
            max16 = PixelData16[j];
          }
          else if (PixelData16[j] < min16) {
            min16 = PixelData16[j];
          }
        }
        // fill missing info
        aDicomInfo[startIndex]->min16 = min16;
        aDicomInfo[startIndex]->max16 = max16;
        free(PixelData16);
        break;
    }
    exec_progress_callback(n, nfiles, 0, 1);
  }
  // fill missing info
  for (i = 0; i < 3; i++) aDicomInfo[startIndex]->ImagePosition[i] = (firstPosition[i] + lastPosition[i]) / 2.;

  switch (bitsAllocated) {
    case 8:
      return ((void *)v8);
      break;
    case 16:
      return ((void *)v16);
      break;
    default:
      return NULL;
  }
}

/*******************************************************
   SortFiles
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: array of file names, number of elements
   output: array of structures DICOMInfo sorted by
           study date and image number, number of
           studies encountered in the list of files
   comments: conventional MR scans can usually be
             sorted by study date, then image number.
             Other data type (i.e. mosaic functional
             data) may have different discrimination
             fields.

*******************************************************/

void SortFiles(char *fNames[], int nFiles, DICOMInfo ***ptrDicomArray, int *nStudies)
{
  int n, npermut;
  BOOL done;

  DICOMInfo **dicomArray, *storage;
  dicomArray = (DICOMInfo **)calloc(nFiles, sizeof(DICOMInfo *));
  *ptrDicomArray = dicomArray;

  for (n = 0; n < nFiles; n++) {
    dicomArray[n] = (DICOMInfo *)calloc(1, sizeof(DICOMInfo));
    if (dicomArray[n] == NULL) {
      printf(
          "DICOM conversion (SortFiles): "
          "cannot allocate %d bytes\n",
          (int)sizeof(DICOMInfo));
      exit(1);
    }
    GetDICOMInfo(fNames[n], dicomArray[n], FALSE, 1);
    // printf("in sort: n=%d, nframes = %d\n",n,dicomArray[n]->NumberOfFrames);
  }

  // sort by acquisition time, then image number
  done = false;
  while (!done) {
    npermut = 0;
    for (n = 0; n < nFiles - 1; n++)
      if (strcmp(dicomArray[n]->AcquisitionTime, dicomArray[n + 1]->AcquisitionTime) > 0) {
        // 2nd time inferior to first
        storage = dicomArray[n];
        dicomArray[n] = dicomArray[n + 1];
        dicomArray[n + 1] = storage;
        npermut++;
      }
    done = (npermut == 0);
  }

  // sort by image number
  done = false;
  while (!done) {
    npermut = 0;
    for (n = 0; n < nFiles - 1; n++)
      if (strcmp(dicomArray[n]->AcquisitionTime, dicomArray[n + 1]->AcquisitionTime) == 0 &&
          dicomArray[n]->ImageNumber > dicomArray[n + 1]->ImageNumber) {
        storage = dicomArray[n];
        dicomArray[n] = dicomArray[n + 1];
        dicomArray[n + 1] = storage;
        npermut++;
      }
    done = (npermut == 0);
  }

  // check number of studies by looking for differences
  // in acquision times
  *nStudies = 1;
  for (n = 0; n < nFiles - 1; n++)
    if (strcmp(dicomArray[n]->AcquisitionTime, dicomArray[n + 1]->AcquisitionTime) != 0) {
      (*nStudies)++;
    }
  if (*nStudies > 1) {
    printf(
        "WARNING: DICOM conversion, %d different acquisition times "
        "have been identified\n",
        *nStudies);
  }
}

/*******************************************************
   IsDICOM
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: DICOM file name
   output: true if file is DICOM part 10 (version 3.0) compliant
*******************************************************/

int IsDICOM(const char *fname)
{
  int d;
  FILE *fp;
  CONDITION cond;
  DCM_OBJECT *object = 0;
  static int yes = 0;           // statically initialized
  static char file[1024] = "";  // statically initialized

  d = 0;
  if (getenv("FS_DICOM_DEBUG")) {
    d = 1;
  }

  // use the cached value if the fname is the same as privious one
  if (d == 0 && !strcmp(fname, file)) {
    return yes;  // used before
  }
  else {
    yes = 0;              // initialize
    strcpy(file, fname);  // save the current filename
  }

  if (d) {
    printf("Entering IsDICOM (%s)\n", fname);
  }
  fflush(stdout);
  fflush(stderr);

  /* DCM_ACCEPTVRMISMATCH appears to be needed with data produced
     by a recent (2/10/03) Siemens upgrade */

  /* check that the file exists */
  fp = fopen(fname, "r");
  if (fp == NULL) {
    return (0);
  }
  fclose(fp);

  /* check that file is not a directory*/
  if (fio_IsDirectory(fname)) {
    return (0);
  }

  if (d) {
    printf("Opening %s as part10\n", fname);
  }
  COND_PopCondition(1); /* Clears the Condition Stack */
  fflush(stdout);
  fflush(stderr);
  cond = DCM_OpenFile(fname, DCM_PART10FILE | DCM_ACCEPTVRMISMATCH, &object);
  fflush(stdout);
  fflush(stderr);
  if (cond != DCM_NORMAL && d) {
    DCMPrintCond(cond);
    COND_DumpConditions();
  }

  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    if (d) {
      printf("Opening as littleendian\n");
    }
    cond = DCM_OpenFile(fname, DCM_ORDERLITTLEENDIAN | DCM_ACCEPTVRMISMATCH, &object);
    if (cond != DCM_NORMAL && d) {
      DCMPrintCond(cond);
    }
  }

  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    if (d) {
      printf("Opening as bigendian\n");
    }
    cond = DCM_OpenFile(fname, DCM_ORDERBIGENDIAN | DCM_ACCEPTVRMISMATCH, &object);
    if (cond != DCM_NORMAL && d) {
      DCMPrintCond(cond);
    }
  }

  if (cond != DCM_NORMAL) {
    DCM_CloseObject(&object);
    if (d) {
      printf("Opening as format conversion\n");
    }
    cond = DCM_OpenFile(fname, DCM_FORMATCONVERSION | DCM_ACCEPTVRMISMATCH, &object);
    if (cond != DCM_NORMAL && d) {
      DCMPrintCond(cond);
    }
  }
  // if(cond == DCM_NORMAL)
  DCM_CloseObject(&object);

  fflush(stdout);
  fflush(stderr);

  if (d) {
    printf("Leaving IsDICOM (%s)\n", fname);
  }

  COND_PopCondition(1); /* Clears the Condition Stack */

  // cache the current value
  if (cond == DCM_NORMAL) {
    yes = 1;
  }
  else {
    yes = 0;
  }

  return (cond == DCM_NORMAL);
}
/*---------------------------------------------------------------*/
static int DCMPrintCond(CONDITION cond)
{
  switch (cond) {
    case DCM_NORMAL:
      printf("DCM_NORMAL\n");
      break;
    case DCM_ILLEGALOPTION:
      printf("DCM_ILLEGALOPTION\n");
      break;
    case DCM_OBJECTCREATEFAILED:
      printf("DCM_OBJECTCREATEFAILED\n");
      break;
    case DCM_FILEOPENFAILED:
      printf("DCM_FILEOPENFAILED\n");
      break;
    case DCM_FILEACCESSERROR:
      printf("DCM_FILEACCESSERROR\n");
      break;
    case DCM_ELEMENTOUTOFORDER:
      printf("DCM_ELEMENTOUTOFORDER\n");
      break;
    default:
      printf("DCMPrintCond: %d unrecognized\n", (int)cond);
      break;
  }
  return (0);
}

/*******************************************************
   ScanDir
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: directory name
   output: array of files listed in directory, and number of files
*******************************************************/

int ScanDir(const char *PathName, char ***FileNames, int *NumberOfFiles)
{
  char **pfn;
  struct dirent **NameList;
  int i, length, pathlength;

  pathlength = strlen(PathName);
  /* select all directory entries, and sort them by name */
  *NumberOfFiles = scandir(PathName, &NameList, 0, alphasort);

  if (*NumberOfFiles < 0) {
    return -1;
  }

  pfn = (char **)calloc(*NumberOfFiles, sizeof(char *));
  for (i = 0; i < *NumberOfFiles; i++) {
    length = pathlength + strlen(NameList[i]->d_name) + 2;  // must be +2
    // printf("i=%d, l=%d, pl=%d, file=%s\n",i,length,
    // pathlength,NameList[i]->d_name);
    pfn[i] = (char *)calloc(length, sizeof(char));
    sprintf(pfn[i], "%s/%s", PathName, NameList[i]->d_name);
  }

  free(NameList);
  *FileNames = pfn;
  return 0;
}

#ifdef SunOS
#ifndef HAVE_SCANDIR
/* added by kteich for solaris, since it doesn't have them by default. */
/* these funcs Copyright (c) Joerg-R. Hill, December 2000 */
int scandir(const char *dir,
            struct_dirent ***namelist,
            int (*select)(const struct_dirent *),
            int (*compar)(const void *, const void *))
{
  DIR *d;
  struct_dirent *entry;
  register int i = 0;
  size_t entrysize;

  if ((d = opendir(dir)) == NULL) {
    return (-1);
  }

  *namelist = NULL;
  while ((entry = readdir(d)) != NULL) {
    if (select == NULL || (select != NULL && (*select)(entry))) {
      *namelist = (struct_dirent **)realloc((void *)(*namelist), (size_t)((i + 1) * sizeof(struct_dirent *)));
      if (*namelist == NULL) {
        return (-1);
      }
      entrysize = sizeof(struct_dirent) - sizeof(entry->d_name) + strlen(entry->d_name) + 1;
      (*namelist)[i] = (struct_dirent *)malloc(entrysize);
      if ((*namelist)[i] == NULL) {
        return (-1);
      }
      memmove((*namelist)[i], entry, entrysize);
      i++;
    }
  }
  if (closedir(d)) {
    return (-1);
  }
  if (i == 0) {
    return (-1);
  }
  if (compar != NULL) {
    qsort((void *)(*namelist), (size_t)i, sizeof(struct_dirent *), compar);
  }

  return (i);
}
#endif
#ifndef HAVE_ALPHASORT
int alphasort(const void *a, const void *b)
{
  struct_dirent **da = (struct_dirent **)a;
  struct_dirent **db = (struct_dirent **)b;
  return (strcmp((*da)->d_name, (*db)->d_name));
}
#endif
#endif

/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
// Start of "new" dicom reader. New one can read everything but siemens
// mosaics.

/*--------------------------------------------------------------
  DICOMRead2() - generic dicom reader. It should be possible to
  use this for everything but siemens mosaics. There is a 
  siemens specific reader called sdcmLoadVolume().
  --------------------------------------------------------------*/
MRI *DICOMRead2(const char *dcmfile, int LoadVolume)
{
  char **FileNames, *dcmdir;
  DICOMInfo RefDCMInfo, TmpDCMInfo, **dcminfo;
  int nfiles, nframes, nslices, r, c, s, f, err, fid;
  int ndcmfiles, nthfile, mritype = 0, IsDWI, IsPhilipsDWI;
  unsigned short *v16 = NULL;
  unsigned char *v08 = NULL;
  DCM_ELEMENT *element;
  double r0, a0, s0, val;
  MRI *mri;
  FSENV *env;
  std::string tmpfile, tmpfilestdout, FileNameUse, cmd;
  int IsCompressed;

  printf("Starting DICOMRead2()\n");

  if (!fio_FileExistsReadable(dcmfile)) {
    printf("ERROR: file %s does not exist.\n", dcmfile);
    exit(1);
  }
  dcmdir = fio_dirname(dcmfile);
  printf("dcmfile = %s\n", dcmfile);
  printf("dcmdir = %s\n", dcmdir);
  if (!IsDICOM(dcmfile)) {
    setenv("FS_DICOM_DEBUG", "1", 1);
    IsDICOM(dcmfile);
    printf("ERROR: %s is not a dicom file or some other problem\n", dcmfile);
    exit(1);
  }
  // Get info from the reference file
  GetDICOMInfo(dcmfile, &RefDCMInfo, FALSE, 1);
  printf("Ref Series No = %d\n", RefDCMInfo.SeriesNumber);
  if (RefDCMInfo.BitsAllocated != 16 && RefDCMInfo.BitsAllocated != 8) {
    printf("ERROR: bits = %d not supported.\n", RefDCMInfo.BitsAllocated);
    printf("Send email to freesurfer@nmr.mgh.harvard.edu\n");
    return (NULL);
  }

  // Scan directory to get a list of all the files
  err = ScanDir(dcmdir, &FileNames, &nfiles);
  if (err) {
    exit(1);
  }
  printf("Found %d files, checking for dicoms\n", nfiles);

  // Go thru each of those to determine which ones are dicom
  // and belong to the same series.
  ndcmfiles = 0;
  for (nthfile = 0; nthfile < nfiles; nthfile++) {
    // printf("%d %s\n",nthfile,FileNames[nthfile]);
    if (!IsDICOM(FileNames[nthfile])) {
      continue;
    }
    GetDICOMInfo(FileNames[nthfile], &TmpDCMInfo, FALSE, 1);
    if (TmpDCMInfo.SeriesNumber != RefDCMInfo.SeriesNumber) {
      continue;
    }
    ndcmfiles++;
  }
  printf("Found %d dicom files in series.\n", ndcmfiles);

  // Go thru each dicom, make sure it belongs to series, load info
  dcminfo = (DICOMInfo **)calloc(ndcmfiles, sizeof(DICOMInfo *));
  ndcmfiles = 0;
  for (nthfile = 0; nthfile < nfiles; nthfile++) {
    // printf("%d %s\n",nthfile,FileNames[nthfile]);
    if (!IsDICOM(FileNames[nthfile])) {
      continue;
    }
    GetDICOMInfo(FileNames[nthfile], &TmpDCMInfo, FALSE, 1);
    if (TmpDCMInfo.SeriesNumber != RefDCMInfo.SeriesNumber) {
      continue;
    }
    dcminfo[ndcmfiles] = (DICOMInfo *)calloc(1, sizeof(DICOMInfo));
    memmove(dcminfo[ndcmfiles], &TmpDCMInfo, sizeof(DICOMInfo));
    ndcmfiles++;
  }

  // Sort twice, 1st NOT using slice direction, 2nd using slice direction
  // First sort will not use it because Vs=0 from GetDICOMInfo()
  printf("First Sorting\n");
  SortDCMFileInfo(dcminfo, ndcmfiles);
  if (Gdiag_no > 0) {
    for (nthfile = 0; nthfile < ndcmfiles; nthfile++) {
      printf("%d %d %6.2f %6.2f %6.2f %s\n",
             nthfile,
             dcminfo[nthfile]->ImageNumber,
             dcminfo[nthfile]->ImagePosition[0],
             dcminfo[nthfile]->ImagePosition[1],
             dcminfo[nthfile]->ImagePosition[2],
             dcminfo[nthfile]->FileName);
    }
  }
  printf("Computing Slice Direction\n");
  DCMSliceDir(dcminfo, ndcmfiles);
  printf("Second Sorting\n");  // Now uses slice direction
  SortDCMFileInfo(dcminfo, ndcmfiles);
  if (Gdiag_no > 0) {
    for (nthfile = 0; nthfile < ndcmfiles; nthfile++) {
      printf("%d %d %6.2f %6.2f %6.2f %s\n",
             nthfile,
             dcminfo[nthfile]->ImageNumber,
             dcminfo[nthfile]->ImagePosition[0],
             dcminfo[nthfile]->ImagePosition[1],
             dcminfo[nthfile]->ImagePosition[2],
             dcminfo[nthfile]->FileName);
    }
  }
  IsDWI = 0;
  for (nthfile = 0; nthfile < ndcmfiles; nthfile++)
    if (dcminfo[nthfile]->bval > 0) IsDWI = 1;
  IsPhilipsDWI = 0;
  if (IsDWI && strcmp(dcminfo[0]->Manufacturer, "Philips Medical Systems ") == 0) IsPhilipsDWI = 1;
  printf("IsDWI = %d, IsPhilipsDWI = %d\n", IsDWI, IsPhilipsDWI);

  printf("Counting frames\n");
  nframes = DCMCountFrames(dcminfo, ndcmfiles);
  nslices = ndcmfiles / nframes;
  printf("nframes = %d\n", nframes);
  printf("nslices = %d\n", nslices);
  printf("ndcmfiles = %d\n", ndcmfiles);
  if (nslices * nframes != ndcmfiles) {
    if ((nslices * nframes <= ndcmfiles) && (getenv("IGNORE_FRAME_COUNT_CHECK"))) {
      printf(
          "WARNING: the number of frames * number of slices is less than\n"
          "the number of dicom files.\n");
      // dont error exit, because we have enough files to complete the deal
    }
    else {
      printf(
          "ERROR: the number of frames * number of slices does\n"
          "not equal the number of dicom files.\n");
      return (NULL);
    }
  }
  // update reference
  memmove(&RefDCMInfo, dcminfo[0], sizeof(DICOMInfo));

  if (IsPhilipsDWI) {
    // Philips puts the mean DWI as the last frame
    nframes = nframes - 1;
    printf("This is a philips DWI, so ignorning the last frame, nframes = %d\n", nframes);
  }

  int DoRescale = 1;
  int RescaleNeeded = 0;
  // Go through all files and determine whether any need to rescale
  for (nthfile = 0; nthfile < ndcmfiles; nthfile++){
    if(dcminfo[nthfile]->RescaleSlope != 1 || dcminfo[nthfile]->RescaleIntercept != 0){
      RescaleNeeded = 1;
      break;
    }
  }

  if(RescaleNeeded){
    // must explicitly set FS_RESCALE_DICOM=0 to turn off rescaling
    if(getenv("FS_RESCALE_DICOM") != NULL && strcmp(getenv("FS_RESCALE_DICOM"),"0")==0){
      printf("INFO: nontrivial rescale factors are present but will not be applied because\n");
      printf("the FS_RESCALE_DICOM environment variable is set to 0.\n");
      printf("If you want to apply rescaling (intercept and slope), then unset that \n");
      printf("environment variable (or set it to non-zero) and re-run\n");
      printf("\n");
      DoRescale = 0;
    }
  }
  else {
    DoRescale = 0;
    printf("INFO: rescale not needed\n");
  }

  if(DoRescale){
    printf("INFO: applying rescale intercept and slope based on (0028,1052) (0028,1053).\n");
    printf("  If you do not want this, then set FS_RESCALE_DICOM to 0 and rerun.\n");
    mritype = MRI_FLOAT;
  }
  else 
    mritype = MRI_SHORT;

  if(LoadVolume)
    mri = MRIallocSequence(RefDCMInfo.Columns, RefDCMInfo.Rows, nslices, mritype, nframes);
  else {
    mri = MRIallocHeader(RefDCMInfo.Columns, RefDCMInfo.Rows, nslices, mritype, nframes);
    mri->nframes = nframes;
  }
  if (mri == NULL) {
    printf("ERROR: mri alloc failed\n");
    return (NULL);
  }

  // Load the MRI header
  mri->tr = RefDCMInfo.RepetitionTime;
  mri->te = RefDCMInfo.EchoTime;
  mri->ti = RefDCMInfo.InversionTime;
  mri->flip_angle = RefDCMInfo.FlipAngle;
  mri->FieldStrength = RefDCMInfo.FieldStrength;
  mri->width = RefDCMInfo.Columns;
  mri->height = RefDCMInfo.Rows;
  mri->depth = nslices;
  mri->nframes = nframes;
  mri->fov = RefDCMInfo.FieldOfView;
  mri->thick = RefDCMInfo.SliceThickness;
  mri->xsize = RefDCMInfo.xsize;
  mri->ysize = RefDCMInfo.ysize;
  mri->zsize = RefDCMInfo.SliceThickness;

  if (IsDWI) {
    mri->bvals = MatrixAlloc(nframes, 1, MATRIX_REAL);
    mri->bvecs = MatrixAlloc(nframes, 3, MATRIX_REAL);
    mri->bvec_space = BVEC_SPACE_VOXEL;
    if (IsPhilipsDWI) mri->bvec_space = BVEC_SPACE_SCANNER;
  }

  if (getenv("FS_FIX_DICOMS")) {
    if (FZERO(mri->xsize)) {
      mri->xsize = 1;
      fprintf(stderr, "!!!!!!! DICOMread: mri->xsize == 0 - setting to 1 !!!!!!!\n");
    }
    if (FZERO(mri->ysize)) {
      mri->ysize = 1;
      fprintf(stderr, "!!!!!!! DICOMread: mri->ysize == 0 - setting to 1 !!!!!!!\n");
    }
    if (FZERO(mri->zsize)) {
      mri->zsize = 1;
      fprintf(stderr, "!!!!!!! DICOMread: mri->zsize == 0 - setting to 1 !!!!!!!\n");
    }
  }

  mri->x_r = -RefDCMInfo.Vc[0];
  mri->x_a = -RefDCMInfo.Vc[1];
  mri->x_s = +RefDCMInfo.Vc[2];
  mri->y_r = -RefDCMInfo.Vr[0];
  mri->y_a = -RefDCMInfo.Vr[1];
  mri->y_s = +RefDCMInfo.Vr[2];
  mri->z_r = -RefDCMInfo.Vs[0];
  mri->z_a = -RefDCMInfo.Vs[1];
  mri->z_s = +RefDCMInfo.Vs[2];
  if (getenv("FS_FIX_DICOMS")) {
    MATRIX *m;
    MRIreInitCache(mri);
    m = MRIgetRasToVoxelXform(mri);
    if (m == NULL) {
      fprintf(stderr, "!!!!!DICOMread: vox2ras not invertible!!!!!!\n");
      mri->x_r = -1;
      mri->x_a = 0;
      mri->x_s = 0;
      mri->y_r = 0;
      mri->y_a = -1;
      mri->y_s = 0;
      mri->z_r = 0;
      mri->z_a = 0;
      mri->z_s = -1;
      MRIreInitCache(mri);
    }
    m = MRIgetVoxelToRasXform(mri);
    if (m == NULL) {
      fprintf(stderr, "!!!!!DICOMread: ras2vox not invertible!!!!!!\n");
      mri->x_r = -1;
      mri->x_a = 0;
      mri->x_s = 0;
      mri->y_r = 0;
      mri->y_a = -1;
      mri->y_s = 0;
      mri->z_r = 0;
      mri->z_a = 0;
      mri->z_s = -1;
      MRIreInitCache(mri);
    }
  }

  // Phase Enc Direction
  if (RefDCMInfo.PhEncDir == NULL) {
    mri->pedir = strcpyalloc("UNKNOWN");
  }
  else {
    mri->pedir = strcpyalloc(RefDCMInfo.PhEncDir);
    str_toupper(mri->pedir);
  }
  printf("PE Dir = %s (dicom read)\n", mri->pedir);

  // RAS of "first" voxel (ie, at col,row,slice=0)
  r0 = -dcminfo[0]->ImagePosition[0];
  a0 = -dcminfo[0]->ImagePosition[1];
  s0 = +dcminfo[0]->ImagePosition[2];
  // Set c_{ras}
  MRIp0ToCRAS(mri, r0, a0, s0);

  // Return here if only reading in the header
  if (!LoadVolume) {
    for (nthfile = 0; nthfile < ndcmfiles; nthfile++) {
      free(dcminfo[nthfile]);
    }
    free(dcminfo);
    return (mri);
  }

  env = FSENVgetenv();

  printf("Loading pixel data\n");

  if (RefDCMInfo.BitsAllocated == 16 || RefDCMInfo.BitsAllocated == 8) {
    nthfile = 0;
    for (s = 0; s < nslices; s++) {
      for (f = 0; f < nframes; f++) {
        if (Gdiag_no > 0) printf("slice %2d, frame %2d\n", s, f);

        /* Handle compression */
        // If changing this code, make sure to change similar code above
        if (strncmp(dcminfo[nthfile]->TransferSyntaxUID, jpegCompressed_UID, 19) == 0 ||
            strncmp(dcminfo[nthfile]->TransferSyntaxUID, rllEncoded_UID, 19) == 0) {
          // setenv DCMDICTPATH /usr/pubsw/packages/dcmtk/current/share/dcmtk/dicom.dic???
          IsCompressed = 1;
	  tmpfile = std::string(env->tmpdir) + "/" + std::string(env->user) + ".tmp.decompressed.dcm.XXXXXX";
	  std::vector<char> tmpfilechars(tmpfile.begin(), tmpfile.end());
	  tmpfilechars.push_back(0); // Null terminate
	  fid = mkstemp(tmpfilechars.data()); // mkstemp updates the "XXXXXX" at the end
	  tmpfilechars.pop_back(); // Don't need the null terminator
	  tmpfile = std::string(tmpfilechars.begin(), tmpfilechars.end()); // Copy back the modified string
          if (fid == -1) {
            printf("ERROR: could not create temp file for decompression %d\n", fid);
            exit(1);
          }
          close(fid);
          if (strncmp(dcminfo[nthfile]->TransferSyntaxUID, jpegCompressed_UID, 19) == 0) {
            printf("JPEG compressed, decompressing\n");
	    tmpfilestdout = tmpfile + ".dcmdjpeg.out";
	    cmd = std::string("fsdcmdecompress --i ")
	      + std::string(dcminfo[nthfile]->FileName) + " --o " + tmpfile + " --jpeg >& " + tmpfilestdout;
          }
          if (strncmp(dcminfo[nthfile]->TransferSyntaxUID, rllEncoded_UID, 19) == 0) {
            printf("RLE compressed, decompressing\n");
	    tmpfilestdout = tmpfile + ".dcmdlrf.out";
	    cmd = std::string("fsdcmdecompress --i ")
	      + std::string(dcminfo[nthfile]->FileName) + " --o " + tmpfile + " --rle >& " + tmpfilestdout;
          }
          printf("cd %s\n", env->cwd);
          printf("%s\n", cmd.c_str());
          err = system(cmd.c_str());
          if (err != 0) {
            printf("ERROR: %d, see %s for more details\n", err, tmpfilestdout.c_str());
            // Should stream tmpfilestdout to terminal
            exit(1);
          }
          FileNameUse = tmpfile;
        }
        else {
          IsCompressed = 0;
          FileNameUse = std::string(dcminfo[nthfile]->FileName);
        }

        element = GetElementFromFile(FileNameUse.c_str(), 0x7FE0, 0x10);
        if (element == NULL) {
          printf("ERROR: reading pixel data from %s\n", dcminfo[nthfile]->FileName);
          MRIfree(&mri);
          exit(1);
        }
        if (IsCompressed) {
          unlink(tmpfile.c_str());
          unlink(tmpfilestdout.c_str());
        }

	/* Dicoms in MRI and PET are usually unsigned. But in CT, they can be signed.
	   When signed, this code forces negative values to 32k. This can be fixed
	   by changing v08 or v16 to be signed -- but that will mess up the unsigned
	   data. To get around this, the PixelRepresentation (0028,0103) needs to be
	   queried; PR=0=unsigned, PR=1=signed. If this turns into a problem, this
	   code can be changed in this way. DNG 6/7/2018.
	 */
	val = 1;
        v08 = (unsigned char *)(element->d.string);
        v16 = (unsigned short *)(element->d.string);
        for (r = 0; r < RefDCMInfo.Rows; r++) {
          for (c = 0; c < RefDCMInfo.Columns; c++) {
            if(RefDCMInfo.BitsAllocated == 8)  val = *(v08++);
            if(RefDCMInfo.BitsAllocated == 16) val = *(v16++);
	    if(DoRescale)
	      val = val*dcminfo[nthfile]->RescaleSlope + dcminfo[nthfile]->RescaleIntercept;
	    MRIsetVoxVal(mri, c, r, s, f, val);
          }
        }
        if (IsDWI) {
          mri->bvals->rptr[f + 1][1] = dcminfo[nthfile]->bval;
          mri->bvecs->rptr[f + 1][1] = dcminfo[nthfile]->bvecx;
          mri->bvecs->rptr[f + 1][2] = dcminfo[nthfile]->bvecy;
          mri->bvecs->rptr[f + 1][3] = dcminfo[nthfile]->bvecz;
        }
        FreeElementData(element);
        free(element);
        nthfile++;
      }                             // frame
      if (IsPhilipsDWI) nthfile++;  // skip the last file
      exec_progress_callback(f, nframes, s, nslices);
    }  // slice
  }    // 16 bit

  for (nthfile = 0; nthfile < ndcmfiles; nthfile++) free(dcminfo[nthfile]);
  free(dcminfo);

  if (IsDWI) DTIbvecChangeSpace(mri, env->desired_bvec_space);

  FSENVfree(&env);

  return (mri);
}

/*------------------------------------------------------------------
  CompareDCMFileInfo() - function to compare two DCM files for use
  with qsort. The order that files will be sorted is actually the
  reverse of order of the criteria listed here. So the adjacent files
  will first differ by Acq Time, then Echo Number, then Image Number,
  then Slice Position, then Series Number. This order is important
  when determining the number of frames and slices (eg, see
  DCMCountFrames()).
  ----------------------------------------------------------*/
int CompareDCMFileInfo(const void *a, const void *b)
{
  DICOMInfo *dcmfi1;
  DICOMInfo *dcmfi2;
  int n;
  float dv[3], dvsum2, dot, actm1, actm2;

  dcmfi1 = *((DICOMInfo **)a);
  dcmfi2 = *((DICOMInfo **)b);

  /* Sort by Series Number */
  if (dcmfi1->SeriesNumber < dcmfi2->SeriesNumber) {
    return (-1);
  }
  if (dcmfi1->SeriesNumber > dcmfi2->SeriesNumber) {
    return (+1);
  }

  /* ------ Sort by Slice Position -------- */
  /* Compute vector from the first to the second */
  dvsum2 = 0;
  for (n = 0; n < 3; n++) {
    dv[n] = dcmfi2->ImagePosition[n] - dcmfi1->ImagePosition[n];
    dvsum2 += (dv[n] * dv[n]);
  }
  if (dvsum2 > 0) {
    for (n = 0; n < 3; n++) {
      dv[n] /= sqrt(dvsum2); /* normalize */
    }
    /* Compute dot product with Slice Normal vector */
    dot = 0;
    for (n = 0; n < 3; n++) {
      dot += (dv[n] * dcmfi1->Vs[n]);
    }
    // Now Sort by slice position.
    if (dot > +0.5) {
      return (-1);
    }
    if (dot < -0.5) {
      return (+1);
    }
  }

  /* Sort by Image Number (Temporal Sequence) */
  if (dcmfi1->ImageNumber < dcmfi2->ImageNumber) {
    return (-1);
  }
  if (dcmfi1->ImageNumber > dcmfi2->ImageNumber) {
    return (+1);
  }

  /* Sort by Echo Number  */
  if (dcmfi1->EchoNumber < dcmfi2->EchoNumber) {
    return (-1);
  }
  if (dcmfi1->EchoNumber > dcmfi2->EchoNumber) {
    return (+1);
  }

  /* Sort by Acquisition Time */
  sscanf(dcmfi1->AcquisitionTime, "%f", &actm1);
  sscanf(dcmfi2->AcquisitionTime, "%f", &actm2);
  if (actm1 < actm2) {
    return (-1);
  }
  if (actm1 > actm2) {
    return (+1);
  }

  printf(
      "WARNING: files are not found to be different and "
      "cannot be sorted\n");
  printf("File1: %s\n", dcmfi1->FileName);
  printf("File2: %s\n", dcmfi2->FileName);

  return (0);
}

/*-----------------------------------------------------------*/
int SortDCMFileInfo(DICOMInfo **dcmfi_list, int nlist)
{
  qsort(dcmfi_list, nlist, sizeof(DICOMInfo **), CompareDCMFileInfo);
  return (0);
}

/*--------------------------------------------------------------------
  DCMCountFrames() - counts the number of frames in a dicom series.
  The dcmfi_list must have been sored with SortDCMFileInfo(). This
  sorts so that all the files that have the same slice position are
  adjacent in the list, so this just counts the number of files until
  the slice position changes from that of the first file. The number
  of slices is then the number of files/nframes.
  -----------------------------------------------------------*/
int DCMCountFrames(DICOMInfo **dcmfi_list, int nlist)
{
  // Assumes they have been sorted
  int nframes, nth, c, stop;
  double ImgPos0[3], d[3];
  DICOMInfo *dcmfi;

  // Keep track of image position of the first file
  for (c = 0; c < 3; c++) {
    ImgPos0[c] = dcmfi_list[0]->ImagePosition[c];
  }

  nframes = 0;
  for (nth = 0; nth < nlist; nth++) {
    dcmfi = dcmfi_list[nth];
    stop = 0;
    // Compare Image Position of the first file to that
    // of the current file. Stop if they are different.
    for (c = 0; c < 3; c++) {
      d[c] = dcmfi->ImagePosition[c] - ImgPos0[c];
      if (fabs(d[c]) > .0001) {
        stop = 1;
      }
    }
    if (stop) {
      break;
    }
    nframes++;
  }
  return (nframes);
}

/*--------------------------------------------------------------------
  DCMSliceDir() - Computes slice direction such that the first slice
  in the dicom series will sort to the first slice in the volume. For
  EPI, this is important for slice timing correction.
  -----------------------------------------------------------*/
int DCMSliceDir(DICOMInfo **dcmfi_list, int nlist)
{
  // Assumes they have been sorted
  int nth, c, stop;
  double ImgPos0[3], d[3], dlength;
  DICOMInfo *dcmfi;
  MATRIX *Vc, *Vr, *Vs;

  d[0] = d[1] = d[2] = 0.0;

  // Keep track of image position of the first file
  for (c = 0; c < 3; c++) {
    ImgPos0[c] = dcmfi_list[0]->ImagePosition[c];
  }

  stop = 0;
  for (nth = 0; nth < nlist; nth++) {
    // Compare Image Position of the first file to that
    // of the current file. Stop if they are different.
    dcmfi = dcmfi_list[nth];
    stop = 0;
    for (c = 0; c < 3; c++) {
      d[c] = dcmfi->ImagePosition[c] - ImgPos0[c];
      if (fabs(d[c]) > .0001) {
        stop = 1;
      }
    }
    if (stop) {
      break;
    }
  }
  if (stop == 0) {
    printf("\n\n WARNING: it appears that the image position for all slices is the same. \n");
    printf(" This is a problem with the DICOM file. Using direction cosines from col and row\n");
    printf(" to compute slice direction, but the output may be misoriented. \n\n");
    dcmfi = dcmfi_list[0];
    Vc = MatrixAlloc(3, 1, MATRIX_REAL);
    Vr = MatrixAlloc(3, 1, MATRIX_REAL);
    for (c = 0; c < 3; c++) {
      Vc->rptr[c + 1][1] = dcmfi->Vc[c];
      Vr->rptr[c + 1][1] = dcmfi->Vr[c];
    }
    Vs = VectorCrossProduct(Vc, Vr, NULL);
    for (c = 0; c < 3; c++) {
      d[c] = Vs->rptr[c + 1][1];
    }
    MatrixFree(&Vc);
    MatrixFree(&Vr);
    MatrixFree(&Vs);
  }
  // Compute the slice direction cosine based on the first file in the
  // series, which should preserve slice order.
  printf("Vs: %g %g %g\n", d[0], d[1], d[2]);
  dlength = sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]);
  if (FZERO(dlength) && (getenv("FS_FIX_DICOMS"))) {
    d[0] = 1.0;
    fprintf(stderr, "!!!!!DICOMread: warning, direction cosines all zero!!!!!!!!\n");
    dlength = 1;
  }
  for (c = 0; c < 3; c++) {
    d[c] /= dlength;
  }

  printf("Vs: %g %g %g\n", d[0], d[1], d[2]);
  for (nth = 0; nth < nlist; nth++) {
    for (c = 0; c < 3; c++) {
      dcmfi_list[nth]->Vs[c] = d[c];
    }
  }

  return (0);
}

// End of "new" dicom reader routines.
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/

/*******************************************************
   CleanFileNames
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: array of file names, number of elements
   output: array of file names that are DICOM 3.0 compliant
*******************************************************/

int CleanFileNames(char **FileNames, int NumberOfDICOMFiles, char ***CleanedFileNames)
{
  char **pfn;
  int i, j, length;

  pfn = (char **)calloc(NumberOfDICOMFiles, sizeof(char *));
  for (i = 0, j = 0; j < NumberOfDICOMFiles; i++) {
    if (IsDICOM(FileNames[i])) {
      length = strlen(FileNames[i]) + 1;
      pfn[j] = (char *)calloc(length, sizeof(char));
      strcpy(pfn[j], FileNames[i]);
      j++;
    }
  }

  printf("Found %d DICOM Files\n", j);

  *CleanedFileNames = pfn;
  return DCM_NOERROR;
}

/* kteich - renamed this to myRound because round() was causing
   versioning problems */
int myRound(double d)
{
  double c, f, ddown, dup;

  c = ceil(d);
  f = floor(d);
  ddown = d - f;
  dup = c - d;
  if (ddown < dup) {
    return (int)f;
  }
  else {
    return (int)c;
  }
}

/*******************************************************
   RASFromOrientation
   Author: E'beth Haley
   Date: 08/01/2002
   input: structure MRI, structure DICOMInfo
   output: fill in MRI structure RAS-related fields
   called-by: DICOMInfo2MRI
*******************************************************/

void DebugPrint(FILE *fp, const char *msg, ...)
{
#ifndef __OPTIMIZE
  va_list ap;
  va_start(ap, msg);
  vfprintf(fp, msg, ap);
  va_end(ap);
#endif
}

int RASFromOrientation(MRI *mri, DICOMInfo *dcm)

/*E*

This code is intended to repair and generalize
RASFromOrientation, which is used for GE DICOM scans, to make
it more like the Doug Greve Siemens DICOM code, in that it
should refer to the direction-cosines matrix and not just try
to match up r/a/s with x/y/z.

meta: Doug Greve notes that we'd eventually like to be able to
load n frames at a time, and in that case, we'd probably have
to rewrite Sebastien's code from scratch - the least reason for
that is S uses "frames" to refer to slices=files.  In the
meantime, I'm just going to reimplement the function
RASFromOrientation, Dougstyle.

Note to self: the DICOMInfo struct is a Sebastienism.  Only
Sebastien code uses it.  The Siemens stuff struct is
SDCMFILEINFO.

Note to self: okay, the DICOMInfo knows the fname it was
created from, so we can use that to get things there are tags
for that don't appear in the DICOMInfo itself.

Siemens version: sdcmLoadVolume returns an MRI which it gets by
passing a dcmfilename to ScanSiemensDCMDir which returns a list
of SDCMFILEINFO's (that it gets from GetSDCMFileInfo which gets
stuff from a dcmfile)

RASFromOrientation is passed a DICOMInfo (that has been partway
filled in by DICOMInfo2MRI() - any other info we need we can
get from dcm->fname), and sets these elements of the mri
struct:

mri->x_r, mri->y_r, mri->z_r,
mri->x_a, mri->y_a, mri->z_a,
mri->x_s, mri->y_s, mri->z_s,
mri->c_r, mri->c_a, mri->c_s,
mri->ras_good_flag=1

By the time that calls RASFromOrientation(), the mri is already
partly filled in by DICOMInfo2MRI().  So it already has:

mri->width = dcm->Columns;
mri->height = dcm->Rows;
mri->depth = dcm->NumberOfFrames;
nvox = dcm->Rows*dcm->Columns*dcm->NumberOfFrames;
mri->fov = dcm->FieldOfView;
mri->thick = dcm->SliceThickness;
mri->xsize = dcm->xsize;
mri->ysize = dcm->ysize;
mri->zsize = dcm->SliceThickness;
mri->nframes = 1;
mri->tr = dcm->RepetitionTime;
mri->te = dcm->EchoTime;
mri->ti = dcm->InversionTime;
mri->flip_angle = dcm->FlipAngle;
mri->imnr0 = 1;
mri->imnr1 = dcm->NumberOfFrames;

We'll need e.g. dcm->NumberOfFrames (NSlices).

We'll more or less use:

T =  [ (Mdc * Delta)    P_o ]
[    0    0    0    1  ]

where Mdc = [ c r s ] is the matrix of direction cosines, which
we'll get from dcmImageDirCos(dcm->FileName, cx,cy,cz,
rx,ry,rz) plus a few lines (that should probably become a
dcmSliceDirCos()) to get the slice column, sx,sy,sz.

Well, okay, we're not using the 4x4 T, but rather
T3x3=Mdc*Delta and P_o separately.

Delta = pixel size in mm = get from dcmGetVolRes() (or use {
dcm->xsize, dcm->ysize, dcm->SliceThickness }.  BTW, this
code's Delta is the vector that the matrix Delta should be the
diag of: diag(Delta_c, Delta_r, Delta_s).

get P_o = [ P_ox P_oy P_oz ] from dcmImagePosition(char
*dcmfile, float *x, float *y, float *z) on first file - or
dcm->FirstImagePosition - wait, those are off by factors of -1,
-1, 1 from each other.

The reason is that DICOM coordinate system is LPS, meanwhile
we use RAS coordinate system. that is why (-1,-1,1) multiplication
occurs.

Note that FirstImagePosition and LastImagePosition are
in LPS coordinate system.  In order to obtain the RAS direction cosine,
you have to change them to RAS by mutliplying (-1,-1,1).

n[3] = [ width, height, depth ] = { dcm->Columns -1, dcm->Rows
-1, dcm->NumberOfFrames -1 }

and finally c_ras = center = T3x3 * n/2

*/
{
  float c[3], r[3], s[3], Delta[3], c_ras[3], Mdc[3][3], T3x3[3][3], rms;  // absIO[6];

  int i, j, n[3];

  // note Doug does the conversion from LPS to RAS and thus
  // these have the correct values
  dcmImageDirCos(dcm->FileName, c, c + 1, c + 2, r, r + 1, r + 2);

  // slice direction is not given in dicom
  //
  // we had to calculate
  // one way is to use the image position to do so
  // if number of slice is greater than one.
  //
  // The image position is given in DICOM coordinate system, which is a
  // LPS (left-posterior-superior) coordinate system.  We are interested in
  // RAS (right-anterior-superior) coordinate system.
  // The first two coordinates sign flips, but not the last.
  //
  // Actually if you used dcmImagePosition(), then this would be the case!
  //
  // Mixing Sebastian routines with Doug routines is a bad thing to do :-(...
  //
  if (dcm->NumberOfFrames > 1)  // means = NSlices > 0
  {
    rms = 0.0;
    for (i = 0; i < 3; i++) {
      if (i < 2)  // convert LP to RA
      {
        s[i] = dcm->FirstImagePosition[i] - dcm->LastImagePosition[i];
      }
      else  // S stays the same
      {
        s[i] = dcm->LastImagePosition[i] - dcm->FirstImagePosition[i];
      }
      rms += s[i] * s[i];
    }
    rms = sqrt(rms);
    for (i = 0; i < 3; i++) {
      s[i] /= rms;
    }
  }
  else {
    // build s[i] from c[i] and r[i]
    // there is an ambiguity in sign due to left-handed
    // or right-handed coords
    // we pick right-handed.  It does not matter anyway.
    // note that c[] and r[] are normalized already and
    // thus no need to calculate
    // the normalization factor.
    s[0] = c[1] * r[2] - c[2] * r[1];
    s[1] = c[2] * r[0] - c[0] * r[2];
    s[2] = c[0] * r[1] - c[1] * r[0];
  }

  for (i = 0; i < 3; i++) {
    Mdc[i][0] = c[i];
    Mdc[i][1] = r[i];
    Mdc[i][2] = s[i];
  }

  // when one slice is given, it is zero
  dcmGetVolRes(dcm->FileName, Delta, Delta + 1, Delta + 2);
  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      T3x3[i][j] = Mdc[i][j] * Delta[j];
    }
  }

  /* fill in n, number */
  n[0] = dcm->Columns - 1;
  n[1] = dcm->Rows - 1;
  n[2] = dcm->NumberOfFrames - 1; /*E* This is Slices :( */

  /*E* The -1 is for how many center-to-center lengths there are.  */
  /*E*
    Doug: c_ras = center = T * n /2.0

    i.e. P_o + Mdc * Delta * n/2

    Okay, that's: take first voxel and go to center by adding on half
    of [ reorder and rescale the nc/nr/ns ] = lengths of the sides of
    the volume.

    dcmImagePosition differs from DICOM file first image position by
    factors of -1 -1 1.  Sebastien code also multiplies by
    diag{-1,-1,1}, so it should be fine to use dcmImagePosition.
  */
  dcmImagePosition(dcm->FileName, c_ras, c_ras + 1, c_ras + 2); /*E* = P_o */

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++) {
      c_ras[i] += T3x3[i][j] * (float)n[j] / 2.;
    }
  }
  // printing a matrix, you can call call MatrixPrint(stdout, m) in gdb
  // no need for #ifdef ..
  /*E* Now set mri values - the direction cosines.  I stole them
    straight from Sebastien's code, but - Are these backwards??  Ok,
    05aug02, transposed them: */

  mri->x_r = Mdc[0][0];
  mri->y_r = Mdc[0][1];
  mri->z_r = Mdc[0][2];
  mri->x_a = Mdc[1][0];
  mri->y_a = Mdc[1][1];
  mri->z_a = Mdc[1][2];
  mri->x_s = Mdc[2][0];
  mri->y_s = Mdc[2][1];
  mri->z_s = Mdc[2][2];
  mri->c_r = c_ras[0];
  mri->c_a = c_ras[1];
  mri->c_s = c_ras[2];
  mri->ras_good_flag = 1;

  return 0;
}

/*******************************************************
   DICOM16To8
   Author: Sebastien Gicquel
   Date: 06/06/2001
   input: array of 16-bit pixels, number of elements
   output: array of 8-bit pixels
*******************************************************/
unsigned char *DICOM16To8(unsigned short *v16, int nvox)
{
  unsigned char *v8;
  int i;
  double min16, max16, min8, max8, ratio;

  min8 = 0;
  max8 = 255;

  v8 = (unsigned char *)calloc(nvox, sizeof(unsigned char));
  if (v8 == NULL) {
    printf("DICOMInfo2MRI: can't allocate %d bytes\n", nvox);
    exit(1);
  }

  for (i = 0, min16 = 65535, max16 = 0; i < nvox; i++) {
    if (v16[i] > max16) {
      max16 = (double)v16[i];
    }
    if (v16[i] < min16) {
      min16 = (double)v16[i];
    }
  }

  ratio = (max8 - min8) / (max16 - min16);
  for (i = 0; i < nvox; i++) {
    v8[i] = (unsigned char)((double)(v16[i]) * ratio);
  }

  return v8;
}
/*******************************************************
   DICOMInfo2MRI
   Author: Sebastien Gicquel
   Date: 06/05/2001
   input: structure DICOMInfo, pixel data
   output: fill in MRI structure, including the image
*******************************************************/

int DICOMInfo2MRI(DICOMInfo *dcm, void *data, MRI *mri)
{
  long n;
  int i, j, k;
  unsigned char *data8;
  unsigned short *data16;

  // fill in the fields
  strcpy(mri->fname, dcm->FileName);

  mri->width = dcm->Columns;
  mri->height = dcm->Rows;
  mri->depth = dcm->NumberOfFrames;
  // nvox = dcm->Rows * dcm->Columns * dcm->NumberOfFrames;

  mri->fov = dcm->FieldOfView;
  mri->thick = dcm->SliceThickness;
  mri->xsize = dcm->xsize;
  mri->ysize = dcm->ysize;
  mri->zsize = dcm->SliceThickness;
  mri->nframes = 1;

  mri->tr = dcm->RepetitionTime;
  mri->te = dcm->EchoTime;
  mri->ti = dcm->InversionTime;
  mri->flip_angle = dcm->FlipAngle;

  mri->imnr0 = 1;
  mri->imnr1 = dcm->NumberOfFrames;

  RASFromOrientation(mri, dcm);

  switch (dcm->BitsAllocated) {
    case 8:
      data8 = (unsigned char *)data;
      for (k = 0, n = 0; k < dcm->NumberOfFrames; k++)
        for (j = 0; j < dcm->Rows; j++)
          for (i = 0; i < dcm->Columns; i++, n++) {
            MRIvox(mri, i, j, k) = data8[n];
          }
      break;
    case 16:
      data16 = (unsigned short *)data;
      for (k = 0, n = 0; k < dcm->NumberOfFrames; k++)
        for (j = 0; j < dcm->Rows; j++)
          for (i = 0; i < dcm->Columns; i++, n++) {
            MRISvox(mri, i, j, k) = data16[n];
          }
      break;
  }

  return 0;
}

/*******************************************************
   DICOMRead
   Author: Sebastien Gicquel
   Date: 06/04/2001
   input: directory name, boolean
   output: MRI structure, including the image if input
           boolean is true
   FileName points to a DICOM filename. All DICOM files
   in same directory will be read

   This routine is used for non-Siemens DICOM files
*******************************************************/

int DICOMRead(const char *FileName, MRI **mri, int ReadImage)
{
  MRI *pmri = NULL;
  const char *c;
  char **CleanedFileNames, **FileNames, PathName[256];
  int i, NumberOfFiles, NumberOfDICOMFiles, nStudies, error;
  int length;
  DICOMInfo **aDicomInfo;
  unsigned char *v8 = NULL;
  unsigned short *v16 = NULL;
  FILE *fp;
  int numFiles;
  int inputIndex;
  int nextIndex;
  int *startIndices;
  int count;

  for (i = 0; i < NUMBEROFTAGS; i++) {
    IsTagPresent[i] = false;
  }

  // find pathname from filename
  c = strrchr(FileName, '/');
  if (c == NULL) {
    PathName[0] = '\0';
  }
  else {
    length = (int)(c - FileName);
    strncpy(PathName, FileName, length);
    PathName[length] = '/';
    PathName[length + 1] = '\0';
  }

  // first check whether this is a siemens dicom file or not
  if (IsSiemensDICOM(FileName)) {
    *mri = sdcmLoadVolume(FileName, 1, -1);
  }

  // scan directory
  error = ScanDir(PathName, &FileNames, &NumberOfFiles);
  if (error) {
    fprintf(stderr, "ScanDir encountered an error\n");
  }

  for (i = 0, NumberOfDICOMFiles = 0; i < NumberOfFiles; i++) {
    if (IsDICOM(FileNames[i])) {
      NumberOfDICOMFiles++;
    }
  }
  // no DICOM files in directory
  if (NumberOfDICOMFiles == 0) {
    fprintf(stderr, "no DICOM files found. Exit\n");
    exit(1);
  }
  else {
    printf("%d DICOM 3.0 files in list\n", NumberOfDICOMFiles);
  }

  // remove non DICOM file names
  CleanFileNames(FileNames, NumberOfDICOMFiles, &CleanedFileNames);

  // sort DICOM files by study date, then image number
  // apparently, number of frames is set here
  printf("sorting\n");
  SortFiles(CleanedFileNames, NumberOfDICOMFiles, &aDicomInfo, &nStudies);

  printf("postsort: nframes = %d\n", aDicomInfo[0]->NumberOfFrames);

  // if more than 1 studies in the file then do more work to create
  // a list of files which to make an output file.   Then use that
  // to print out the information.
  inputIndex = 0;
  count = 0;
  if (nStudies > 1) {
    // create an array of starting indices
    startIndices = (int *)calloc(nStudies, sizeof(int));

    printf("Generating log file dicom.log\n");
    fp = fopen("dicom.log", "w");
    if (fp == NULL) {
      printf("Can not create file dicom.log\n");
    }
    else {
      printf("Starting filenames for the studies\n");
      for (i = 0; i < NumberOfDICOMFiles; i++) {
        fprintf(fp, "%s\t%d\n", aDicomInfo[i]->FileName, aDicomInfo[i]->ImageNumber);
        // find the index for the input file
        if (!strcmp(FileName, aDicomInfo[i]->FileName)) {
          inputIndex = i;
        }
        // print starting filename for the studies
        if (i == 0) {
          startIndices[count] = i;
          count++;
          printf(" %s\n", aDicomInfo[i]->FileName);
        }
        // get a different studies
        else if (i != 0 && strcmp(aDicomInfo[i]->AcquisitionTime, aDicomInfo[i - 1]->AcquisitionTime) != 0) {
          startIndices[count] = i;
          count++;
          printf(" %s\n", aDicomInfo[i]->FileName);
        }
      }
      fclose(fp);
    }
    // ok found the Dicominfo for the input file and its image number.
    // look for the starting index for this selected list
    nextIndex = inputIndex;
    for (i = 0; i < count; i++) {
      // startIndices is bigger than inputIndex, then it is not the
      // starting index of the selected file.
      if (startIndices[i] > inputIndex) {
        if (i > 0) {
          inputIndex = startIndices[i - 1];
        }
        // inputIndex is the previous one
        nextIndex = startIndices[i];  // next series starting here
        break;
      }
    }
    // could end up nextIndex = inputIndex
    if (inputIndex == nextIndex) {
      nextIndex = inputIndex + 1;
    }

    free((void *)startIndices);
  }
  else  // only 1 studies
  {
    nextIndex = inputIndex + NumberOfDICOMFiles;
  }

  // make sure that it is within the range
  if (inputIndex < 0 || inputIndex > NumberOfDICOMFiles - 1) {
    fprintf(stderr,
            "Failed to find the right starting "
            "index for the image %s.  Exit.\n",
            FileName);
    exit(1);
  }

  // no longer relies on the image number.
  // verify ImageNumber to be 1
  // if (aDicomInfo[inputIndex]->ImageNumber != 1)
  // {
  //   fprintf(stderr, "The image number of the first
  //        image was not 1 for the list.  Exit.\n");
  //   exit(1);
  // }
  if (nStudies > 1) {
    printf("Use %s as the starting filename for this extraction.\n", aDicomInfo[inputIndex]->FileName);
    printf(
        "If like to use the different one, please pick "
        "one of the files listed above.\n");
  }
  numFiles = nextIndex - inputIndex;

  // verify with aDicomInfo->NumberOfFrames;
  if (aDicomInfo[inputIndex]->NumberOfFrames != numFiles) {
    printf(
        "WARNING: NumberOfFrames %d != Found Count of slices %d.\n", aDicomInfo[inputIndex]->NumberOfFrames, numFiles);
  }
  // find the total number of files belonging to this list
  // the list has been sorted and thus look for the end.
  // numFiles = aDicomInfo[inputIndex]->NumberOfFrames; this is not reliable
  // use the ImageNumber directly
  //
  // Watch out!!!!!!
  // Note that some dicom image consists of many images put into one
  // the increment of image number can be 30
  // If the image number increament is greater than one, we bail out
  //
  // non siemens dicom file cannot handle mosaic image
  // nextIndex - inputIndex is the number of files.
  // last image number - first image number > numfiles then must be mosaic
  if ((aDicomInfo[nextIndex - 1]->ImageNumber - aDicomInfo[inputIndex]->ImageNumber) > numFiles) {
    fprintf(stderr,
            "Currently non-Siemens mosaic image "
            "cannot be handled.  Exit\n");
    exit(1);
  }

  printf("about to read in image\n");
  switch (aDicomInfo[inputIndex]->BitsAllocated) {
    case 8:
      printf("ReadDICOMImage2 8 bit\n");
      v8 = (unsigned char *)ReadDICOMImage2(numFiles, aDicomInfo, inputIndex);
      pmri = MRIallocSequence(aDicomInfo[inputIndex]->Columns, aDicomInfo[inputIndex]->Rows, numFiles, MRI_UCHAR, 1);
      DICOMInfo2MRI(aDicomInfo[inputIndex], (void *)v8, pmri);
      free(v8);
      break;
    case 16:
      printf("ReadDICOMImage2 16 bit\n");
      v16 = (unsigned short *)ReadDICOMImage2(numFiles, aDicomInfo, inputIndex);
      pmri = MRIallocSequence(aDicomInfo[inputIndex]->Columns, aDicomInfo[inputIndex]->Rows, numFiles, MRI_SHORT, 1);
      DICOMInfo2MRI(aDicomInfo[inputIndex], (void *)v16, pmri);
      free(v16);
      break;
  }
  printf("postread: nframes = %d\n", aDicomInfo[inputIndex]->NumberOfFrames);

  // display only first DICOM header
  PrintDICOMInfo(aDicomInfo[inputIndex]);

  *mri = pmri;

  return 0;
}

/*
  \fn int dcmGetDWIParams()
  \brief Extracts DWI parameters (bvals and bvecs) from a dicom file
  (Siemens or GE).  The gradient direction should be in voxel
  coordinates. Returns non-zero if there is a problem. If it cannot find
  the bvalues and bvects, then zero is still returned but bvals and bvecs
  are set to 0.
 */
int dcmGetDWIParams(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec)
{
  DCM_ELEMENT *e;
  CONDITION cond;
  DCM_TAG tag;
  unsigned int rtnLength;
  void *Ctx = NULL;
  int err;
  double rms;

  *pbval = 0;
  *pxbvec = 0;
  *pybvec = 0;
  *pzbvec = 0;

  Ctx = NULL;
  tag = DCM_MAKETAG(0x8, 0x70);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    free(e);
    return (10);
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  if (cond != DCM_NORMAL) {
    free(e);
    return (10);
  }

  if (strcmp(e->d.string, "SIEMENS") == 0 || strcmp(e->d.string, "SIEMENS ") == 0) {
    if (Gdiag_no > 0) printf("Attempting to get DWI Parameters from Siemens DICOM\n");
    err = dcmGetDWIParamsSiemens(dcm, pbval, pxbvec, pybvec, pzbvec);
    if (err) return (err);
  }
  else if (strcmp(e->d.string, "GE MEDICAL SYSTEMS") == 0) {
    if (Gdiag_no > 0) printf("Attempting to get DWI Parameters from GE DICOM\n");
    err = dcmGetDWIParamsGE(dcm, pbval, pxbvec, pybvec, pzbvec);
    if (err) return (err);
  }
  else if (strcmp(e->d.string, "Philips Medical Systems ") == 0) {
    // Note: need space at the end of 'Systems' above
    if (Gdiag_no > 0) printf("Attempting to get DWI Parameters from Phlips DICOM\n");
    err = dcmGetDWIParamsPhilips(dcm, pbval, pxbvec, pybvec, pzbvec);
    if (err) return (err);
  }
  else {
    printf("ERROR: don't know how to get DWI parameters from --%s--\n", e->d.string);
    fflush(stdout);
    return (1);
  }
  FreeElementData(e);

  rms = sqrt((*pxbvec) * (*pxbvec) + (*pybvec) * (*pybvec) + (*pzbvec) * (*pzbvec));
  if (Gdiag_no > 0) printf("%lf %lf %lf %lf %lf\n", *pbval, *pxbvec, *pybvec, *pzbvec, rms);

  if (*pbval != 0 && (fabs(*pxbvec) > 1.0 || fabs(*pybvec) > 1.0 || fabs(*pzbvec) > 1.0 || fabs(rms - 1) > .001) &&
      fabs(rms) > .001) {
    printf("%lf %lf %lf %lf %lf \n", *pbval, *pxbvec, *pybvec, *pzbvec, rms);
    printf("WARNING: These don't look like reasonable DWI params.\n");
  }

  if (*pbval == 0) {
    *pxbvec = 0;
    *pybvec = 0;
    *pzbvec = 0;
  }

  return (0);
}

int dcmGetDWIParamsPhilips(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec)
{
  DCM_ELEMENT *e;
  CONDITION cond;
  DCM_TAG tag;
  unsigned int rtnLength;
  void *Ctx = NULL;

  if (Gdiag_no > 0) printf("Entering dcmGetDWIParamsPhilips()\n");

  // since it's difficult to check whether Philip's DWI info exists, this function should
  // return 0 after a read failure instead of returning an error (assuming that the DWI
  // data doesn't exist). However, if the user has specified 'FS_LOAD_DWI', we should assume
  // that DWI data definitely exists, and a read failure should throw an error
  char *pc;
  int forceDWI = 1;
  pc = getenv("FS_LOAD_DWI");
  if (pc == NULL) {
    forceDWI = 0;
  } else if (strcmp(pc, "0") == 0) {
    forceDWI = 0;
  }

  *pbval = 0;
  *pxbvec = 0;
  *pybvec = 0;
  *pzbvec = 0;

  // bvalue
  Ctx = NULL;
  tag = DCM_MAKETAG(0x18, 0x9087);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    free(e);
    if (forceDWI) {
      return(7);
    } else {
      return(0);
    }
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  if (cond != DCM_NORMAL) {
    free(e);
    return (8);
  }
  *pbval = e->d.fd[0];
  free(e);

  // gradient vector
  Ctx = NULL;
  tag = DCM_MAKETAG(0x18, 0x9089);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    free(e);
    return (9);
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  if (cond != DCM_NORMAL) {
    free(e);
    return (10);
  }
  // Not sure about the sign here. LZ gave DNG some bvecs that appear to be
  // (+1,-1,-1), but changing
  (*pxbvec) = -1 * e->d.fd[0];  // ???
  (*pybvec) = -1 * e->d.fd[1];  // ???
  (*pzbvec) = +1 * e->d.fd[2];  // ???
  free(e);

  if (Gdiag_no > 0) printf("%lf %lf %lf %lf\n", *pbval, *pxbvec, *pybvec, *pzbvec);

  if (*pbval == 0) {
    *pxbvec = 0;
    *pybvec = 0;
    *pzbvec = 0;
  }

  return (0);
}

/*
  \fn int dcmGetDWIParamsGE()
  \brief Extracts DWI info from a GE file. No transformation is
   applied to the gradients because it is assumed that they are in
   voxel space already. However, gradients are converted to RAS.
   Returns 0 if no error.
*/
int dcmGetDWIParamsGE(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec)
{
  DCM_ELEMENT *e;
  CONDITION cond;
  DCM_TAG tag;
  unsigned int rtnLength, n;
  void *Ctx = NULL;

  if (Gdiag_no > 0) printf("Entering dcmGetDWIParamsGE()\n");

  *pbval = 0;
  *pxbvec = 0;
  *pybvec = 0;
  *pzbvec = 0;

  // bvalue
  Ctx = NULL;
  tag = DCM_MAKETAG(0x43, 0x1039);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    free(e);
    return (7);
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  if (cond != DCM_NORMAL) {
    free(e);
    return (7);
  }
  for (n = 0; n < strlen(e->d.string); n++)
    if (e->d.string[n] == '\\') e->d.string[n] = ' ';
  sscanf(e->d.string, "%lf", pbval);
  free(e);

  // x part of gradient vector
  Ctx = NULL;
  tag = DCM_MAKETAG(0x19, 0x10bb);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    free(e);
    return (7);
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  if (cond != DCM_NORMAL) {
    free(e);
    return (7);
  }
  sscanf(e->d.string, "%lf", pxbvec);
  (*pxbvec) *= -1;  // convert from LPS to RAS
  free(e);

  // y part of gradient vector
  Ctx = NULL;
  tag = DCM_MAKETAG(0x19, 0x10bc);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    free(e);
    return (7);
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  if (cond != DCM_NORMAL) {
    free(e);
    return (7);
  }
  sscanf(e->d.string, "%lf", pybvec);
  (*pybvec) *= -1;  // convert from LPS to RAS
  free(e);

  // z part of gradient vector
  Ctx = NULL;
  tag = DCM_MAKETAG(0x19, 0x10bd);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    free(e);
    return (7);
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  if (cond != DCM_NORMAL) {
    free(e);
    return (7);
  }
  sscanf(e->d.string, "%lf", pzbvec);
  free(e);

  if (Gdiag_no > 0) printf("%lf %lf %lf %lf\n", *pbval, *pxbvec, *pybvec, *pzbvec);

  if (*pbval == 0) {
    *pxbvec = 0;
    *pybvec = 0;
    *pzbvec = 0;
  }

  return (0);
}

/*
  \fn int dcmGetDWIParamsSiemens()
  \brief Extracts DWI info from a Siemens file, transforms gradients
  to voxel coordinates. First, looks for bval tag 0x19 0x100c. If this
  exists, then gets bvec from 0x100e. If not, then uses an alternative
  method. Returns 0 if everything ok. Gradients transformed to RAS.
 */
int dcmGetDWIParamsSiemens(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec)
{
  DCM_ELEMENT *e;
  CONDITION cond;
  DCM_TAG tag;
  unsigned int rtnLength;
  void *Ctx = NULL;
  int err;
  double Vcx, Vcy, Vcz, Vrx, Vry, Vrz, Vsx, Vsy, Vsz;
  MATRIX *Mdc, *G, *G2, *iMdc;

  if (Gdiag_no > 0) printf("Entering dcmGetDWIParamsSiemens()\n");

  *pbval = 0;
  *pxbvec = 0;
  *pybvec = 0;
  *pzbvec = 0;

  // Get the bvalue
  Ctx = NULL;
  tag = DCM_MAKETAG(0x19, 0x100c);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    // The bvalue tag does not exist, try alternative method
    if (Gdiag_no > 0) printf("  bval 0x19,0x100c does not exist, try alternative method\n");
    err = dcmGetDWIParamsSiemensAlt(dcm, pbval, pxbvec, pybvec, pzbvec);
    if (err) return (err);
  }
  else {
    // The bvalue tag does exist, get gradients
    if (Gdiag_no > 0) printf("  bval 0x19,0x100c does exist, getting gradients\n");
    AllocElementData(e);
    cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
    if (cond != DCM_NORMAL) {
      free(e);
      return (6);
    }
    sscanf(e->d.string, "%lf", pbval);
    free(e);

    if (*pbval > 0) {
      Ctx = NULL;
      tag = DCM_MAKETAG(0x19, 0x100e);
      e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
      cond = DCM_GetElement(&dcm, tag, e);
      if (cond != DCM_NORMAL) {
        if (Gdiag_no > 0) printf("  bvec 0x19,0x100e does not exist, returning 7\n");
        free(e);
        return (7);
      }
      AllocElementData(e);
      cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
      if (cond != DCM_NORMAL) {
        if (Gdiag_no > 0) printf("  could not get bvec 0x19,0x100e returning 8\n");
        free(e);
        return (8);
      }
      *pxbvec = e->d.fd[0];
      *pybvec = e->d.fd[1];
      *pzbvec = e->d.fd[2];
      (*pxbvec) *= -1;  // convert from LPS to RAS
      (*pybvec) *= -1;  // convert from LPS to RAS
      free(e);
    }
  }

  err = dcmImageDirCosObject(dcm, &Vcx, &Vcy, &Vcz, &Vrx, &Vry, &Vrz);
  if (err) return (9);

  if (*pbval == 0) {
    *pxbvec = 0;
    *pybvec = 0;
    *pzbvec = 0;
  }

  /* Mdc maps from vox coords to scanner coords.
     ImageDirCos2Slice() just computes the slice DC using a cross
     product of the in-plane DC. This makes the sign of slice DC is
     arbitrary, but it should not matter (?)*/
  Mdc = ImageDirCos2Slice(Vcx, Vcy, Vcz, Vrx, Vry, Vrz, &Vsx, &Vsy, &Vsz);
  iMdc = MatrixInverse(Mdc, NULL);  // iMdc maps from scanner coords to vox coords
  G = MatrixAlloc(3, 1, MATRIX_REAL);
  G->rptr[1][1] = *pxbvec;
  G->rptr[2][1] = *pybvec;
  G->rptr[3][1] = *pzbvec;
  G2 = MatrixMultiplyD(iMdc, G, NULL);
  if (Gdiag_no > 0) {
    printf("Transforming gradient directions into voxel coordinates\n");
    printf("DC: %f %f %f \n%f %f %f\n%f %f %f\n", Vcx, Vcy, Vcz, Vrx, Vry, Vrz, Vsx, Vsy, Vsz);
    printf("iMdc = \n");
    MatrixPrint(stdout, iMdc);
    printf("Before: %lf %lf %lf %lf\n", *pbval, *pxbvec, *pybvec, *pzbvec);
  }
  *pxbvec = G2->rptr[1][1];
  *pybvec = G2->rptr[2][1];
  *pzbvec = G2->rptr[3][1];
  MatrixFree(&Mdc);
  MatrixFree(&iMdc);
  MatrixFree(&G);
  MatrixFree(&G2);

  if (Gdiag_no > 0) printf("After:  %lf %lf %lf %lf\n", *pbval, *pxbvec, *pybvec, *pzbvec);
  return (0);
}

/*
  \fn int dcmGetDWIParamsSiemensAlt()
  \brief This is an alternative method to extract DWI info from a
  Siemens file. It looks in 0x29 0x1010. This is a nasty string with
  all kinds of control characters, etc.  It looks for key words in
  this string and blindly extracts values based on proximity to the
  key word. The gradients are NOT transformed to image
  space. Gradients are transformed to RAS. This alternate method is
  probably less reliable than the main method and is probably only
  needed for older scanners.  To report non-zero values for the 4
  params, the bvalue must be present and > 0, there must be 3 bvecs
  found, the abs(bvecs) must be < 1, and the rms must be close to 1.
  This routine is unreliable and quite finicky because it is just
  looking for key words. The string may contain the key words even if
  the sequence was not diffusion weighted. This makes it hard to
  do error checking.
 */
int dcmGetDWIParamsSiemensAlt(DCM_OBJECT *dcm, double *pbval, double *pxbvec, double *pybvec, double *pzbvec)
{
  DCM_ELEMENT *e;
  CONDITION cond;
  DCM_TAG tag;
  unsigned int rtnLength, m;
  void *Ctx = NULL;
  int k, bval_flag, bvec_flag, Allow;
  char c, tmpstr[2000], *pc;
  double val, rms;

  if (Gdiag_no > 0) printf("Entering dcmGetDWIParamsSiemensAlt()\n");

  *pbval = 0;
  *pxbvec = 0;
  *pybvec = 0;
  *pzbvec = 0;

  Allow = 1;
  pc = getenv("FS_ALLOW_DWI_SIEMENS_ALT");
  if (pc == NULL)
    Allow = 0;
  else if (strcmp(pc, "0") == 0)
    Allow = 0;
  if (!Allow) {
    if (Gdiag_no > 0) {
      printf("Not attempting dcmGetDWIParamsSiemensAlt(), to allow");
      printf("  setenv FS_ALLOW_DWI_SIEMENS_ALT 1\n");
    }
    return (0);
  }

  // Get the nasty string
  Ctx = NULL;
  tag = DCM_MAKETAG(0x29, 0x1010);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  // If string is not there, assume it is not DWI
  if (cond != DCM_NORMAL) {
    free(e);
    return (0);
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  // If string has no value, assume it is not DWI
  if (cond != DCM_NORMAL) {
    free(e);
    return (0);
  }

  // Scroll through the nasty string and find the keywords
  bval_flag = 0;
  bvec_flag = 0;
  for (unsigned int n = 0; n < e->length; n++) {
    c = e->d.string[n];
    if (c == 'B') {
      sscanf(&(e->d.string[n]), "%s", tmpstr);
      if (strcmp(tmpstr, "B_value") == 0) {
        for (m = n; m < e->length; m++) {
          c = e->d.string[m];
          if (isdigit(c)) {
            sscanf(&(e->d.string[m]), "%lf", pbval);
            bval_flag = 1;
            if (Gdiag_no > 0) printf("bval = %lf\n", *pbval);
            break;
          }  // if
        }    // for m
        // Skip past the number just read
        while (isdigit(c) || c == '.') {
          if (m >= e->length) break;
          m = m + 1;
          c = e->d.string[m];
        }
        n = m;
      }  // B_value
    }    // B
    if (c == 'D') {
      sscanf(&(e->d.string[n]), "%s", tmpstr);
      if (strcmp(tmpstr, "DiffusionGradientDirection") == 0) {
        for (k = 0; k < 3; k++) {
          for (m = n; m < e->length; m++) {
            c = e->d.string[m];
            if (isdigit(c) || c == '-' || c == '+' || c == '.') {
              sscanf(&(e->d.string[m]), "%lf", &val);
              if (Gdiag_no > 0) printf("k=%d, val = %lf\n", k, val);
              if (k == 0) *pxbvec = -val;  // convert from LPS to RAS
              if (k == 1) *pybvec = -val;  // convert from LPS to RAS
              if (k == 2) *pzbvec = val;
              bvec_flag++;
              break;
            }  // if
          }    // m
          // Skip past the number just read
          while (isdigit(c) || c == '-' || c == '+' || c == '.') {
            if (m >= e->length) break;
            m = m + 1;
            c = e->d.string[m];
          }
          n = m;
        }  // k
      }    // Gradient
    }      // D
  }        // loop over characters

  rms = sqrt((*pxbvec) * (*pxbvec) + (*pybvec) * (*pybvec) + (*pzbvec) * (*pzbvec));
  if (Gdiag_no > 0) printf("%lf %lf %lf %lf %lf\n", *pbval, *pxbvec, *pybvec, *pzbvec, rms);

  if (!bval_flag || bvec_flag != 3) {
    *pbval = 0;
    *pxbvec = 0;
    *pybvec = 0;
    *pzbvec = 0;
    return (0);
  }

  // Have to check that the values are reasonable because we are just grabbing
  // numbers from a file.
  if (*pbval == 0 || fabs(*pxbvec) > 1.0 || fabs(*pybvec) > 1.0 || fabs(*pzbvec) > 1.0 || fabs(rms - 1) > .001) {
    *pbval = 0;
    *pxbvec = 0;
    *pybvec = 0;
    *pzbvec = 0;
    if (Gdiag_no > 0) {
      printf("%lf %lf %lf %lf %lf\n", *pbval, *pxbvec, *pybvec, *pzbvec, rms);
      printf("These don't look like reasonable DWI params, so I'm setting to 0\n");
    }
  }

  return (0);
}

/*
  \fn int dcmImageDirCosObject()
  \brief Extracts the direction cosine for the image from the DCM object.
  Differs from dcmImageDirCos() in that it uses an object instead of
  opening a file. Returns 1 if error, 0 otherwise.
 */
int dcmImageDirCosObject(DCM_OBJECT *dcm, double *Vcx, double *Vcy, double *Vcz, double *Vrx, double *Vry, double *Vrz)
{
  DCM_ELEMENT *e;
  CONDITION cond;
  DCM_TAG tag;
  unsigned int rtnLength;
  void *Ctx = NULL;
  unsigned int n;
  int nbs;
  char *s;
  double rms;

  /* Load the direction cosines - this is a string of the form:
     Vcx\Vcy\Vcz\Vrx\Vry\Vrz */
  Ctx = NULL;
  tag = DCM_MAKETAG(0x0020, 0x0037);
  e = (DCM_ELEMENT *)calloc(1, sizeof(DCM_ELEMENT));
  cond = DCM_GetElement(&dcm, tag, e);
  if (cond != DCM_NORMAL) {
    free(e);
    return (1);
  }
  AllocElementData(e);
  cond = DCM_GetElementValue(&dcm, e, &rtnLength, &Ctx);
  if (cond != DCM_NORMAL) {
    free(e);
    return (1);
  }
  s = e->d.string;

  /* replace back slashes with spaces */
  nbs = 0;
  for (n = 0; n < strlen(s); n++) {
    if (s[n] == '\\') {
      s[n] = ' ';
      nbs++;
    }
  }
  if (nbs != 5) return (1);

  sscanf(s, "%lf %lf %lf %lf %lf %lf ", Vcx, Vcy, Vcz, Vrx, Vry, Vrz);

  /* Convert Vc from LPS to RAS and Normalize */
  rms = sqrt((*Vcx) * (*Vcx) + (*Vcy) * (*Vcy) + (*Vcz) * (*Vcz));
  (*Vcx) /= -rms;
  (*Vcy) /= -rms;
  (*Vcz) /= +rms;

  /* Convert Vr from LPS to RAS and Normalize */
  rms = sqrt((*Vrx) * (*Vrx) + (*Vry) * (*Vry) + (*Vrz) * (*Vrz));
  (*Vrx) /= -rms;
  (*Vry) /= -rms;
  (*Vrz) /= +rms;

  FreeElementData(e);
  free(e);

  return (0);
}

/*
  \fn MATRIX *ImageDirCos2Slice()
  \brief Computes the slice direction cosine given the image direction
  cosine. Uses a cross-product. The sign of the slice direction cosine
  will be arbitrary.
 */
MATRIX *ImageDirCos2Slice(
    double Vcx, double Vcy, double Vcz, double Vrx, double Vry, double Vrz, double *Vsx, double *Vsy, double *Vsz)
{
  VECTOR *Vc, *Vr, *Vs;
  MATRIX *Mdc;

  Vc = MatrixAlloc(3, 1, MATRIX_REAL);
  Vr = MatrixAlloc(3, 1, MATRIX_REAL);

  Vc->rptr[1][1] = Vcx;
  Vc->rptr[1][2] = Vcy;
  Vc->rptr[1][3] = Vcz;
  Vr->rptr[1][1] = Vrx;
  Vr->rptr[1][2] = Vry;
  Vr->rptr[1][3] = Vrz;

  Vs = VectorCrossProduct(Vc, Vr, NULL);
  *Vsx = Vs->rptr[1][1];
  *Vsy = Vs->rptr[2][1];
  *Vsz = Vs->rptr[3][1];

  MatrixFree(&Vc);
  MatrixFree(&Vr);
  MatrixFree(&Vs);

  Mdc = MatrixAlloc(3, 3, MATRIX_REAL);
  Mdc->rptr[1][1] = Vcx;
  Mdc->rptr[2][1] = Vcy;
  Mdc->rptr[3][1] = Vcz;
  Mdc->rptr[1][2] = Vrx;
  Mdc->rptr[2][2] = Vry;
  Mdc->rptr[3][2] = Vrz;
  Mdc->rptr[1][3] = *Vsx;
  Mdc->rptr[2][3] = *Vsy;
  Mdc->rptr[3][3] = *Vsz;

  return (Mdc);
}


/*!
  \fn int DCMcheckInterceptSlope(DCM_OBJECT *object)
  This prints out a warning if intercept slope are present and not equal to 0,1
*/
int DCMcheckInterceptSlope(DCM_OBJECT *object)
{
  CONDITION cond;
  DCM_TAG tag;
  int ret = 0;
  char *strtmp;

  tag = DCM_MAKETAG(0x28, 0x1052);
  cond = GetString(&object, tag, &strtmp);
  if(cond == DCM_NORMAL) {
    double RescaleIntercept;
    sscanf(strtmp, "%lf", &RescaleIntercept);
    free(strtmp);
    if(RescaleIntercept != 0.0){
      printf("\n\n");
      printf("WARNING: RescaleIntercept = %lf but will not be applied\n",RescaleIntercept);
      printf("\n\n");
    ret = 1;
    }
  }
  tag = DCM_MAKETAG(0x28, 0x1053);
  cond = GetString(&object, tag, &strtmp);
  if(cond == DCM_NORMAL) {
    double RescaleSlope;
    sscanf(strtmp, "%lf", &RescaleSlope);
    free(strtmp);
    if(RescaleSlope != 1.0){
      printf("\n\n");
      printf("WARNING: RescaleSlope = %lf but will not be applied\n",RescaleSlope);
      printf("\n\n");
      ret = 1;
    }
  }
  return(ret);
}
