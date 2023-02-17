/**
 * @brief Identify MRI volume format based on filename extension.
 *
 */
/*
 * Original Author: Christian Haselgrove
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <errno.h>
#include <libgen.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>

#include "mri.h"

#include "Bruker.h"
#include "DICOMRead.h"
#include "analyze.h"
#include "fio.h"
#include "machine.h"
#include "minc.h"
#include "mrisurf.h"
#include "proto.h"
#include "signa.h"
#include "utils.h"

#include "mri_identify.h"

extern int errno;

#ifdef SunOS
int stricmp(char *str1, char *str2);
#endif

char *type_to_string(int type)
{
  char *typestring;
  const char *tmpstr;
  int lentmp;

  switch (type) {
    case MRI_CORONAL_SLICE_DIRECTORY:
      tmpstr = "COR";
      break;
    case MRI_GCA_FILE:
      tmpstr = "GCA";
      break;
    case MRI_MINC_FILE:
      tmpstr = "MINC";
      break;
    case MRI_ANALYZE_FILE:
      tmpstr = "analyze3d";
      break;
    case MRI_ANALYZE4D_FILE:
      tmpstr = "analyze4d";
      break;
    case MRI_MGH_FILE:
      tmpstr = "MGH";
      break;
    case GENESIS_FILE:
      tmpstr = "genesis";
      break;
    case GE_LX_FILE:
      tmpstr = "gelx";
      break;
    case SIEMENS_FILE:
      tmpstr = "siemens";
      break;
    case DICOM_FILE:
      tmpstr = "dicom";
      break;
    case SIEMENS_DICOM_FILE:
      tmpstr = "siemens_dicom";
      break;
    case BRIK_FILE:
      tmpstr = "brik";
      break;
    case BSHORT_FILE:
      tmpstr = "bshort";
      break;
    case BFLOAT_FILE:
      tmpstr = "bfloat";
      break;
    case BHDR:
      tmpstr = "bhdr";
      break;
    case SDT_FILE:
      tmpstr = "varian";
      break;
    case OTL_FILE:
      tmpstr = "outline";
      break;
    case GDF_FILE:
      tmpstr = "gdf";
      break;
    case BRUKER_FILE:
      tmpstr = "bruker";
      break;
    case XIMG_FILE:
      tmpstr = "ximg";
      break;
    case NIFTI1_FILE:
      tmpstr = "nifti1";
      break;
    case NII_FILE:
      tmpstr = "nii";
      break;
    case NRRD_FILE:
      tmpstr = "nrrd";
      break;
    case MRI_CURV_FILE:
      tmpstr = "curv";
      break;
    case GIFTI_FILE:
      tmpstr = "gii";
      break;
    case MGH_ANNOT:
      tmpstr = "annot";
      break;
    case NPY_FILE:
      tmpstr = "npy";
      break;
    default:
      tmpstr = "unknown";
      break;
  }

  lentmp = strlen(tmpstr);
  typestring = (char *)calloc(lentmp + 1, sizeof(char));
  memmove(typestring, tmpstr, lentmp);
  return (typestring);
}

int string_to_type(const char *string)
{
  int type = MRI_VOLUME_TYPE_UNKNOWN;
  char ls[STRLEN];

  // no extension, then return
  if (strlen(string) == 0) return MRI_VOLUME_TYPE_UNKNOWN;

  strcpy(ls, string);
  StrLower(ls);
  // is this compressed?
  if (strcmp(ls, "gz") == 0) type = MRI_GZIPPED;
  else if ((strcmp(ls, "ae") == 0) || (strcmp(ls, "sae") == 0) || (strcmp(ls, "csae") == 0)) type = MGH_AUTOENCODER;
  else if ((strcmp(ls, "m3z") == 0) || (strcmp(ls, "m3d") == 0)) type = MGH_MORPH;
  else if (strcmp(ls, "cor") == 0) type = MRI_CORONAL_SLICE_DIRECTORY;
  else if (strcmp(ls, "minc") == 0 || strcmp(ls, "mnc") == 0) type = MRI_MINC_FILE;
  else if (strcmp(ls, "spm") == 0 || strcmp(ls, "analyze") == 0 || strcmp(ls, "analyze3d") == 0) type = MRI_ANALYZE_FILE;
  else if (strcmp(ls, "analyze4d") == 0 || strcmp(ls, "img") == 0)  // img used to be 3d
    type = MRI_ANALYZE4D_FILE;
  else if (strcmp(ls, "mgh") == 0 || strcmp(ls, "mgz") == 0) type = MRI_MGH_FILE;
  else if (strcmp(ls, "gca") == 0) type = MRI_GCA_FILE;
  else if (strcmp(ls, "signa") == 0) type = SIGNA_FILE;
  else if (strcmp(ls, "ge") == 0 || strcmp(ls, "genesis") == 0) type = GENESIS_FILE;
  else if (strcmp(ls, "gelx") == 0 || strcmp(ls, "lx") == 0) type = GE_LX_FILE;
  else if (strcmp(ls, "bshort") == 0) type = BSHORT_FILE;
  else if (strcmp(ls, "bfloat") == 0) type = BFLOAT_FILE;
  else if (strcmp(ls, "bhdr") == 0) type = BHDR;
  else if (strcmp(ls, "siemens") == 0 || strcmp(ls, "ima") == 0) type = SIEMENS_FILE;
  else if (strcmp(ls, "dicom") == 0) type = DICOM_FILE;
  else if (strcmp(ls, "siemens_dicom") == 0) type = SIEMENS_DICOM_FILE;
  else if (strcmp(ls, "brik") == 0 || strcmp(ls, "afni") == 0) type = BRIK_FILE;
  else if (strcmp(ls, "sdt") == 0 || strcmp(ls, "varian") == 0) type = SDT_FILE;
  else if (strcmp(ls, "otl") == 0 || strcmp(ls, "outline") == 0) type = OTL_FILE;
  else if (strcmp(ls, "gdf") == 0) type = GDF_FILE;
  else if (strcmp(ls, "bruker") == 0) type = BRUKER_FILE;
  else if (strcmp(ls, "ximg") == 0) type = XIMG_FILE;
  else if (strcmp(ls, "nifti1") == 0) type = NIFTI1_FILE;
  else if (strcmp(ls, "nii") == 0) type = NII_FILE;
  else if (strcmp(ls, "nrrd") == 0) type = NRRD_FILE;
  // check for IMAGE file
  else if (!strcmp(ls, "mat") || !strcmp(ls, "tif") || !strcmp(ls, "tiff") || !strcmp(ls, "jpg") || !strcmp(ls, "jpeg") ||
      !strcmp(ls, "pgm") || !strcmp(ls, "ppm") || !strcmp(ls, "pbm") || !strcmp(ls, "rgb"))
    type = IMAGE_FILE;
  else if (strcmp(ls, "curv") == 0) type = MRI_CURV_FILE;
  else if (strcmp(ls, "gii") == 0) type = GIFTI_FILE;
  else if (strcmp(ls, "vtk") == 0) type = VTK_FILE;
  else if (strcmp(ls, "label") == 0) type = MGH_LABEL_FILE;
  else if (strcmp(ls, "annot") == 0) type = MGH_ANNOT;
  else if (strcmp(ls, "npy") == 0) type = NPY_FILE;

  return (type);

} /* end string_to_type() */

// why the routine does not check ANALYZE4D?
int mri_identify(const char *fname_passed)
{
  char fname[STRLEN];
  int type = 0;
  char *ext;

  // Before coming in here, @ and # have been processed
  // remove @ and # strings
  MRIgetVolumeName(fname_passed, fname);

  // check whether this is a directory
  if (!fio_IsDirectory(fname))  // if not a directory, then do extension check
  {
    // now get the extension
    ext = strrchr(fname, '.');
    if (ext) {  // if found a dot (.)
      ++ext;    // now points to extension
      // first use the extension to identify
      type = string_to_type(ext);

      if (type == MRI_GZIPPED) {
        if (strstr(fname, ".mgh.gz")) type = MRI_MGH_FILE;
        if (strstr(fname, ".nii.gz"))
          type = NII_FILE;
        else {
          type = MRI_VOLUME_TYPE_UNKNOWN;
          printf(
              "INFO: Currently supports gzipped mgh or nifti file "
              "(.mgz, .mgh.gz, or nii.gz) only.\n");
          printf("fname = %s\n", fname);
          return type;
        }
      }

      ///////////////////////////////////////////////
      // if type is found then verify
      // IMAGE file uses only extension
      if (type == IMAGE_FILE || type == NPY_FILE) return type;

      if (type != MRI_VOLUME_TYPE_UNKNOWN) {
        switch (type) {
          case MRI_GCA_FILE:
            return (type);
            break;
          case BRUKER_FILE:  // this cannot be identified by extension
            if (is_bruker(fname)) return type;
            break;
          case MRI_CORONAL_SLICE_DIRECTORY:
            if (is_cor(fname)) return type;
            break;
          case BHDR:
            return type;
            break;
          case BSHORT_FILE:
            if (is_bshort(fname)) return type;
            break;
          case BFLOAT_FILE:
            if (is_bfloat(fname)) return type;
            break;
          case SIEMENS_DICOM_FILE:
            if (IsSiemensDICOM(fname)) return type;
            break;
          case DICOM_FILE:
            if (IsDICOM(fname)) return type;
            break;
          case GENESIS_FILE:
            if (is_genesis(fname)) return type;
            break;
          case SIGNA_FILE:
            if (is_signa(fname)) return type;
            break;
          case GE_LX_FILE:
            if (is_ge_lx(fname)) return type;
            break;
          case SDT_FILE:
            if (is_sdt(fname)) return type;
            break;
          case MRI_MGH_FILE:
            if (is_mgh(fname)) return type;
            break;
          case MRI_MINC_FILE:
            if (is_mnc(fname)) return type;
            break;
          case NII_FILE:
            if (is_nii(fname)) return type;
            break;
          case NIFTI1_FILE:
            if (is_nifti1(fname)) return type;
            break;
          case MRI_ANALYZE_FILE:
            // Need to check nifti1 here because it has .img/.hdr
            // like analyze but has a different header
            if (is_nifti1(fname)) return NIFTI1_FILE;
            if (is_analyze(fname)) return type;
            break;
          case MRI_ANALYZE4D_FILE:
            // must add is_analyze4d().  I have no idea what to do thus return
            return type;
            break;
          case SIEMENS_FILE:
            if (is_siemens(fname)) return type;
            break;
          case BRIK_FILE:
            if (is_brik(fname)) return type;
            break;
          case OTL_FILE:
            if (is_otl(fname)) return type;
            break;
          case GDF_FILE:
            if (is_gdf(fname)) return type;
            break;
          case XIMG_FILE:
            if (is_ximg(fname)) return type;
            break;
          case NRRD_FILE:
            if (is_nrrd(fname)) return type;
            break;
          case GIFTI_FILE:
            return type;
            break;
          case MGH_ANNOT:
            return type;
            break;
          case VTK_FILE:
            return type;
            break;
	case MGH_LABEL_FILE:
	  return type ;
	  break ;
          default:
            break;
        }
      }
    }
  }
  //////////////////////////////////////////////////////////////
  // using extension to find type failed or verification failed
  //////////////////////////////////////////////////////////////
  // if type cannot be found then go through the list again
  if (is_bruker(fname))
    return (BRUKER_FILE);
  else if (is_cor(fname))
    return (MRI_CORONAL_SLICE_DIRECTORY);
  else if (is_bshort(fname))
    return (BSHORT_FILE);
  else if (is_bfloat(fname))
    return (BFLOAT_FILE);
  else if (IsSiemensDICOM(fname))
    return (SIEMENS_DICOM_FILE);
  else if (IsDICOM(fname))
    return (DICOM_FILE);
  else if (is_genesis(fname))
    return (GENESIS_FILE);
  else if (is_signa(fname))
    return (SIGNA_FILE);
  else if (is_ge_lx(fname))
    return (GE_LX_FILE);
  else if (is_sdt(fname))
    return (SDT_FILE);
  else if (is_mgh(fname))
    return (MRI_MGH_FILE);
  else if (is_mnc(fname))
    return (MRI_MINC_FILE);
  else if (is_nifti1(fname))  // must appear before ANALYZE
    return (NIFTI1_FILE);
  else if (is_nii(fname))
    return (NII_FILE);
  else if (is_analyze(fname))
    return (MRI_ANALYZE_FILE);
  else if (is_siemens(fname))
    return (SIEMENS_FILE);
  else if (is_brik(fname))
    return (BRIK_FILE);
  else if (is_otl(fname))
    return (OTL_FILE);
  else if (is_gdf(fname))
    return (GDF_FILE);
  else if (is_ximg(fname))
    return (XIMG_FILE);
  else if (is_nrrd(fname))
    return (NRRD_FILE);
  else if (IDisCurv(fname))
    return (MRI_CURV_FILE);
  else if (type == MGH_MORPH)
    return (type);
  else if (type == MGH_AUTOENCODER)
    return (type);
  else
    return (MRI_VOLUME_TYPE_UNKNOWN);
} /*  end mri_identify()  */

/*
  \fn int IDtypeFromStem(char *stem)
  \brief Returns the file type from the stem by looking on disk.
     This can be dependent on order, if more than one file type
     exists with the same stem.
*/
int IDtypeFromStem(const char *stem)
{
  char tmpstr[2000];

  sprintf(tmpstr, "%s.nii", stem);
  if (fio_FileExistsReadable(tmpstr)) return (NII_FILE);

  sprintf(tmpstr, "%s.nii.gz", stem);
  if (fio_FileExistsReadable(tmpstr)) return (NII_FILE);

  sprintf(tmpstr, "%s.mgh", stem);
  if (fio_FileExistsReadable(tmpstr)) return (MRI_MGH_FILE);

  sprintf(tmpstr, "%s.mgz", stem);
  if (fio_FileExistsReadable(tmpstr)) return (MRI_MGH_FILE);

  sprintf(tmpstr, "%s_000.bfloat", stem);
  if (fio_FileExistsReadable(tmpstr)) return (BFLOAT_FILE);

  sprintf(tmpstr, "%s_000.bshort", stem);
  if (fio_FileExistsReadable(tmpstr)) return (BSHORT_FILE);

  return (MRI_VOLUME_TYPE_UNKNOWN);
}

/*
  \fn int IDstemFromName(char *stem)
  \brief Returns the file stem by looking for various extensions.
*/
char *IDstemFromName(const char *name)
{
  char *stem = NULL;
  const char *ext;
  int len;

  len = strlen(name);
  if (len < 3) return (NULL);  // cant be right

  // Try extensions of length 1
  ext = &(name[len - 2]);
  if (!strcmp(ext, ".w")) {
    stem = (char *)calloc(len - 2 + 1, sizeof(char));
    memmove(stem, name, len - 2);
    return (stem);
  }

  // Try extensions of length 3
  if (len < 5) return (NULL);  // cant be right
  ext = &(name[len - 4]);
  if (!strcmp(ext, ".nii") || !strcmp(ext, ".mgz") || !strcmp(ext, ".mgh") || !strcmp(ext, ".img")) {
    stem = (char *)calloc(len - 4 + 1, sizeof(char));
    memmove(stem, name, len - 4);
    return (stem);
  }

  // Try extensions of length 4
  if (len < 6) return (NULL);  // cant be right
  ext = &(name[len - 5]);
  if (!strcmp(ext, ".bhdr")) {
    stem = (char *)calloc(len - 5 + 1, sizeof(char));
    memmove(stem, name, len - 5);
    return (stem);
  }

  // Try extensions of length 6
  if (len < 8) return (NULL);  // cant be right
  ext = &(name[len - 7]);
  if (!strcmp(ext, ".nii.gz")) {
    stem = (char *)calloc(len - 7 + 1, sizeof(char));
    memmove(stem, name, len - 7);
    return (stem);
  }

  // Try _000.bfloat and _000.short
  if (len < 12) return (NULL);  // cant be right
  ext = &(name[len - 11]);
  if (!strcmp(ext, "_000.bfloat") || !strcmp(ext, "_000.bshort")) {
    stem = (char *)calloc(len - 11 + 1, sizeof(char));
    memmove(stem, name, len - 11);
    return (stem);
  }

  // Quietly return
  // printf("ERROR: cannot determine stem from %s\n",name);

  return (NULL);
}

/*
  \fn int IDextensionFromName(char *stem)
  \brief Returns the file extension by looking for various extensions.
*/
char *IDextensionFromName(const char *name)
{
  int len;
  char *ext;

  len = strlen(name);
  if (len < 3) return (NULL);  // cant be right

  ext = (char *)calloc(20, sizeof(char));

  // Does not do .w

  // Try extensions of length 3
  if (len < 5) {free(ext); return (NULL);}  // cant be right
  ext = strncpy(ext, &(name[len - 3]), 3);
  if (!strcmp(ext, "nii") || !strcmp(ext, "mgz") || !strcmp(ext, "mgh") || !strcmp(ext, "img")) return (ext);

  // Try extensions of length 4
  if (len < 6) {free(ext); return (NULL);}  // cant be right
  ext = strncpy(ext, &(name[len - 4]), 4);
  if (!strcmp(ext, "bhdr")) return (ext);

  // Try extensions of length 6
  if (len < 8) {free(ext); return (NULL);}  // cant be right
  ext = strncpy(ext, &(name[len - 6]), 6);
  if (!strcmp(ext, "nii.gz")) return (ext);

  // Try _000.bfloat and _000.short
  if (len < 12) {free(ext); return (NULL);}  // cant be right
  ext = strncpy(ext, &(name[len - 11]), 11);
  if (!strcmp(ext, "_000.bfloat")) {
    sprintf(ext, "bfloat");
    return (ext);
  }
  if (!strcmp(ext, "_000.bshort")) {
    sprintf(ext, "bshort");
    return (ext);
  }

  free(ext);
  return (NULL);
}

/*
  \fn char *IDnameFromStem(char *stem)
  \brief Returns the full file name from the stem by looking on disk.
     This can be dependent on order, if more than one file type
     exists with the same stem.
*/
char *IDnameFromStem(const char *stem)
{
  char tmpstr[2000];

  sprintf(tmpstr, "%s.nii", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  sprintf(tmpstr, "%s.nii.gz", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  sprintf(tmpstr, "%s.mgh", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  sprintf(tmpstr, "%s.mgz", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  sprintf(tmpstr, "%s.img", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  sprintf(tmpstr, "%s.bhdr", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  sprintf(tmpstr, "%s_000.bfloat", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  sprintf(tmpstr, "%s_000.bshort", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  sprintf(tmpstr, "%s.w", stem);
  if (fio_FileExistsReadable(tmpstr)) return (strcpyalloc(tmpstr));

  return (NULL);
}

int is_cor(const char *fname)
{
  struct stat stat_buf;
  char *fname2, *base;
  int iscor = 0;

  if (stat(fname, &stat_buf) < 0) return (0);

  /* if it's a directory, it's a COR dir. */
  if (S_ISDIR(stat_buf.st_mode)) iscor = 1;

  /* if the first four letters are COR- */
  fname2 = strdup(fname);
  base = basename(fname2);
  if (strncmp(base, "COR-", 4) == 0) iscor = 1;

  free(fname2);

  return (iscor);
  ;

} /*  end is_cor()  */

int is_brik(const char *fname)
{
  const char *dot;

  dot = strrchr(fname, '.');
  if (dot) {
    dot++;
    if (!stricmp(dot, "BRIK")) return (1);
  }

  return (0);

} /* end is_brik() */

int is_siemens(const char *fname)
{
  FILE *fp;
  char string[4];
  const char *dot;

  dot = strrchr(fname, '.');
  if (dot) {
    if (!stricmp(dot + 1, "ima")) return (1);
  }

  if ((fp = fopen(fname, "r")) == NULL) {
    errno = 0;
    return (0);
  }

  fseek(fp, 5790, SEEK_SET);
  if (fread(string, 2, 1, fp) < 1) {
    errno = 0;
    fclose(fp);
    return (0);
  }
  string[2] = '\0';
  if (strcmp(string, "SL")) {
    fclose(fp);
    return (0);
  }
  fseek(fp, 5802, SEEK_SET);
  if (fread(string, 2, 1, fp) < 1) {
    errno = 0;
    fclose(fp);
    return (0);
  }
  string[2] = '\0';
  if (strcmp(string, "SP")) {
    fclose(fp);
    return (0);
  }
  fseek(fp, 5838, SEEK_SET);
  if (fread(string, 3, 1, fp) < 1) {
    errno = 0;
    fclose(fp);
    return (0);
  }
  string[3] = '\0';
  if (strcmp(string, "FoV")) {
    fclose(fp);
    return (0);
  }

  fclose(fp);

  return (1);

} /* end is_siemens() */

int is_genesis(const char *fname)
{
  FILE *fp;
  long32 magic;
  const char *dot;

  if (!strncmp(fname, "I.", 2)) return (1);

  dot = strrchr(fname, '.');
  if (dot) {
    if (!strcmp(dot + 1, "MR")) return (1);
  }

  if ((fp = fopen(fname, "r")) == NULL) {
    errno = 0;
    return (0);
  }

  if (fread(&magic, 4, 1, fp) < 1) {
    errno = 0;
    fclose(fp);
    return (0);
  }
  magic = orderLong32Bytes(magic);

  fclose(fp);

  if (magic == GE_MAGIC) return (1);

  return (0);

} /*  end is_genesis()  */

int is_ge_lx(const char *fname)
{
  FILE *fp;
  long32 magic;

  if ((fp = fopen(fname, "r")) == NULL) {
    errno = 0;
    return (0);
  }

  fseek(fp, 3228, SEEK_CUR);
  if (fread(&magic, 4, 1, fp) != 1) {
    errno = 0;
    fclose(fp);
    return (0);
  }
  magic = orderLong32Bytes(magic);

  fclose(fp);

  if (magic == GE_MAGIC) return (1);

  return (0);

} /*  end is_ge_lx()  */

int is_analyze(const char *fname)
{
  FILE *fp;
  dsr hdr;
  char hfname[STRLEN];
  const char *dot;
  long hdr_length;

  strcpy(hfname, fname);

  dot = strrchr(fname, '.');
  if (dot) {
    if (!stricmp(dot + 1, "img")) return (1);
    return (0);
  }

  strcat(hfname, ".hdr");

  if ((fp = fopen(hfname, "r")) == NULL) {
    errno = 0;
    return (0);
  }

  if (fread(&hdr, sizeof(hdr), 1, fp) < 1) {
    errno = 0;
    fclose(fp);
    return (0);
  }

  fseek(fp, 0, SEEK_END);
  hdr_length = ftell(fp);
  fclose(fp);

  if (hdr_length != orderIntBytes(hdr.hk.sizeof_hdr)) return (0);
  if (orderIntBytes(hdr.hk.extents) != 16384) return (0);
  if (hdr.hk.regular != 'r') return (0);

  return (1);

} /*  end is_analyze()  */

int is_mnc(const char *fname)
{
  char buf[3];
  FILE *fp;
  const char *dot;

  dot = strrchr(fname, '.');
  if (dot) {
    if (!stricmp(dot + 1, "mnc")) return (1);

    if (!stricmp(dot + 1, "gz")) {
      /* --- get the next to last dot or the beginning of the file name --- */
      for (dot--; *dot != '.' && dot > fname; dot--)
        ;
      if (stricmp(dot, ".mnc.gz") == 0) return (1);
    }
  }

  if ((fp = fopen(fname, "r")) == NULL) {
    errno = 0;
    return (0);
  }

  if (fread(buf, 1, 3, fp) < 3) {
    errno = 0;
    fclose(fp);
    return (0);
  }

  fclose(fp);

  if (strncmp(buf, "CDF", 3) == 0) return (1);

  return (0);

} /*  end is_mnc()  */

int is_mgh(const char *fname)
{
  FILE *fp;
  int width, height, depth, nframes;
  // int version, type, dof;

  if (strstr(fname, ".mgh") || strstr(fname, ".mgz") || strstr(fname, ".mgh.gz")) return 1;

  if ((fp = fopen(fname, "r")) == NULL) {
    errno = 0;
    return (0);
  }

  // version =
  freadInt(fp);
  width = freadInt(fp);
  height = freadInt(fp);
  depth = freadInt(fp);
  nframes = freadInt(fp);
  // type =
  freadInt(fp);
  // dof =
  freadInt(fp);

  fclose(fp);

  /* my estimates (ch) */
  if (width < 64 || height < 64 || width > 1024 || height > 1024) return (0);
  if (depth > 2000) return (0);
  if (nframes > 2000) return (0);

  return (1);

} /*  end is_mgh()  */
/*--------------------------------------*/
int is_bshort(const char *fname)
{
  const char *dot;
  dot = strrchr(fname, '.');
  if (dot) {
    dot++;
    if (!strcmp(dot, "bshort")) return (1);
  }
  return (0);

} /*  end is_bshort()  */
/*--------------------------------------*/
int is_bfloat(const char *fname)
{
  const char *dot;
  dot = strrchr(fname, '.');
  if (dot) {
    dot++;
    if (!strcmp(dot, "bfloat")) return (1);
  }
  return (0);
} /*  end is_bfloat()  */
/*--------------------------------------*/
int is_bhdr(const char *fname)
{
  const char *dot;
  dot = strrchr(fname, '.');
  if (dot) {
    dot++;
    if (!strcmp(dot, "bhdr")) return (1);
  }
  return (0);
} /*  end is_bhdr()  */
/*--------------------------------------
  bhdr_stem(). Given fname = stem.bhdr,
  returns stem.
--------------------------------------*/
char *bhdr_stem(const char *fname)
{
  char *stem;
  int i, len;

  if (!is_bhdr(fname)) return (NULL);
  len = strlen(fname);
  stem = (char *)calloc(len + 1, sizeof(char));
  memmove(stem, fname, len);
  i = len - 1;
  while (stem[i] != '.') {
    stem[i] = '\0';
    i--;
  }
  stem[i] = '\0';

  return (stem);
}
/*------------------------------------------------
  bhdr_firstslicefname(). Given fname = stem.bhdr,
  finds and returns stem_000.precision on disk,
  where precision is either bfloat or bshort.
  This can then be used as input to MRIread().
  --------------------------------------------*/
char *bhdr_firstslicefname(const char *fname)
{
  char *stem, *firstslicefname;
  int len;

  if (!is_bhdr(fname)) return (NULL);
  stem = bhdr_stem(fname);

  len = strlen(stem) + 12;
  firstslicefname = (char *)calloc(len, sizeof(char));

  sprintf(firstslicefname, "%s_000.bfloat", stem);
  if (fio_FileExistsReadable(firstslicefname)) {
    free(stem);
    return (firstslicefname);
  }

  sprintf(firstslicefname, "%s_000.bshort", stem);
  if (fio_FileExistsReadable(firstslicefname)) {
    free(stem);
    return (firstslicefname);
  }

  free(stem);
  return (NULL);
}
/*-------------------------------------------------------------
  bhdr_precisionstring(). Given fname = stem.bhdr, finds
  stem_000.precision on disk, where precision is either bfloat or
  bshort.
  --------------------------------------------------------------*/
char *bhdr_precisionstring(const char *fname)
{
  char *stem, *precision;
  char tmpstr[2000];

  if (!is_bhdr(fname)) return (NULL);
  stem = bhdr_stem(fname);

  precision = (char *)calloc(7, sizeof(char));

  sprintf(precision, "bfloat");
  sprintf(tmpstr, "%s_000.%s", stem, precision);
  if (fio_FileExistsReadable(tmpstr)) {
    free(stem);
    return (precision);
  }

  sprintf(precision, "bshort");
  sprintf(tmpstr, "%s_000.%s", stem, precision);
  if (fio_FileExistsReadable(tmpstr)) {
    free(stem);
    return (precision);
  }

  free(stem);
  return (NULL);
}

/*----------------------------------------------------------------------
  bhdr_precision(). Given fname = stem.bhdr, finds stem_000.precision
  on disk, where precision is either bfloat or bshort. If bfloat, then
  returns MRI_FLOAT. If bshort, then returns  MRI_SHORT.
  ---------------------------------------------------------------------*/
int bhdr_precision(const char *fname)
{
  char *stem;
  char tmpstr[2000];

  if (!is_bhdr(fname)) return (0);
  stem = bhdr_stem(fname);

  sprintf(tmpstr, "%s_000.bfloat", stem);
  if (fio_FileExistsReadable(tmpstr)) {
    free(stem);
    return (MRI_FLOAT);
  }

  sprintf(tmpstr, "%s_000.bshort", stem);
  if (fio_FileExistsReadable(tmpstr)) {
    free(stem);
    return (MRI_SHORT);
  }

  free(stem);
  return (MRI_VOLUME_TYPE_UNKNOWN);
}

/*--------------------------------------*/
int is_sdt(const char *fname)
{
  char header_fname[STR_LEN];
  char *dot;
  FILE *fp;

  if ((fp = fopen(fname, "r")) == NULL) {
    errno = 0;
    return (0);
  }

  fclose(fp);

  strcpy(header_fname, fname);

  if ((dot = strrchr(header_fname, '.')))
    sprintf(dot + 1, "spr");
  else
    strcat(header_fname, ".spr");

  if ((fp = fopen(header_fname, "r")) == NULL) {
    errno = 0;
    return (0);
  }

  fclose(fp);

  return (1);

} /* end is_sdt() */

int is_gdf(const char *fname)
{
  const char *dot;

  dot = strrchr(fname, '.');
  if (dot != NULL) {
    if (strcmp(dot, ".gdf") == 0) return (TRUE);
  }

  return (FALSE);

} /* end is_gdf() */

int is_otl(const char *fname)
{
  const char *dot;

  dot = strrchr(fname, '.');
  if (dot != NULL) {
    if (strcmp(dot, ".otl") == 0) return (TRUE);
  }

  return (FALSE);

} /* end is_otl() */

int is_ximg(const char *fname) { return (FALSE); } /* end is_ximg() */

int is_nifti1(const char *fname)
{
  char fname_stem[STRLEN];
  std::string hdr_fname;
  char *dot;
  FILE *fp;
  char magic[4];

  // printf("Checking NIFTI1\n");

  strcpy(fname_stem, fname);
  dot = strrchr(fname_stem, '.');
  if (dot != NULL) {
    if (strcmp(dot, ".img") == 0 || strcmp(dot, ".hdr") == 0) *dot = '\0';
  }
  hdr_fname = std::string(fname_stem) + ".hdr";

  fp = fopen(hdr_fname.c_str(), "r");
  if (fp == NULL) {
    errno = 0;
    return (FALSE);
  }

  if (fseek(fp, 344, SEEK_SET) == -1) {
    errno = 0;
    fclose(fp);
    return (FALSE);
  }

  if (fread(magic, 1, 4, fp) != 4) {
    errno = 0;
    fclose(fp);
    return (FALSE);
  }

  fclose(fp);

  // printf("Checking magic number %s %s\n",magic,NIFTI1_MAGIC);

  if (memcmp(magic, NIFTI1_MAGIC, 4) != 0) return (FALSE);

  return (TRUE);

} /*  end is_nifti1()  */

int is_nii(const char *fname)
{
  const char *dot;
  FILE *fp;
  char magic[4];

  dot = strrchr(fname, '.');

  if (dot != NULL) {
    if (strcmp(dot, ".nii") == 0) return (TRUE);
    if (strcmp(dot, ".gz") == 0) {
      if (strstr(fname, ".nii.gz")) return (TRUE);
    }
  }

  fp = fopen(fname, "r");
  if (fp == NULL) {
    errno = 0;
    return (FALSE);
  }

  if (fseek(fp, 344, SEEK_SET) == -1) {
    errno = 0;
    fclose(fp);
    return (FALSE);
  }

  if (fread(magic, 1, 4, fp) != 4) {
    errno = 0;
    fclose(fp);
    return (FALSE);
  }

  fclose(fp);

  if (memcmp(magic, NII_MAGIC, 4) != 0) return (FALSE);

  return (TRUE);

} /*  end is_nii()  */

int is_nrrd(const char *fname)
{
  const char *dot;
  FILE *fp;
  char magic[4];

  // Check that the extension is .nrrd
  dot = strrchr(fname, '.');

  if (dot != NULL) {
    if ((strcmp(dot, ".nrrd") == 0) && ((signed)strlen(fname) == dot - fname + 5)) {
      return (TRUE);
    }
  }
  // TODO: add check for .nhdr, or do that in separate function

  // Check for the Nrrd magic "NRRD", ignoring the next 4 chars
  // which specify the format version. Files w/ correct magic,
  // but w/out correct extension will be recognized.

  // should it read as binary instead?
  fp = fopen(fname, "r");
  if (fp == NULL) {
    errno = 0;
    return (FALSE);
  }

  if (fread(magic, 1, 4, fp) != 4) {
    errno = 0;
    fclose(fp);
    return (FALSE);
  }
  fclose(fp);

  if (memcmp(magic, NRRD_MAGIC, 4) != 0) return (FALSE);

  return (TRUE);

} /*  end is_nrrd()  */

/*----------------------------------------
  c++ version of is_nrrd():
  ifstream file;
  string magic;
  string filename(fname);

  // Check that the extension is .nrrd
  string::size_type nrrdPos = filename.rfind(".nrrd");
  if((nrrdPos != string::npos) && (nrrdPos == filename.length() - 5)){
    return(TRUE);
  }

  // TODO: add check for .nhdr, or do that in separate function

  // Check for the Nrrd magic "NRRD", ignoring the next 4 chars
  // which specify the format version. Files w/ correct magic,
  // but w/out correct extension will be recognized.

  //should it read as binary instead?
  file.open(fname, ifstream::in);
  if(file.fail()){
    errno = 0;
    return(FALSE);
  }
  file >> magic;
  file.close();

  magic = magic.substr(0, 4);
  if(magic.compare(NRRD_MAGIC) != 0){
    errno = 0;
    return(FALSE);
  }

  return(TRUE);

END c++ version of is_nrrd()
 ----------------------------------------*/

/*----------------------------------------
  IDisCurv() - surface curve file format
  ----------------------------------------*/
int IDisCurv(const char *curvfile)
{
  int magno;
  FILE *fp;

  fp = fopen(curvfile, "r");
  if (fp == NULL) return (0);  // return quietly
  fread3(&magno, fp);
  fclose(fp);
  if (magno != NEW_VERSION_MAGIC_NUMBER) return (FALSE);
  return (TRUE);
}

/* EOF */
