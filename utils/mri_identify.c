#include <stdio.h>
#include <string.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>
#include <libgen.h>
#include "mri.h"
#include "proto.h"
#include "mri_identify.h"
#include "analyze.h"
#include "volume_io.h"
#include "machine.h"
#include "signa.h"
#include "fio.h"
#include "DICOMRead.h"
#include "Bruker.h"

extern int errno;

#ifdef SunOS
int stricmp(char *str1, char *str2) ;
#endif

char *type_to_string(int type)
{
  char *typestring;
  char *tmpstr;
  int lentmp;

  switch(type){
  case MRI_CORONAL_SLICE_DIRECTORY: tmpstr = "COR"; break;
  case MRI_MINC_FILE:      tmpstr = "MINC"; break;
  case MRI_ANALYZE_FILE:   tmpstr = "analyze3d"; break;
  case MRI_ANALYZE4D_FILE: tmpstr = "analyze4d"; break;
  case MRI_MGH_FILE: tmpstr = "MGH"; break;
  case GENESIS_FILE: tmpstr = "genesis"; break;
  case GE_LX_FILE:   tmpstr = "gelx"; break;
  case SIEMENS_FILE: tmpstr = "siemens"; break;
  case DICOM_FILE:   tmpstr = "dicom"; break;
  case SIEMENS_DICOM_FILE: tmpstr = "siemens_dicom"; break;
  case BRIK_FILE:   tmpstr = "brik"; break;
  case BSHORT_FILE: tmpstr = "bshort"; break;
  case BFLOAT_FILE: tmpstr = "bfloat"; break;
  case SDT_FILE:    tmpstr = "varian"; break;
  case OTL_FILE:    tmpstr = "outline"; break;
  case GDF_FILE:    tmpstr = "gdf"; break;
  case BRUKER_FILE: tmpstr = "bruker"; break;
  default: tmpstr = "unknown"; break;
  }

  lentmp = strlen(tmpstr);
  typestring = (char *)calloc(lentmp+1,sizeof(char));
  memcpy(typestring,tmpstr,lentmp);
  return(typestring);
}


int string_to_type(char *string)
{

  int type = MRI_VOLUME_TYPE_UNKNOWN;
  char ls[STRLEN];

  strcpy(ls, string);
  StrLower(ls);

  if(strcmp(ls, "cor") == 0)
    type = MRI_CORONAL_SLICE_DIRECTORY;
  if(strcmp(ls, "minc") == 0 || strcmp(ls, "mnc") == 0)
    type = MRI_MINC_FILE;
  if(strcmp(ls, "spm") == 0 || strcmp(ls, "analyze") == 0 ||
     strcmp(ls, "analyze3d") == 0)
    type = MRI_ANALYZE_FILE;
  if(strcmp(ls, "analyze4d") == 0)
    type = MRI_ANALYZE4D_FILE;
  if(strcmp(ls, "mgh") == 0)
    type = MRI_MGH_FILE;
  if(strcmp(ls, "signa") == 0)
    type = SIGNA_FILE;
  if(strcmp(ls, "ge") == 0 || strcmp(ls, "genesis") == 0)
    type = GENESIS_FILE;
  if(strcmp(ls, "gelx") == 0 || strcmp(ls, "lx") == 0)
    type = GE_LX_FILE;
  if(strcmp(ls, "bshort") == 0)
    type = BSHORT_FILE;
  if(strcmp(ls, "bfloat") == 0)
    type = BFLOAT_FILE;
  if(strcmp(ls, "siemens") == 0 || strcmp(ls, "ima") == 0)
    type = SIEMENS_FILE;
  if(strcmp(ls, "dicom") == 0)
    type = DICOM_FILE;
  if(strcmp(ls, "siemens_dicom") == 0)
    type = SIEMENS_DICOM_FILE;
  if(strcmp(ls, "brik") == 0 || strcmp(ls, "afni") == 0)
    type = BRIK_FILE;
  if(strcmp(ls, "sdt") == 0 || strcmp(ls, "varian") == 0)
    type = SDT_FILE;
  if(strcmp(ls, "otl") == 0 || strcmp(ls, "outline") == 0)
    type = OTL_FILE;
  if(strcmp(ls, "gdf") == 0)
    type = GDF_FILE;
  if(strcmp(ls, "bruker")==0)
    type = BRUKER_FILE;
  return(type);

} /* end string_to_type() */

int mri_identify(char *fname_passed)
{

  char fname[STRLEN];

  MRIgetVolumeName(fname_passed, fname);

  if (is_bruker(fname))
    return(BRUKER_FILE);
  else if(is_cor(fname))
    return(MRI_CORONAL_SLICE_DIRECTORY);
  else if(is_bshort(fname))
    return(BSHORT_FILE);
  else if(is_bfloat(fname))
    return(BFLOAT_FILE);
  else if (IsSiemensDICOM(fname))
    return(SIEMENS_DICOM_FILE);
  else if (IsDICOM(fname))
    return(DICOM_FILE);
  else if(is_genesis(fname))
    return(GENESIS_FILE);
  else if(is_signa(fname))
    return(SIGNA_FILE);
  else if(is_ge_lx(fname))
    return(GE_LX_FILE);
  else if(is_sdt(fname))
    return(SDT_FILE);
  else if(is_mgh(fname))
    return(MRI_MGH_FILE);
  else if(is_mnc(fname))
    return(MRI_MINC_FILE);
  else if(is_analyze(fname))
    return(MRI_ANALYZE_FILE);
  else if(is_siemens(fname))
    return(SIEMENS_FILE);
  else if(is_brik(fname))
    return(BRIK_FILE);
  else if(is_otl(fname))
    return(OTL_FILE);
  else if(is_gdf(fname))
    return(GDF_FILE);
  else
    return(MRI_VOLUME_TYPE_UNKNOWN);

}  /*  end mri_identify()  */

int is_cor(char *fname)
{

  struct stat stat_buf;
  char *fname2, *base;
  int iscor = 0;

  if(stat(fname, &stat_buf) < 0)
    return(0);

  /* if it's a directory, it's a COR dir. */
  if(S_ISDIR(stat_buf.st_mode))
    iscor = 1;

  /* if the first four letters are COR- */
  fname2 = strdup(fname);
  base = basename(fname2);
  if(strncmp(base,"COR-",4) == 0)
    iscor = 1;

  free(fname2);

  return(iscor);;

}  /*  end is_cor()  */

int is_bruker(char *fname)
{
  struct stat stat_buf;
  char methodFile[512];
  char acqpFile[512];
  char dataFile[512];
  char d3procFile[512];

  if(stat(fname, &stat_buf) < 0)
    return(0);

  /* if it's a directory, it's a COR dir. */
  if(!S_ISDIR(stat_buf.st_mode))
    return 0;
  // must check all these files exist or not
  return checkBrukerFiles(fname, methodFile, acqpFile, dataFile, d3procFile, 0);
}

int is_brik(char *fname)
{

  char *dot;

  dot = strrchr(fname, '.');
  if(dot)
  {
    dot++;
    if(!stricmp(dot, "BRIK"))
      return(1);
  }

  return(0);

} /* end is_brik() */

int is_siemens(char *fname)
{

  FILE *fp;
  char string[4];
  char *dot;

  dot = strrchr(fname, '.');
  if(dot)
  {
    if(!stricmp(dot+1, "ima"))
      return(1);
  }

  if((fp = fopen(fname, "r")) == NULL)
  {
    errno = 0;
    return(0);
  }

  fseek(fp, 5790, SEEK_SET);
  if(fread(string, 2, 1, fp) < 1)
  {
    errno = 0;
    fclose(fp);
    return(0);
  }
  string[2] = '\0';
  if(strcmp(string, "SL"))
    {
    fclose(fp);
    return(0);
    }
  fseek(fp, 5802, SEEK_SET);
  if(fread(string, 2, 1, fp) < 1)
  {
    errno = 0;
    fclose(fp);
    return(0);
  }
  string[2] = '\0';
  if(strcmp(string, "SP"))
    {
    fclose(fp);
    return(0);
    }
  fseek(fp, 5838, SEEK_SET);
  if(fread(string, 3, 1, fp) < 1)
  {
    errno = 0;
    fclose(fp);
    return(0);
  }
  string[3] = '\0';
  if(strcmp(string, "FoV"))
    {
    fclose(fp);
    return(0);
    }

  fclose(fp);

  return(1);

} /* end is_siemens() */

int is_genesis(char *fname)
{

  FILE *fp;
  long magic;
  char *dot;

  if(!strncmp(fname, "I.", 2))
    return(1);

  dot = strrchr(fname, '.');
  if(dot)
  {
    if(!strcmp(dot+1, "MR"))
      return(1);
  }

  if((fp = fopen(fname, "r")) == NULL)
  {
    errno = 0;
    return(0);
  }

  if(fread(&magic, 4, 1, fp) < 1)
  {
    errno = 0;
    fclose(fp);
    return(0);
  }
  magic = orderLongBytes(magic);

  fclose(fp);

  if(magic == GE_MAGIC)
    return(1);

  return(0);

}  /*  end is_genesis()  */

int is_ge_lx(char *fname)
{

  FILE *fp;
  long magic;

  if((fp = fopen(fname, "r")) == NULL)
  {
    errno = 0;
    return(0);
  }

  fseek(fp, 3228, SEEK_CUR);
  if(fread(&magic, 4, 1, fp) < 0)
  {
    errno = 0;
    fclose(fp);
    return(0);
  }
  magic = orderLongBytes(magic);

  fclose(fp);

  if(magic == GE_MAGIC)
    return(1);

  return(0);

}  /*  end is_ge_lx()  */

int is_analyze(char *fname)
{

  FILE *fp;
  dsr hdr;
  char hfname[STRLEN], *dot;
  long hdr_length;

  strcpy(hfname, fname);

  dot = strrchr(fname, '.') ;
  if (dot)
  {
    if (!stricmp(dot+1, "img"))
      return(1) ;
    return(0);
  }

  dot = '\0';
  sprintf(hfname, "%s.hdr", hfname);

  if((fp = fopen(hfname, "r")) == NULL)
  {
    errno = 0;
    return(0);
  }

  if(fread(&hdr, sizeof(hdr), 1, fp) < 1)
  {
    errno = 0;
    fclose(fp);
    return(0);
  }

  fseek(fp, 0, SEEK_END);
  hdr_length = ftell(fp);
  fclose(fp);

  if(hdr_length != orderIntBytes(hdr.hk.sizeof_hdr))
    return(0);
  if(orderIntBytes(hdr.hk.extents) != 16384)
    return(0);
  if(hdr.hk.regular != 'r')
    return(0);

  return(1);

}  /*  end is_analyze()  */

int is_mnc(char *fname)
{

  char buf[3];
  FILE *fp;
  char *dot;

  dot = strrchr(fname, '.') ;
  if (dot)
  {

    if(!stricmp(dot+1, "mnc"))
      return(1) ;

    if(!stricmp(dot+1, "gz"))
    {
      /* --- get the next to last dot or the beginning of the file name --- */
      for(dot--;*dot != '.' && dot > fname;dot--);
      if(stricmp(dot, ".mnc.gz") == 0)
        return(1);
    }

  }

  if((fp = fopen(fname, "r")) == NULL)
  {
    errno = 0;
    return(0);
  }

  if(fread(buf, 1, 3, fp) < 3)
  {
    errno = 0;
    fclose(fp);
    return(0);
  }

  fclose(fp);

  if(strncmp(buf, "CDF", 3) == 0)
    return(1);

  return(0);

}  /*  end is_mnc()  */

int is_mgh(char *fname)
{
  FILE *fp;
  int version, width, height, depth, nframes, type, dof;
  char *dot ;

  dot = strrchr(fname, '.') ;
  if (dot)
  {
    if (!stricmp(dot+1, "mgh"))
      return(1) ;
  }

  if((fp = fopen(fname, "r")) == NULL)
  {
    errno = 0;
    return(0);
  }

  version = freadInt(fp) ;
  width = freadInt(fp) ;
  height = freadInt(fp) ;
  depth =  freadInt(fp) ;
  nframes = freadInt(fp) ;
  type = freadInt(fp) ;
  dof = freadInt(fp) ;

  fclose(fp);

/* my estimates (ch) */
  if(width < 64 || height < 64 || width > 1024 || height > 1024)
    return(0);
  if(depth > 2000)
    return(0);
  if(nframes > 2000)
    return(0);

  return(1);

}  /*  end is_mgh()  */

int is_bshort(char *fname)
{

  char *dot;

  dot = strrchr(fname, '.');
  if(dot)
  {
    dot++;
    if(!strcmp(dot, "bshort"))
      return(1);
  }

  return(0);

}  /*  end is_bshort()  */

int is_bfloat(char *fname)
{

  char *dot;

  dot = strrchr(fname, '.');
  if(dot)
  {
    dot++;
    if(!strcmp(dot, "bfloat"))
      return(1);
  }

  return(0);

}  /*  end is_bfloat()  */

int is_sdt(char *fname)
{

  char header_fname[STR_LEN];
  char *dot;
  FILE *fp;

  if((fp = fopen(fname, "r")) == NULL)
  {
    errno = 0;
    return(0);
  }

  fclose(fp);

  strcpy(header_fname, fname);

  if((dot = strrchr(header_fname, '.')))
    sprintf(dot+1, "spr");
  else
    strcat(header_fname, ".spr");

  if((fp = fopen(header_fname, "r")) == NULL)
  {
    errno = 0;
    return(0);
  }

  fclose(fp);

  return(1);

} /* end is_sdt() */

int is_gdf(char *fname)
{

  char *dot;

  dot = strrchr(fname, '.');
  if(dot != NULL)
  {
    if(strcmp(dot, ".gdf") == 0)
      return(TRUE);
  }

  return(FALSE);

} /* end is_gdf() */

int is_otl(char *fname)
{

  char *dot;

  dot = strrchr(fname, '.');
  if(dot != NULL)
  {
    if(strcmp(dot, ".otl") == 0)
      return(TRUE);
  }

  return(FALSE);

} /* end is_otl() */

/* EOF */
