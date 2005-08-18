//
// mris_info.cpp
//

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>

extern "C" {
#include "fio.h"
#include "mri.h"
#include "utils.h"
#include "gcsa.h"
#include "colortab.h"
#include "transform.h"
#include "mrisurf.h"

  char *Progname = "mris_info";
}

// copied from mrisurf.c
#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_QUAD_FILE_MAGIC_NUMBER  (-3 & 0x00ffffff)

using namespace std;

int main(int argc, char *argv[])
{
	char ext[STRLEN] ;
  vector<string> type;
  type.push_back("MRIS_BINARY_QUADRANGLE_FILE");
  type.push_back("MRIS_ASCII_TRIANGLE_FILE");
  type.push_back("MRIS_GEO_TRIANGLE_FILE");
  type.push_back("MRIS_TRIANGULAR_SURFACE=MRIS_ICO_SURFACE");
  type.push_back("MRIS_ICO_FILE");
  type.push_back("MRIS_VTK_FILE");
  
  if (argc < 1)
  {
    cout << "Usage: mris_info <surface>" << endl;
    return -1;
  }

	if (!stricmp(FileNameExtension(argv[1], ext), "gcs"))
	{
		GCSA *gcsa = GCSAread(argv[1]) ;
		if (!gcsa)
		{
			cerr << "could not open " << argv[1] << endl;
			return -1;
		}

		printf("GCSA file %s opened\n", argv[1]) ;
		if (gcsa->ct != NULL)
		{
			printf("color table:\n") ;
			CTABprint(stdout, gcsa->ct) ;
		}
		return(0) ;
	}

  MRIS *mris = MRISread(argv[1]);
  if (!mris)
  {
    cerr << "could not open " << argv[1] << endl;
    return -1;
  }

  cout << "SURFACE INFO ================================================== " << endl;
  cout << "type        : " << type[mris->type].c_str() << endl;
  if (mris->type == MRIS_BINARY_QUADRANGLE_FILE)
  {
    FILE *fp = fopen(argv[1], "rb") ;
    int magic;
    fread3(&magic, fp) ;
    if (magic == QUAD_FILE_MAGIC_NUMBER) 
    {
      cout << "              QUAD_FILE_MAGIC_NUMBER" << endl;
    }
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) 
    {
      cout << "              NEW_QUAD_FILE_MAGIC_NUMBER" << endl;
    }
    fclose(fp);
  }
  cout << "num vertices: " << mris->nvertices << endl;
  cout << "num faces   : " << mris->nfaces << endl;
  cout << "num stripgs : " << mris->nstrips << endl;
  cout << "ctr         : (" << mris->xctr << ", " << mris->yctr << ", " << mris->zctr << ")" << endl;
  cout << "vertex locs : " << (mris->useRealRAS ? "scannerRAS" : "surfaceRAS") << endl;
  if (mris->lta)
  {
    cout << "talairch.xfm: " << endl;
    MatrixPrint(stdout, mris->lta->xforms[0].m_L);
    cout << "surfaceRAS to talaraiched surfaceRAS: " << endl;
    MatrixPrint(stdout, mris->SRASToTalSRAS_);
    cout << "talairached surfaceRAS to surfaceRAS: " << endl;
    MatrixPrint(stdout, mris->TalSRASToSRAS_);
  }
  vg_print(&mris->vg); 

	{
		int i ;
		for (i = 0 ; i < mris->ncmds ; i++)
			printf("cmd[%d]: %s\n", i, mris->cmdlines[i]) ;
	}
  MRISfree(&mris);
}
