#include <stdio.h>
#include <stdlib.h>
#include <sys/utsname.h>
#include <sys/stat.h>

#include "kvlAtlasMesh.h"
#include "kvlAtlasMeshCollection.h"

#include "version.h"
#include "cmdargs.h"

#include "mrisurf.h"
#include "matrix.h"
#include "mri.h"

static int  parse_commandline(int argc, char **argv);
static void check_options();
static void print_usage();
static void usage_exit();
static void print_help();
static void print_version();
static void dump_options();

static int  __mesh2Surf(const char *mesh, const MRI *niiVol, LTA *lta, const char *oSurf, const char *overlayMri);
static void __getTemplateVolgeom(const MRI *niiVol, LTA *lta, int invert, VOL_GEOM *template_vg);

struct utsname uts;
char *cmdline, cwd[2000];

int checkoptsonly = 0;
const char *Progname = NULL;
const char *atlasMesh = NULL;
const char *templateVol = NULL;
const char *outSurf = NULL;
const char *priorsmgz = NULL;
const char *ltaFile = NULL;
int invertLTA = 0;

/*
NAME
	samsegmesh2surf

DESCRIPTION
	This program can be used to generate freesurfer surface from a samseg 
	atlas mesh file. It can also generate priors at each vertex as overlay
	MRI volume (nVertices x 1 x 1 x nClasses).
	       

REQUIRED ARGUMENTS
	Exactly one input is required:

	--atlasmesh    atlas-mesh-collection-file
		input samseg atlas mesh collection file

	At least one input is required:

	--template     atlas-template-volume
		input atlas template volume

	--lta          lta-transform
		input LTA transform to be applied to surface. If both --lta 
		and --template are specified, automatic detection of which 
		direction the LTA goes by looking at which side matches the 
		atlas template volume geomery. Automatically invert if 
		necessary. --invert will not be applied in this case.

	At least one output is required:

	--osurf        output-freesufer-surface
		output freesurfer surface

	--opriors output-priors-as-mri.mgz
		output priors as mri volume

OPTIONAL ARGUMENTS
	--invert
		inverts LTA transform

EXAMPLE 1
	use given template volume as source image, output surface is aligned 
	with template volume:
	
	  samsegmesh2surf 
	    --atlasmesh atlas_level1.txt.gz 
	    --template template.nii 
	    --osurf out.surf 
	    --opriors priors.mgz
	 

EXAMPLE 2
	use LTA src volume as source image, apply the LTA matrix to align 
	output surface with LTA dst volume:
	
	  samsegmesh2surf 
	    --atlasmesh atlas_level1.txt.gz 
	    --lta template.lta
	    --osurf out.surf  
	    --opriors priors.mgz
	 

EXAMPLE 3
	invert LTA, use LTA dst volume as source image, apply the LTA matrix 
	to align output surface with LTA src volume:
	
	  samsegmesh2surf 
	    --atlasmesh atlas_level1.txt.gz 
	    --lta template.lta
	    --invert
	    --osurf out.surf  
	    --opriors priors.mgz
 */
int main(int argc, char** argv)
{
  int nargs;

  nargs = handleVersionOption(argc, argv, "samsegmesh2surf");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);

  /* read template, calculate vox2tkras matrix
   * /usr/local/freesurfer/dev/average/samseg/20Subjects_smoothing2_down2_smoothingForAffine2/template.nii
   * dimensions: 181 x 217 x 181
   *
   * >/usr/local/freesurfer/dev/bin/mri_info template.nii --vox2ras
   *   -1.00000   -0.00000    0.00000   90.00000 
   *   -0.00000    1.00000   -0.00000 -126.00000 
   *    0.00000    0.00000    1.00000  -72.00000 
   *    0.00000    0.00000    0.00000    1.00000 
   * >/usr/local/freesurfer/dev/bin/mri_info template.nii --vox2ras-tkr
   *   -1.00000    0.00000    0.00000   90.50000 
   *    0.00000    0.00000    1.00000  -90.50000 
   *    0.00000   -1.00000    0.00000  108.50000 
   *    0.00000    0.00000    0.00000    1.00000 
   */
  MRI *niiVol = NULL;
  if (templateVol != NULL)
    niiVol = MRIread(templateVol);

  LTA *lta = NULL;
  if (ltaFile != NULL)
    lta = LTAread(ltaFile);

  __mesh2Surf(atlasMesh, niiVol, lta, outSurf, priorsmgz);

  return 0;
}


/* --------------------------------------------- */
static int __mesh2Surf(const char *meshfile, const MRI *niiVol, LTA *lta, const char *oSurf, const char *overlayMri)
{
  MRIS *Surf = NULL;
  MRI  *Vol  = NULL;
  MATRIX *vox2tkras = NULL;

  // read mesh collection
  struct stat buffer;   
  if (stat(meshfile, &buffer) != 0) // error
  {
    /* perror automatically print error message to stderr, use strerror() instead
       char errmsg[256] = {'\0'};
       perror(errmsg);
    */
    printf("Error: failed to get %s status - %s\n", meshfile, strerror(errno));
    return -1;
  }

  kvl::AtlasMeshCollection::Pointer  meshStartCollection = kvl::AtlasMeshCollection::New();
  if (!meshStartCollection->Read(meshfile))
  {
    printf("ERROR: Couldn't read atlas mesh %s\n", meshfile);
    return -1;
  } 
  else
  {
    printf("reading atlas mesh %s\n", meshfile);
  }

  kvl::AtlasMesh::ConstPointer constMesh = meshStartCollection->GetReferenceMesh();
  const int numClasses = constMesh->GetPointData()->Begin().Value().m_Alphas.Size();
  const int numMeshPoints = constMesh->GetPoints()->Size();

  // count triangle faces
  kvl::AtlasMesh::Pointer mutableMesh = const_cast< kvl::AtlasMesh* >( constMesh.GetPointer() );
  kvl::AtlasMesh::CellsContainer::ConstIterator cellsIt    = constMesh->GetCells()->Begin();
  kvl::AtlasMesh::CellsContainer::ConstIterator cellsItEnd = constMesh->GetCells()->End();
  int num_faces = 0;
  for ( ;  cellsIt != cellsItEnd; ++cellsIt )
    if ( cellsIt.Value()->GetType() == kvl::AtlasMesh::CellType::TRIANGLE_CELL )
      num_faces++;

  // the surface codes are borrowed from mrisReadGIFTIdanum()
  // create MRIS structure with given number of vertices and faces
  const int num_vertices = numMeshPoints;  //constMesh->GetPoints()->Size();
  printf("num_vertices = %d, num_faces = %d\n", num_vertices, num_faces);

  if (oSurf != NULL)
  {
    Surf = MRISalloc(num_vertices, num_faces);  //MRISoverAlloc(nVertices, 0, nVertices, 0);
    Surf->type = MRIS_TRIANGULAR_SURFACE; // defined in mrisurf.h
    Surf->hemisphere = NO_HEMISPHERE;
    Surf->useRealRAS = 0;

    VOL_GEOM template_vg;
    __getTemplateVolgeom(niiVol, lta, invertLTA, &template_vg);

    vox2tkras = template_vg.get_Vox2TkregRAS();

    /****** begin setting surface volume geometry ******/
    // Surf->vg & Surf->xctr, Surf->yctr, Surf->zctr
    Surf->vg.width  = template_vg.width;
    Surf->vg.height = template_vg.height;
    Surf->vg.depth  = template_vg.depth;
    Surf->vg.xsize  = template_vg.xsize; Surf->vg.ysize = template_vg.ysize; Surf->vg.zsize = template_vg.zsize;
    Surf->vg.x_r = template_vg.x_r; Surf->vg.x_a = template_vg.x_a; Surf->vg.x_s = template_vg.x_s;
    Surf->vg.y_r = template_vg.y_r; Surf->vg.y_a = template_vg.y_a; Surf->vg.y_s = template_vg.y_s;
    Surf->vg.z_r = template_vg.z_r; Surf->vg.z_a = template_vg.z_a; Surf->vg.z_s = template_vg.z_s;
    Surf->vg.c_r = template_vg.c_r; Surf->vg.c_a = template_vg.c_a; Surf->vg.c_s = template_vg.c_s;

    Surf->vg.valid = 1;
    /****** end setting surface volume geometry ******/
  }

  // create MRI structure for overlay
  if (overlayMri != NULL)
    Vol = new MRI({numMeshPoints, 1, 1, numClasses}, MRI_FLOAT);

  kvl::AtlasMesh::PointsContainer::Iterator pointsIt    = mutableMesh->GetPoints()->Begin(); 
  kvl::AtlasMesh::PointsContainer::Iterator pointsItEnd = mutableMesh->GetPoints()->End();
  for (; pointsIt != pointsItEnd; pointsIt++)
  {
    kvl::AtlasMesh::PointIdentifier pointid   = pointsIt.Index();

    if (Surf != NULL)
    {
      // (c, r, s)
      kvl::AtlasMesh::PointType       thispoint = pointsIt.Value();
      int c = floor(thispoint[0]);
      int r = floor(thispoint[1]);
      int s = floor(thispoint[2]); 

      MATRIX *CRS = MatrixAlloc(4, 1, MATRIX_REAL);
      CRS->rptr[1][1] = c;
      CRS->rptr[2][1] = r;
      CRS->rptr[3][1] = s;
      CRS->rptr[4][1] = 1;

      // Convert the CRS to tkRAS
      MATRIX *tkRAS = MatrixAlloc(4, 1, MATRIX_REAL);
      tkRAS->rptr[4][1] = 1;
      tkRAS = MatrixMultiply(vox2tkras, CRS, tkRAS);

      // set surface vertex xyz
      MRISsetXYZ(Surf, pointid, tkRAS->rptr[1][1], tkRAS->rptr[2][1], tkRAS->rptr[3][1]);
      Surf->vertices[pointid].origarea = -1;
    }

    // assign alphas to overlay volume
    if (Vol != NULL)
    {
      for (int nClass = 0; nClass < numClasses; nClass++)
      {
        float alpha = constMesh->GetPointData()->ElementAt(pointid).m_Alphas[nClass];
        MRIsetVoxVal(Vol, pointid, 0, 0, nClass, alpha);
      }
    }
  }  

  if (Vol != NULL)
  {
    printf("Writing Overlay Volume to %s ...\n", overlayMri);
    int err = MRIwrite(Vol, overlayMri);
    if (err)
    {
      printf("ERROR outputing volume %s\n", overlayMri);
      return 1;
    }
  }

  if (Surf == NULL)
    return 0;

  mrisComputeSurfaceDimensions(Surf);

  cellsIt = mutableMesh->GetCells()->Begin();
  for (int face_index = 0; cellsIt != cellsItEnd; ++cellsIt)
  {
    if ( cellsIt.Value()->GetType() == kvl::AtlasMesh::CellType::TRIANGLE_CELL )
    {
      kvl::AtlasMesh::CellType::PointIdIterator  pit = cellsIt.Value()->PointIdsBegin();
      kvl::AtlasMesh::PointIdentifier  pointId0 = *pit;
      pit++;
      kvl::AtlasMesh::PointIdentifier  pointId1 = *pit;
      pit++;
      kvl::AtlasMesh::PointIdentifier  pointId2 = *pit;

      Surf->faces[face_index].v[0] = pointId0;
      Surf->vertices_topology[pointId0].num++;
      Surf->faces[face_index].v[1] = pointId1;
      Surf->vertices_topology[pointId1].num++;
      Surf->faces[face_index].v[2] = pointId2;
      Surf->vertices_topology[pointId2].num++;

      face_index++;
    } // if (kvl::AtlasMesh::CellType::TRIANGLE_CELL)
  }

  // each vertex has a face list (faster than face list in some operations)
  for (int vertex_index = 0; vertex_index < num_vertices; vertex_index++)
  {
    Surf->vertices_topology[vertex_index].f = (int   *)calloc(Surf->vertices_topology[vertex_index].num, sizeof(int));
    Surf->vertices_topology[vertex_index].n = (uchar *)calloc(Surf->vertices_topology[vertex_index].num, sizeof(uchar));
    Surf->vertices_topology[vertex_index].num = 0;  // this gets re-calc'd next...
  }

  for (int face_index = 0; face_index < Surf->nfaces; face_index++)
  {
    FACE *face = &Surf->faces[face_index];
    int n;
    for (n = 0; n < VERTICES_PER_FACE; n++)
      Surf->vertices_topology[face->v[n]].f[Surf->vertices_topology[face->v[n]].num++] = face_index;  // note that .num is auto-incremented!
  }

  for (int vertex_index = 0; vertex_index < num_vertices; vertex_index++) {
    int n, m;
    for (n = 0; n < Surf->vertices_topology[vertex_index].num; n++) {
      for (m = 0; m < VERTICES_PER_FACE; m++) {
        if (Surf->faces[Surf->vertices_topology[vertex_index].f[n]].v[m] == vertex_index) {
          Surf->vertices_topology[vertex_index].n[n] = m;
        }
      }
    }
  }

  mrisCompleteTopology(Surf);

  MRIScomputeNormals(Surf);

  UpdateMRIS(Surf, NULL);

  // apply LTA transform if --lta <> is specified
  if (lta != NULL)
  {
    /* 
     * It is commented in mrisurf_metricProperties.cpp::MRISltaMultiply():
     *
     * 1. Determine which direction the LTA goes by looking at which side
     *    matches the surface volume geometry. Invert if necessary
     * 2. lta needs to be in regdat space to apply it to the surface
     * 3. lta is changed to type REGISTER_DAT before applying.
     *    In REGISTER_DAT format, the matrix actually goes from target/dst
     *    to mov/src so have to invert the matrix.  Note can't use
     *    LTAinvert() because it will also reverse the src and dst vol geometries.
     * 4. multiply the surf coords by the matrix
     * 5. Reverse the faces if the reg has neg determinant
     * 6. Copy the volume geometry of the destination volume
     *
     */
    int err = MRISltaMultiply(Surf, lta);
    if (err)
    {
      printf("ERROR apply LTA transform to surface\n");
      return 1;
    }
  }

  printf("Writing Surface to %s ...\n", oSurf);
  int err = MRISwrite(Surf, oSurf);
  if (err)
  {
    printf("ERROR outputing surface %s\n", oSurf);
    return 1;
  }

  return 0;
}


/* --------------------------------------------- */
static void __getTemplateVolgeom(const MRI *niiVol, LTA *lta, int invert, VOL_GEOM *template_vg)  
{
  if (niiVol != NULL)
  {
    // use template as source image
    //copyVolGeom(niiVol, template_vg);
    // niiVol is MRI*, this will copy only vol geometry
    *template_vg = *niiVol;
  }
  else
  {
    if (invert)
    {
      printf("INFO: invert LTA\n");
      LTAinvert(lta, lta);
    }

    // get it from LTA src volume geom
    printf("INFO: Use LTA src volume geometry\n");

    //copyVolGeom(&(lta->xforms[0].src), template_vg);
    *template_vg = lta->xforms[0].src;
  }
}


/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))             print_help() ;
    else if (!strcasecmp(option, "--version"))     print_version() ;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--atlasmesh"))
    {
      if (nargc < 1) CMDargNErr(option,1);
      atlasMesh = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--template"))
    {
      if (nargc < 1) CMDargNErr(option,1);
      templateVol = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--osurf"))
    {
      if (nargc < 1) CMDargNErr(option,1);
      outSurf = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--opriors"))
    {
      if (nargc < 1) CMDargNErr(option,1);
      priorsmgz = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--lta"))
    {
      if (nargc < 1) CMDargNErr(option,1);
      ltaFile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--invert"))
    {
      invertLTA = 1;
    }
    else 
    {
      printf("ERROR: Option %s unknown\n", option);
      if (CMDsingleDash(option))
        printf("       Did you really mean -%s ?\n\n", option);
      print_help();
      exit(1);
    }

    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* --------------------------------------------- */
static void check_options()
{
  dump_options();

  // atlasMesh is required
  if (atlasMesh == NULL)
  {
    printf("ERROR: input samseg atlas mesh is required: --atlasmesh <atlas-mesh-collection-file>\n");
    exit(1);
  }

  // either outSurf or priorsmgz needs to be specified
  if (outSurf == NULL && priorsmgz == NULL)
  {
    printf("ERROR: At least one output is required: --osurf <output-freesufer-surface> --opriors <output-priors-as-mri.mgz>\n");
    exit(1);
  }

  //
  if (templateVol == NULL && ltaFile == NULL)
  {
    printf("ERROR: At least one input is required: --template <atlas-template-volume> --lta <lta-transform>\n");
    exit(1);   
  }

  return;
}
/* ------------------------------------------------------ */
#include "samsegmesh2surf.help.xml.h"
static void print_usage()
{
  outputHelpXml(samsegmesh2surf_help_xml, samsegmesh2surf_help_xml_len);
}
/* ------------------------------------------------------ */
static void print_help()
{
  print_usage();
  exit(1) ;
}
/* ------------------------------------------------------ */
static void usage_exit() {
  print_usage();
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version() {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void dump_options() {
  printf("\n");
  printf("%s\n", getVersion().c_str());
  printf("cwd %s\n", cwd);
  printf("cmdline %s\n", cmdline);
  printf("sysname  %s\n", uts.sysname);
  printf("hostname %s\n", uts.nodename);
  printf("machine  %s\n", uts.machine);
  printf("\n");

  // If both --lta and --template are specified, 
  // automatic detection of which direction the LTA goes by looking at which side matches the atlas template volume geomery. 
  // Automatically invert if necessary. --invert will not be applied in this case
  printf("atlas mesh collection      : %s\n", atlasMesh);
  printf("atlas template volume      : %s\n", (templateVol != NULL) ? templateVol : "null");
  printf("LTA transform matrix       : %s\n", (ltaFile != NULL) ? ltaFile : "null");
  printf("invert LTA                 : %s\n", (invertLTA) ? ((templateVol != NULL && ltaFile != NULL) ? "Ignored" : "Yes") : "No");
  printf("output surface             : %s\n", outSurf);
  printf("output priors as mri volume: %s\n", priorsmgz);
  printf("\n");
}
