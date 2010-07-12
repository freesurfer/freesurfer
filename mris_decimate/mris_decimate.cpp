/*
 * Original Author: Dan Ginsburg (@ Children's Hospital Boston)
 * CVS Revision Info:
 *    $Author: ginsburg $
 *    $Date: 2010/07/12 19:38:12 $
 *    $Revision: 1.1 $
 *
 * Copyright (C) 2010,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

///
/// \file mris_decimate.cpp
///
/// \brief Brief description
/// Reduce the number of vertices and faces in a surface. 
///
/// \b NAME
///
/// mris_decimate 
///
/// \b SYNPOSIS
///
/// mris_decimate [options] <input surface> <output surface>
///
/// \b DESCRIPTION
///
/// This tool reduces the number of triangles in a surface using the
/// the vtkDecimatePro class documented at:
///
///            http://www.vtk.org/doc/release/5.2/html/a00310.htmls 
///
/// Please see the VTK documentation for details on the decimation algorithm.
/// mris_decimate will read in an existing surface and write out a new one
/// with less triangles.  The decimation level and other options can be provided
/// on the command-line.
///
/// By default, this tool is command-line driven.  If you wish to use a VTK-based
/// interactive visualization of the surface, uncomment the '#define BUILT_IN_VISUALIZER'
/// below.
///


// $Id: mris_decimate.cpp,v 1.1 2010/07/12 19:38:12 ginsburg Exp $


#ifdef BUILT_IN_VISUALIZER
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkCommand.h>
#include <vtkWidgetEvent.h>
#include <vtkCallbackCommand.h>
#include <vtkWidgetEventTranslator.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkSliderWidget.h>
#include <vtkSliderRepresentation2D.h>
#include <vtkProperty2D.h>
#include <vtkTextProperty.h>
#include <vtkTextActor.h>
#include <vtkWindowToImageFilter.h>
#include <vtkPNGWriter.h>
#endif // BUILT_IN_VISUALIZER

#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkIdList.h>
#include <vtkCellArray.h>
#include <vtkPolyData.h>
#include <vtkDecimatePro.h>
#include <vtkSmartPointer.h>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include <sstream>

#include "mris_decimate.h"

extern "C"
{
#include "macros.h"
#include "utils.h"
#include "fio.h"
#include "version.h"
#include "cmdargs.h"
#include "error.h"
#include "diag.h"
}



///
//  Function Prototypes
//
static int  get_option(int argc, char **argv);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int decimateSurface(MRI_SURFACE *mris, const DECIMATION_OPTIONS &decimationOptions, char *outFilePath);
std::string getStatsAsString(vtkPolyData *decimatedMesh, vtkDecimatePro *decimate);
int main(int argc, char *argv[]) ;

///
//  Global Variables
//
static char vcid[] = "$Id: mris_decimate.cpp,v 1.1 2010/07/12 19:38:12 ginsburg Exp $";
char *Progname = NULL;
char *cmdline;
int debug=0;
int checkoptsonly=0;
bool gSavePNGImage = false;
char *gSavePNGImagePath = NULL;
struct utsname uts;
DECIMATION_OPTIONS gDecimationOptions;

///////////////////////////////////////////////////////////////////////////////////
//  
//  Public Functions
//
//

#ifdef BUILT_IN_VISUALIZER

///
/// \class vtkSliderCallback
/// \brief Callback to handle adjustment to 3D slider in visualizer mode
///
class vtkSliderCallback : public vtkCommand
{
public:
    static vtkSliderCallback *New() 
    {
        return new vtkSliderCallback;
    }

    virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
        vtkSliderWidget *sliderWidget = reinterpret_cast<vtkSliderWidget*>(caller);
        this->Decimate->SetTargetReduction(static_cast<vtkSliderRepresentation *>(sliderWidget->GetRepresentation())->GetValue());
        Decimate->Update();

        vtkSmartPointer<vtkPolyData> decimatedMesh = Decimate->GetOutput();
 
        std::cout << "After decimation" << std::endl << "------------" << std::endl;
 
        std::cout << "There are " << decimatedMesh->GetNumberOfPolys() << " triangles." << std::endl;
        std::cout << "There are " << decimatedMesh->GetNumberOfPoints() << " points." << std::endl;

        TextActor->SetInput(getStatsAsString(decimatedMesh, Decimate).c_str());
    }
    vtkSliderCallback():Decimate(0) {}
    vtkDecimatePro *Decimate;
    vtkTextActor *TextActor;
};
#endif // BUILT_IN_VISUALIZER

///
/// \fn int decimateSurface(MRI_SURFACE *mris, const DECIMATION_OPTIONS &decimationOptions, char *outFilePath)
/// \brief This function performs decimation on the input surface and outputs the new surface to a file.
/// \param mris Input loaded MRI_SURFACE to decimate
/// \param decimationOptions Options controlling the decimation arguments (see DECIMATION_OPTIONS)
/// \param outFilePath Full path to output file to write decimated surface to
/// \return 0 on success, 1 on failure
///
int decimateSurface(MRI_SURFACE *mris, const DECIMATION_OPTIONS &decimationOptions, char *outFilePath)
{    

    // ---------------------------------------------------------------------------
    // Copy the vertex and face buffer to VTK data structures
    // ---------------------------------------------------------------------------
    vtkSmartPointer<vtkPolyData> mesh = vtkSmartPointer<vtkPolyData>::New();
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
    vtkSmartPointer<vtkCellArray> polys = vtkSmartPointer<vtkCellArray>::New();

    for( int v = 0; v < mris->nvertices; v++)
    {
        points->InsertPoint(v, mris->vertices[v].x,
                               mris->vertices[v].y,
                               mris->vertices[v].z);
    }

    for (int f = 0; f < mris->nfaces; f++)
    {
        vtkIdType pts[3];
        pts[0] =  mris->faces[f].v[0];
        pts[1] =  mris->faces[f].v[1];
        pts[2] =  mris->faces[f].v[2];

        polys->InsertNextCell(3, pts);
    }

    mesh->SetPoints(points);
    mesh->SetPolys(polys);

    std::cout << std::endl << "Before decimation" << std::endl << "------------" << std::endl;
 
    std::cout << "There are " << mesh->GetNumberOfPolys() << " triangles." << std::endl;
    std::cout << "There are " << mesh->GetNumberOfPoints() << " points." << std::endl;


    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
    decimate->SetInput(mesh);
    
    // ---------------------------------------------------------------------------
    // Set Decimation options from the command-line
    // ---------------------------------------------------------------------------
    decimate->SetTargetReduction(decimationOptions.reductionLevel);
    if(decimationOptions.setPreserveTopology)
        decimate->SetPreserveTopology(decimationOptions.preserveTopology);
    if(decimationOptions.setSplitting)
        decimate->SetSplitting(decimationOptions.splitting);
    if(decimationOptions.setFeatureAngle)
        decimate->SetFeatureAngle(decimationOptions.featureAngle);
    if(decimationOptions.setBoundaryVertexDeletion)
        decimate->SetBoundaryVertexDeletion(decimationOptions.boundaryVertexDeletion);
    if(decimationOptions.setDegree)
        decimate->SetDegree(decimationOptions.degree);    

    // ---------------------------------------------------------------------------
    // Perform the decimation
    // ---------------------------------------------------------------------------
    decimate->Update();

    vtkSmartPointer<vtkPolyData> decimatedMesh = decimate->GetOutput();
 
    std::cout << "After decimation" << std::endl << "------------" << std::endl;
 
    std::cout << "There are " << decimatedMesh->GetNumberOfPolys() << " triangles." << std::endl;
    std::cout << "There are " << decimatedMesh->GetNumberOfPoints() << " points." << std::endl;

    
#ifdef BUILT_IN_VISUALIZER

    // ---------------------------------------------------------------------------
    // If compiled with BUILT_IN_VISUALIZER, display a simple GUI using VTK
    // that can adjust the decimation level interactively.
    // ---------------------------------------------------------------------------
    
    // Now we'll look at it.
    vtkSmartPointer<vtkPolyDataMapper> polyMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
      polyMapper->SetInput(decimatedMesh);
      polyMapper->SetScalarRange(0,7);
    vtkSmartPointer<vtkActor> polyActor = vtkSmartPointer<vtkActor>::New();
      polyActor->SetMapper(polyMapper);

    // The usual rendering stuff.
    vtkSmartPointer<vtkCamera> camera = vtkSmartPointer<vtkCamera>::New();
      camera->SetPosition(-1,0,0);
      camera->SetFocalPoint(0,0,0);
      camera->Roll(90.0);

    vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
    vtkSmartPointer<vtkRenderWindow> renWin = vtkSmartPointer<vtkRenderWindow>::New();
    renWin->AddRenderer(renderer);
    renWin->SetWindowName("mris_decimate_gui");

    vtkSmartPointer<vtkRenderWindowInteractor> iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    iren->SetRenderWindow(renWin);

    renderer->AddActor(polyActor);
      renderer->SetActiveCamera(camera);
      renderer->ResetCamera();
      renderer->SetBackground(0,0,0);

    renWin->SetSize(1024,768);

    vtkSmartPointer<vtkSliderRepresentation2D> sliderRep =
        vtkSmartPointer<vtkSliderRepresentation2D>::New();
    sliderRep->SetMinimumValue(0.0);
    sliderRep->SetMaximumValue(1.0);
    sliderRep->SetValue(decimationOptions.reductionLevel);
    sliderRep->SetTitleText("Decimation Level");

    //set color properties:
    //change the color of the knob that slides
    sliderRep->GetSliderProperty()->SetColor(1,0,0);//red

    //change the color of the text indicating what the slider controls
    sliderRep->GetTitleProperty()->SetColor(1,0,0);//red

    //change the color of the text displaying the value
    sliderRep->GetLabelProperty()->SetColor(1,0,0);//red

    //change the color of the knob when the mouse is held on it
    sliderRep->GetSelectedProperty()->SetColor(0,1,0);//green

    //change the color of the bar
    sliderRep->GetTubeProperty()->SetColor(1,1,0);//yellow

    //change the color of the ends of the bar
    sliderRep->GetCapProperty()->SetColor(1,1,0);//yellow

    sliderRep->GetPoint1Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint1Coordinate()->SetValue(.25 ,0.1);
    sliderRep->GetPoint2Coordinate()->SetCoordinateSystemToNormalizedDisplay();
    sliderRep->GetPoint2Coordinate()->SetValue(.75, 0.1);
 
    vtkSmartPointer<vtkSliderWidget> sliderWidget =
        vtkSmartPointer<vtkSliderWidget>::New();
    sliderWidget->SetInteractor(iren);
    sliderWidget->SetRepresentation(sliderRep);
    sliderWidget->SetAnimationModeToJump();
    sliderWidget->EnabledOn();

    // Text widget
    vtkSmartPointer<vtkTextActor> textActor = 
        vtkSmartPointer<vtkTextActor>::New();    

    textActor->SetInput(getStatsAsString(decimatedMesh, decimate).c_str());
    textActor->GetTextProperty()->SetColor(1.0, 0.0, 0.0);
    textActor->GetTextProperty()->SetJustificationToLeft();
    textActor->GetTextProperty()->SetVerticalJustificationToBottom();
    textActor->GetPositionCoordinate()->SetCoordinateSystemToNormalizedViewport();
    textActor->GetPosition2Coordinate()->SetCoordinateSystemToNormalizedViewport();
    textActor->GetPositionCoordinate()->SetValue(0.0, 0.01);
    textActor->GetPosition2Coordinate()->SetValue(0.5, 0.35);
    renderer->AddActor2D(textActor);
    
    vtkSmartPointer<vtkSliderCallback> callback =
        vtkSmartPointer<vtkSliderCallback>::New();
    callback->Decimate = decimate;
    callback->TextActor = textActor;
 
    sliderWidget->AddObserver(vtkCommand::InteractionEvent,callback);

    iren->Initialize();

    // interact with data
    renWin->Render();

    if (gSavePNGImage)
    {
        textActor->VisibilityOff();
        sliderRep->VisibilityOff();
        vtkSmartPointer<vtkWindowToImageFilter> w2if = vtkSmartPointer<vtkWindowToImageFilter>::New();
        vtkSmartPointer<vtkPNGWriter> pngWriter = vtkSmartPointer<vtkPNGWriter>::New();

        w2if->SetInput(renWin);
        pngWriter->SetInput(w2if->GetOutput());
        pngWriter->SetFileName(gSavePNGImagePath);
        pngWriter->Write();
        return 0;
    }
    iren->Start();
#endif // BUILT_IN_VISUALIZER

    // ---------------------------------------------------------------------------
    // Get the decimated mesh and modify the 'mris' to use this mesh rather
    // than the original
    // ---------------------------------------------------------------------------
    decimatedMesh = decimate->GetOutput();
    
    // Free the vertex and face buffers
    for (int vno = 0 ; vno < mris->nvertices ; vno++)
    {
        if (mris->vertices[vno].f)
        {
            free(mris->vertices[vno].f) ;
            mris->vertices[vno].f = NULL ;
        }
        if (mris->vertices[vno].n)
        {
            free(mris->vertices[vno].n) ;
            mris->vertices[vno].n = NULL ;
        }
        if (mris->vertices[vno].dist)
        {
            free(mris->vertices[vno].dist) ;
            mris->vertices[vno].dist = NULL ;
        }
        if (mris->vertices[vno].dist_orig)
        {
            free(mris->vertices[vno].dist_orig) ;
            mris->vertices[vno].dist_orig = NULL ;
        }
        if (mris->vertices[vno].v)
        {
            free(mris->vertices[vno].v) ;
            mris->vertices[vno].v = NULL ;
        }
    }
    free(mris->vertices);
    free(mris->faces);

    mris->nvertices = decimatedMesh->GetNumberOfPoints();
    mris->nfaces = decimatedMesh->GetNumberOfPolys();

    // Allocate at the new size       
    mris->vertices = (VERTEX *)calloc(mris->nvertices, sizeof(VERTEX)) ;
    if (!mris->vertices)
    {
        ErrorExit(ERROR_NO_MEMORY,
                  "decimateSurface(%d, %d): could not allocate vertices",
                  mris->nvertices, sizeof(VERTEX));
    }
    mris->faces = (FACE *)calloc(mris->nfaces, sizeof(FACE)) ;
    if (!mris->faces)
    {
        ErrorExit(ERROR_NO_MEMORY,
                  "decimateSurface(%d, %d): could not allocate faces",
                  mris->nfaces, sizeof(FACE));
    }

    // Copy the data into the vertex and face buffers
    for(int v = 0; v < mris->nvertices; v++)
    {
        double *vert = decimatedMesh->GetPoint(v);
        mris->vertices[v].x = (float)vert[0];
        mris->vertices[v].y = (float)vert[1];
        mris->vertices[v].z = (float)vert[2];
    }
    vtkIdType *decimatedFaces = decimatedMesh->GetPolys()->GetPointer();
    for (int f = 0; f < mris->nfaces; f++)        
    {
        mris->faces[f].v[0] = decimatedFaces[f * 4 + 1];
        mris->faces[f].v[1] = decimatedFaces[f * 4 + 2];
        mris->faces[f].v[2] = decimatedFaces[f * 4 + 3];
    }

    // Finally, write out the final mesh and free the MRI_SURFACE
    MRISwrite(mris, outFilePath);
    MRISfree(&mris);

    return 0;
}


///////////////////////////////////////////////////////////////////////////////////
//  
//  Private Functions
//
//

///
/// Get the stats from a mesh and decimation as a string for output
/// \param decimatedMesh vtkPolyData for decimated mesh
/// \param decimate vtkDecimatePro that is used for decimation
///
std::string getStatsAsString(vtkPolyData *decimatedMesh, 
                             vtkDecimatePro *decimate)
{
    std::stringstream textStr;
    textStr << "Triangles: " << decimatedMesh->GetNumberOfPolys() << std::endl;
    textStr << "Vertices: " << decimatedMesh->GetNumberOfPoints() << std::endl;
    textStr << "Reduction Level: " << decimate->GetTargetReduction() << std:: endl;
    textStr << "Preserve Topology: " << decimate->GetPreserveTopology() << std::endl;
    textStr << "Splitting: " << decimate->GetSplitting() << std::endl;
    textStr << "Feature Angle: " << decimate->GetFeatureAngle() << std::endl;
    textStr << "Boundary Vertex Deletion: " << decimate->GetBoundaryVertexDeletion() << std::endl;
    textStr << "Degree: " << decimate->GetDegree();

    return textStr.str();
}

///
/// \fn Main entrypoint for mris_decimate
/// \return 0 on succesful run, 1 on error
///
int main(int argc, char *argv[]) 
{
    // Initialize Decimation options
    memset(&gDecimationOptions, 0, sizeof(DECIMATION_OPTIONS));
    gDecimationOptions.reductionLevel = 0.5; // Default decimation level if not specified

    char **av;
    char *in_fname, out_fpath[STRLEN] ;
    int ac, nargs;
    MRI_SURFACE *mris ;

    nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
    if (nargs && argc - nargs == 1) 
        exit (0);
    Progname = argv[0] ;
    argc -= nargs;
    cmdline = argv2cmdline(argc,argv);
    uname(&uts);

    if (argc < 3)
        usage_exit();

    ErrorInit(NULL, NULL, NULL) ;
    DiagInit(NULL, NULL, NULL) ;
    
    ac = argc ;
    av = argv ;
    for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) 
    {
        nargs = get_option(argc, argv) ;
        argc -= nargs ;
        argv += nargs ;
    }

    if (argc < 3)
        usage_exit();

    in_fname = argv[1] ;
    FileNameAbsolute(argv[2], out_fpath);

    dump_options(stdout);

    mris = MRISfastRead(in_fname) ;
    if (!mris)
    {
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                  Progname, in_fname) ;
    }

    decimateSurface(mris, gDecimationOptions, out_fpath);

    return 0;
}

///
/// \fn int get_option(int argc, char **argv)
/// \brief Parses a command-line argument
/// \param argc - number of command line arguments
/// \param argv - pointer to a character pointer
///
static int get_option(int argc, char *argv[]) 
{
    int  nargs = 0 ;
    char *option ;

    option = argv[1] + 1 ;            /* past '-' */
    if (!stricmp(option, "-help"))
    {
        print_help() ;
    }
    else if (!stricmp(option, "-version"))
    {
        print_version() ;
    } 
    else switch (toupper(*option)) 
    {
    case 'R':
        gDecimationOptions.reductionLevel = atof(argv[2]) ;
        printf("using reduction = %2.2f\n", gDecimationOptions.reductionLevel) ;
        nargs = 1 ;
        break ;
    case 'T':
        gDecimationOptions.setPreserveTopology = true;
        gDecimationOptions.preserveTopology = atoi(argv[2]);
        printf("using preserveTopology = %d\n", gDecimationOptions.preserveTopology) ;
        nargs = 1;
        break;
    case 'S':
        gDecimationOptions.setSplitting = true;
        gDecimationOptions.splitting = atoi(argv[2]);
        printf("using splitting = %d\n", gDecimationOptions.splitting) ;
        nargs = 1;
        break;
    case 'F':
        gDecimationOptions.setFeatureAngle = true;
        gDecimationOptions.featureAngle = atof(argv[2]);
        printf("using featureAngle = %2.2f\n", gDecimationOptions.featureAngle) ;
        nargs = 1;
        break;
    case 'B':
        gDecimationOptions.setBoundaryVertexDeletion = true;
        gDecimationOptions.boundaryVertexDeletion = atoi(argv[2]);
        printf("using boundaryVertexDeletion = %d\n", gDecimationOptions.boundaryVertexDeletion) ;
        nargs = 1;
        break;
    case 'D':
        gDecimationOptions.setDegree = true;
        gDecimationOptions.degree = atoi(argv[2]);
        printf("using degree = %d\n", gDecimationOptions.degree) ;
        nargs = 1;
        break;
    case 'P':
        gSavePNGImage = true;
        gSavePNGImagePath = argv[2];
        printf("save screenshot to '%s'\n", gSavePNGImagePath);
        nargs = 1;
        break;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


///
/// \fn static void usage_exit(void)
/// \brief Prints usage and exits
///
static void usage_exit(void) 
{
    print_usage() ;
    exit(1) ;
}

///
/// \fn static void print_usage(void)
/// \brief Prints usage and returns (does not exit)
///
static void print_usage(void) 
{
    // Create an object to get default values
    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();

    printf("USAGE: %s [options] <input surface> <output surface>\n",Progname) ;
    printf("\n");
    printf("This program reduces the number of triangles in a surface and\n");
    printf("outputs the new surface to a file.  All of the command line\n");
    printf("options are used to set parameters of the decimation algorithm\n");
    printf("which is provided by the VTK class 'vtkDecimatePro' documented at:\n\n");
    printf("    http://www.vtk.org/doc/release/5.2/html/a00310.htmls\n\n");
    printf("For more information on what the options do, please consult the VTK\n");
    printf("documentation for 'vtkDecimatePro'\n");
    printf("\nValid options are:\n\n");
    printf("   -r reductionLevel\n");
    printf("         percentage to reduce triangles in surface\n");
    printf("         by (value between 0<-->1.0, default: 0.5) \n");
    printf("   -t preserveTopology\n");
    printf("         whether or not to preserve topology (0 or 1, default: %d)\n", decimate->GetPreserveTopology());
    printf("   -s splitting\n");
    printf("         whether to allow splitting of the mesh (0 or 1, default: %d)\n", decimate->GetSplitting());
    printf("   -f featureAngle\n");
    printf("         specify the angle used to define what an edge is (default: %2.2f)\n", decimate->GetFeatureAngle());
    printf("   -b boundaryVertexDeletion\n");
    printf("         whether to allow deletion of vertices on the boundary \n");
    printf("         of a mesh (default: %d)\n", decimate->GetBoundaryVertexDeletion());
    printf("   -d degree\n");
    printf("         if the degree of vertex is greater than degree then\n");
    printf("         it will be split (default: %d)\n", decimate->GetDegree());
#ifdef BUILT_IN_VISUALIZER
    printf("   -p png_file\n");
    printf("         save a screenshot to <png_file>\n");
#endif
    printf("\n");
    printf("\n");
    printf("   --help      print out information on how to use this program\n");
    printf("   --version   print out version and exit\n");
    printf("\n");
    printf("\n");
}

///
/// \fn static void print_help(void)
/// \brief Prints help and exits
///
static void print_help(void) 
{
    print_usage() ;
    exit(1) ;
}

///
/// \fn static void print_version(void)
/// \brief Prints version and exits
///
static void print_version(void) 
{
    printf("%s\n", vcid) ;
    exit(1) ;
}


///
/// \fn static void dump_options(FILE *fp)
/// \brief Prints command-line options to the given file pointer
/// \param FILE *fp - file pointer
///
static void dump_options(FILE *fp) 
{
    // Create an object to get default values
    vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
    
    fprintf(fp,"\n");
    fprintf(fp,"%s\n",vcid);
    fprintf(fp,"cmdline %s\n",cmdline);
    fprintf(fp,"sysname  %s\n",uts.sysname);
    fprintf(fp,"hostname %s\n",uts.nodename);
    fprintf(fp,"machine  %s\n",uts.machine);
    fprintf(fp,"user     %s\n",VERuser());

    fprintf(fp,"\nreductionLevel            %f\n", gDecimationOptions.reductionLevel);
    if (gDecimationOptions.setPreserveTopology)
        fprintf(fp,"preserveTopology          %d\n", gDecimationOptions.preserveTopology) ;
    else
        fprintf(fp,"preserveTopology          %d\n", decimate->GetPreserveTopology()) ;
    

    if (gDecimationOptions.setSplitting)
        fprintf(fp,"splitting                 %d\n", gDecimationOptions.splitting) ;
    else
        fprintf(fp,"splitting                 %d\n", decimate->GetSplitting()) ;


    if (gDecimationOptions.setFeatureAngle)
        fprintf(fp,"featureAngle              %2.2f\n", gDecimationOptions.featureAngle) ;
    else
        fprintf(fp,"featureAngle              %2.2f\n", decimate->GetFeatureAngle()) ;

    if (gDecimationOptions.setBoundaryVertexDeletion)
        fprintf(fp,"boundaryVertexDeletion    %d\n", gDecimationOptions.boundaryVertexDeletion) ;
    else
        fprintf(fp,"boundaryVertexDeletion    %d\n", decimate->GetBoundaryVertexDeletion()) ;


    if (gDecimationOptions.setDegree)
        fprintf(fp,"degree                    %d\n", gDecimationOptions.degree) ;
    else
        fprintf(fp,"degree                    %d\n", decimate->GetDegree()) ;

    if (gSavePNGImage)
    {
        fprintf(fp, "Saving screenshot to '%s'\n", gSavePNGImagePath);
    }

    return;
}
