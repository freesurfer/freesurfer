/**
 * @brief the 'main' for freeview
 *
 */
/*
 * Original Author: Ruopeng Wang
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
#include "MainApplication.h"
#include "MainWindow.h"
#include "MyCmdLineParser.h"
#include "MyUtils.h"
#include <QStringList>
#include <QDateTime>
#include <QMessageBox>
#include <QDir>
#include <QDebug>
#include <QStyleFactory>
#include "CursorFactory.h"
#include "vtkObject.h"
#include "LineProf.h"
#include <stdio.h>
#include <stdlib.h>
#include <fenv.h>
#include <QFile>
#include <QSurfaceFormat>

#if VTK_MAJOR_VERSION > 7
#include <QVTKOpenGLWidget.h>
#endif

#include "fsinit.h"
#include "log.h"
#include "chklc.h"
#include "utils.h"

const char* Progname;

#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
void myMessageOutput(QtMsgType type, const QMessageLogContext &context, const QString &msg)
{
  Q_UNUSED(context);

  if (msg.contains("sRGB profile") || msg.contains("QWidget::create")
      || msg.contains("unregister timer") || msg.contains("QPainter::") || msg.contains("Empty filename passed") ||
          msg.contains("posix collation implementation"))
    return;

  switch ((int)type)
  {
  case QtDebugMsg:
    fprintf(stdout, "%s\n", qPrintable(msg));
    fflush(0);
    break;
  case QtWarningMsg:
    fprintf(stderr, "%s\n", qPrintable(msg));
    fflush(0);
    break;
  case QtCriticalMsg:
    fprintf(stderr, "%s\n", qPrintable(msg));
    fflush(0);
    break;
  case QtFatalMsg:
    fprintf(stderr, "%s\n", qPrintable(msg));
    fflush(0);
    abort();
  }
}
#else
void myMessageOutput(QtMsgType type, const char *msg)
{
  QString qmsg = msg;
  if (qmsg.contains("sRGB profile") || qmsg.contains("QWidget::create") || qmsg.contains("unregister timer"))
    return;

  switch (type)
  {
  case QtDebugMsg:
    fprintf(stdout, "%s\n", msg);
    fflush(stdout);
    break;
  case QtWarningMsg:
    fprintf(stderr, "%s\n", msg);
    fflush(0);
    break;
  case QtCriticalMsg:
    fprintf(stderr, "%s\n", msg);
    fflush(0);
    break;
  case QtFatalMsg:
    fprintf(stderr, "%s\n", msg);
    fflush(0);
    abort();
  }
}
#endif

int main(int argc, char *argv[])
{
  Progname = argv[0];
  throwExceptions(true);

  putenv((char*)"SURFER_FRONTDOOR=");
  if (getenv("FS_DISABLE_LANG") == NULL)
    putenv((char*)"LANG=en_US");
#if (QT_VERSION >= QT_VERSION_CHECK(5, 0, 0))
  qInstallMessageHandler(myMessageOutput);
#else
  qInstallMsgHandler(myMessageOutput);
#endif

  LineProf::InitializePetsc(true);
  FSinit();
  setRandomSeed(-1L);

  CmdLineEntry cmdLineDesc[] =
  {
    CmdLineEntry( CMD_LINE_OPTION, "v", "volume", "<FILE>...", "Load one or multiple volume files. Available sub-options are: \n\n"
    "':colormap=name' Set colormap for display. Valid names are grayscale/lut/heat/jet/gecolor/nih/pet/binary. \n\n"
    "':grayscale=min,max' Set grayscale window values.\n\n"
    "':heatscale=min,mid,max' Set heat scale values. If only 2 values given, min and mid will be set to the same value.\n\n"
    "':heatscaleoptions=option1[,option2]' Set heat scale options. Options can be 'truncate','invert', or both.\n\n"
    "':colorscale=min,max' Set generic colorscale values for jet/gecolor/nih/pet.\n\n"
    "':lut=name' Set lookup table to the given name. Name can be the name of a stock color table or the filename of a color table file.\n\n"
    "':binary_color=color' Set the color of the non-zero voxels. Colormap must be set to binary. Color can be a color name such as 'red' or 3 values as RGB components of the color, e.g., '255,0,0'.\n\n"
    "':percentile=flag' Set min/mid/max thresholds as percentile.\n\n"
    "':vector=flag' Display 3 frame volume as vectors. flag can be 'yes', 'true' or '1'.\n\n"
    "':tensor=flag' Display 9 frame volume as tensors. flag can be 'yes', 'true' or '1'.\n\n"
    "':vector_width=width' Set line width of the vectors.\n\n"
    "':vector_norm_threshold=value' Set norm threshold for the vector display.\n\n"
    "':vector_skip=n' Skip n voxels in vector display. Default is 0 (no skipping).\n\n"
    "':render=flag' When displaying as vectors or tensors, render the glyph in the given form. For vector, flag can be 'line' as simple line or 'bar' as 3D bar (might be slow). For tensor, flag can be 'boxoid' or 'ellipsoid' (slow!).\n\n"
    "':inversion=flag' When displaying as vectors or tensors, invert the given component of the vectors. Valid flags are 'x', 'y' and 'z'.\n\n"
    "':outline=flag' Display labels as outline only. flag can be '1', 'yes' or 'true'.\n\n"
    "':reg=filename' or ':transform=filename' Set registration/transform file for the volume. filename can contain relative path to the volume file.\n\n"
    "':resample=method' or 'interpolation=method' Set the interpolation method when resampling is necessary. method can be 'nearest' (default), 'trilinear' or 'cubic.\n\n"
    "':opacity=value' Set the opacity of the volume layer. value ranges from 0 to 1.\n\n"
    "':mask=volume_name' Use the given volume to as mask for display. The mask volume must be loaded first.\n\n"
    "':isosurface=low_threshold,high_threshold' Set 3D display as isosurface. High_threshold is optional. If no threshold or simply 'on' is given, threshold will be either automatically determined or retrieved from the save previously settings.\n\n"
    "':isosurface_color=color' Set the color of the isosurface. Color can be a color name such as 'red' or 3 values as RGB components of the color, e.g., '255,0,0'.\n\n"
    "':isosurface_smooth=iterations' Set the number of smooth iterations. Default value is 5.\n\n"
    "':extract_all_regions=flag' Set isosurface to extract all regions. default setting is on.\n\n"
    "':isosurface_output=filename' Save isosurface to file. Extension can be .vtk or .stl.\n\n"
    "':surface_region=file' Load isosurface region(s) from the given file. isosurface display will automatically be turned on.\n\n"
    "':name=display_name' Set the display name of the volume.\n\n"
    "':lock=lock_status' Lock the volume layer so it will not be moved in the layer stack. Status can be '1' or 'true'.\n\n"
    "':visible=visibility' Set the initial visibility of the volume. Visibility can be '1' or '0' or 'true' or 'false'.\n\n"
    "':smoothed=flag' Set smoothed display for the volume. Flag can be '1' or '0' or 'true' or 'false'.\n\n"
    "':rgb=flag' Display 3-frame volume in RGB color. Voxel values must be in the range of 0~255. Flag can be '1' or '0' or 'true' or 'false'.\n\n"
    "':structure=name_or_value' Move the slice in the main viewport to where it has the most of the given structure.\n\n"
    "':ignore_header=flag' Ignore header information. Use the existing volume's header info. Flag can be '1' or '0' or 'true' or 'false'.\n\n"
    "':frame=number' Set active frame (0 based).\n\n"
    "Example:\nfreeview -v T1.mgz:colormap=heatscale:heatscale=10,100,200\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_SWITCH, "r", "resample", "", "Resample oblique data to standard RAS." ),
    CmdLineEntry( CMD_LINE_SWITCH, "conform", "conform", "", "Conform the volume to the first loaded volume." ),
    CmdLineEntry( CMD_LINE_SWITCH, "trilinear", "trilinear", "", "Use trilinear as the default resample method." ),
    CmdLineEntry( CMD_LINE_SWITCH, "cubic", "cubic", "", "Use cubic as the default resample method." ),
    CmdLineEntry( CMD_LINE_SWITCH, "smoothed", "smoothed", "", "Use smoothed display as the default display method for volumes." ),
    CmdLineEntry( CMD_LINE_OPTION, "colormap", "colormap", "<TYPE>", "Use the give colormap type as the colormap for all the volumes to be loaded after.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "dti", "dti", "<VECTOR> <FA>...", "Load one or more dti volumes. Need two files for each dti volume. First one is vector file. Second one is FA (brightness) file.", 2, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "tv", "tract-volume", "<FILE>...", "Load one or more tract volumes.Available sub-options are:\n\n"
    "':frame=frame_number' Show only a single tract by frame number.\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "f", "surface", "<FILE>...", "Load one or multiple surface files. Available sub-options are:\n\n"
    "':curvature=curvature_filename' Load curvature data from the given curvature file. By default .curv file will be loaded if available.\n\n"
    "':curvature_method=method' method to display curvature. available options are 'off', 'binary' or 'threshold'. default is 'threshold'.\n\n"
    "':curvature_setting=midpoint,slope' method to set curvature midpoint and slope.\n\n"
    "':overlay=overlay_filename' Load overlay data from file.\n\n"
    "':overlay_reg=overlay_registration_filename' Apply registration when loading overlay data.\n\n"
    "':overlay_method=method_name' Set overlay method. Valid names are 'linear', 'linearopaque' and 'piecewise'.\n\n"
    "':overlay_color=colorscale,settings' Set overlay color setttings. Valid names are 'colorwheel', 'truncate' and 'inverse'. Use comma to apply more than one.\n\n"
    "':overlay_opacity=opacity' Set opacity of overlay.\n\n"
    "':overlay_threshold=low,(mid,)high(,percentile)' Set overlay threshold values, separated by comma. When overlay method is linear or linearopaque, only 2 numbers (low and high) are needed. When method is piecewise, 3 numbers are needed. If last element is 'percentile', use the give numbers as percentile.\n\n"
    "':overlay_mask=filename(,invert)' Use given label file as mask for overlay. If invert is specified, use the inverted mask.\n\n"
    "':overlay_frame=frame_number' Set active frame of multi-frame overlay.\n\n"
    "':overlay_smooth=smooth_steps' Set smooth steps for overlay.\n\n"
    "':overlay_zorder=number' Set z-order for rendering overlay.\n\n"
    "':mrisp=filename' load parameterization data from file as overlay.\n\n"
    "':sphere=filename' load given sphere file.\n\n"
    "':patch=patch_filename' Load given patch file.\n\n"
    "':correlation=correlation_filename' Load correlation data from file. Correlation data is treated as a special kind of overlay data.\n\n"
    "':color=colorname' Set the base color of the surface. Color can be a color name such as 'red' or 3 values as RGB components of the color, e.g., '255,0,0'.\n\n"
    "':opacity=value' Set surface opacity. Value ranges from 0 to 1.\n\n"
    "':edgecolor=colorname' Set the color of the slice intersection outline on the surface. If set to 'overlay', will use overlay color\n\n"
    "':edgethickness=thickness' Set the thickness of the slice intersection outline on the surface. set 0 to hide it.\n\n"
    "':annot=filenames' Set annotation files to load.\n\n"
    "':annot_outline=flag' Show surface annotation as outline. flag can be 'true', 'yes' or '1'.\n\n"
    "':annot_zorder=number' Set z-order for rendering annotation.\n\n"
    "':name=display_name' Set the display name of the surface.\n\n':offset=x,y,z' Set the position offset of the surface. Useful for connectivity display.\n\n"
    "':visible=visibility' Set the initial visibility of the surface. Visibility can be '1' or '0' or 'true' or 'false'.\n\n':vector=filename' Load a vector file for display.\n\n"
    "':target_surf=filename' Load a target surface file for vectors to project on for 2D display.\n\n"
    "':label=filename' Load a surface label file.\n\n"
    "':label_outline=flag' Show surface labels as outline. flag can be 'true', 'yes' or '1'.\n\n"
    "':label_color=colorname' Set the color of the surface label.\n\n"
    "':label_centroid=flag' Move 3D cursor to the centroid of the label. flag can be 'true', 'yes' or '1'.\n\n"
    "':label_visible=flag' Set label visibility.\n\n"
    "':label_opacity=value' Set label opacity. Value ranges from 0 to 1.\n\n"
    "':label_zorder=number' Set z-order for rendering labels.\n\n"
    "':spline=filename' Load a spline file for display.\n\n"
    "':vertex=flag' Show surface vertices on both 2D and 3D views. flag can be 'true', 'on' or '1'.\n\n"
    "':vertexcolor=colorname' Set color of the vertices. Color can be a color name such as 'red' or 3 values as RGB components of the color, e.g., '255,0,0'.\n\n"
    "':goto=vertex_number' Move the cursor to the given vertex and center the view around it.\n\n"
    "':hide_in_3d=flag' Hide surface in 3D view. flag can be 'true', 'on' or '1'.\n\n"
    "':all=flag' Indicate to load all available surfaces. flag can be 'true', 'yes' or '1'.\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "edgecolor", "edgecolor", "<color>", "Use the give color as the edge color for all the surfaces to be loaded after.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "l", "label", "<FILE>...", "Load one or multiple label(ROI) files. Available sub-options are:\n\n"
    "':ref=ref_volume' Enter the name of the reference volume for this label file. The volume is one of the volumes given by -v option. \n\n"
    "':color=name' Set color of the label. Name can be a generic color name such as 'red' or 'lightgreen', or three integer values as RGB values ranging from 0 to 255. For example '255,0,0' is the same as 'red'.\n\n"
    "':opacity=value' Set the opacity of the label. value ranges from 0 to 1. \n\n"
    "':threshold=value' Set the threshold of the label.\n\n"
    "':centroid=flag' Move the cursor to the centroid of the label. flag can be '1', 'true' or 'yes'.\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "w", "way-points", "<FILE>...", "Load one or multiple way points files. Available sub-options are:\n\n"
    "':color=name' Set color of the way points. Name can be a generic color name such as 'red' or 'lightgreen', or three integer values as RGB values ranging from 0 to 255. For example '255,0,0' is the same as 'red'.\n\n"
    "':splinecolor=name' Set color of the spline.\n\n"
    "':splineheatmap=filename,min,mid,max' Set heatmap data file to load and heatmap settings for the spline.\n\n"
    "':radius=value' Set radius of the way points.\n\n"
    "':splineradius=value' Set radius of the spline tube.\n\n"
    "':name=display_name' Set the display name of the way points.\n\n"
    "':visible=visibility' Set the initial visibility of the way points. Visibility can be '1' or '0' or 'true' or 'false'.\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "c", "control-points", "<FILE>...", "Load one or multiple control points files. Available sub-options are:\n\n"
    "':color=name' Set color of the control points. Name can be a generic color name such as 'red' or 'lightgreen', or three integer values as RGB values ranging from 0 to 255. For example '255,0,0' is the same as 'red'.\n\n"
    "':radius=value' Set radius of the control points.\n\n"
    "':name=display_name' Set the display name of the control points.\n\n"
    "':visible=visibility' Set the initial visibility of the control points. Visibility can be '1' or '0' or 'true' or 'false'.\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "p-labels", "p-labels", "<FILES>...", "Load multiple p-label volume files.\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "p-prefix", "p-prefix", "<PREFIX>...", "Set the file name prefix for p-label volume. program will use this to figure out label name from file name.\n", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "p-lut", "p-lut", "<NAME>...", "Set the look up table name to use for p-label display. name can be the name of a stock lookup table or the file name of a lookup table file. default is the default freesurfer look up table.\n", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "cmat", "connectome-matrix", "<CMAT_FILE> <PARCELLATION_FILE>", "Load connectome matrix data files. Requires a cmat file and a parcellation volume file. Available sub-options are:\n\n'lut=color_table' Enter the name or file name of the color table to be used. Default is the FreeSurfer default LUT.\n", 2, 2 ),
    CmdLineEntry( CMD_LINE_OPTION, "fcd", "fcd", "<SUBJECT_DIR> <SUBJECT> [SUFFIX]", "Load FCD data. Requires subject directory and subject. Suffix is optional.\n", 2, 3 ),
    CmdLineEntry( CMD_LINE_OPTION, "t", "tract", "<FILE>...", "Load one or more tract files.\n\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "tc", "tract-cluster", "<DIRECTORY>", "Load tract cluster data from given directory.\n", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "recon", "recon", "<SUBJECT_NAME>...", "Load a series of pre-defined volumes and surfaces of given subject(s).\n", 1, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "lineprofile", "lineprofile", "<OUTPUT_FILE>", "Compute the thickness of layers along line profiles and export them to given csv file. Initial lines (waypoints) must be loaded in order with waypoint options. Available sub-options are:\n\n"
    "':spacing=value' Set spacing of the line profiles. Default value is 1.0.\n\n"
    "':resolution=value' Set resolution of the line profiles. Default value is 1.0.\n\n"
    "':offset=value' Set the offset to compute line profiles. Default value is 5.\n\n"
    "':segments=value' Set the number of segments on the line profiles. Default value is 100.\n", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "ss", "screenshot", "<FILE> <MAGIFICATION_FACTOR> <AUTO_TRIM>", "Take a screen shot of the main viewport and then quit the program. Default value for magnification factor is 1. AUTO_TRIM can be 'autotrim', 'true' or '1'. NOTE: AUTO_TRIM option is only available on Linux.", 1, 3 ),
    //    CmdLineEntry( CMD_LINE_OPTION, "fly", "fly-through", "<START_SLICE_NUMBER> <END_SLICE_NUMBER> <PREFIX>", "Fly through slices and take screenshot of each slice", 1, 3 ),
    CmdLineEntry( CMD_LINE_OPTION, "layout", "layout", "<STYLE>", "Set layout of the view panels as given. Accepted styles are 1, 2, 3 & 4. 1 is single panel. The rest are 3 different 4-panel styles.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "view", "view", "<VIEW>", "Set the 3D view as given. Accepted views are 'left', 'right', 'lateral', 'medial', 'anterior', 'posterior', 'inferior' and 'superior'. 'lateral' and 'medial' will only work when there is visible surface.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "viewport", "viewport", "<NAME>", "Set the main viewport as given. Accepted names are 'sagittal' or 'x', 'coronal' or 'y', 'axial' or 'z' and '3d'.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "viewsize", "viewsize", "<width> <height>", "Set the size of the main viewport. The size of the whole window will be changed accordingly.", 2, 2 ),
    CmdLineEntry( CMD_LINE_OPTION, "zoom", "zoom", "<FACTOR>", "Set zoom factor of the main viewport.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "cam", "camera", "<OPERATION1> <FACTOR1> <OPERATION2> <FACTOR2>...", "Set a series of camera operations for the 3D view. Valid operations are:\n\n"
    "'Azimuth' Rotate the camera about the view up vector centered at the focal point. The result is a horizontal rotation of the camera.\n\n"
    "'Dolly' Divide the camera's distance from the focal point by the given dolly value. Use a value greater than one to dolly-in toward the focal point, and use a value less than one to dolly-out away from the focal point.\n\n"
    "'Elevation' Rotate the camera about the cross product of the negative of the direction of projection and the view up vector, using the focal point as the center of rotation. The result is a vertical rotation of the scene.\n\n"
    "'Roll' Rotate the camera about the direction of projection. This will spin the camera about its axis.\n\n"
    "'Zoom' Same as 'Dolly'.\n\n"
    "Note that the order matters!\n\n"
    "For example: '-cam dolly 1.5 azimuth 30' will zoom in the camera by 1.5 times and then rotate it along the view up vector by 30 degrees.\n", 2, 1000 ),
    CmdLineEntry( CMD_LINE_OPTION, "ras", "ras", "<X> <Y> <Z> [Option]", "Set cursor location at the given RAS coordinate. If option is given as 'tkreg', the entered coordinate is in tkreg coordinate.", 3, 4 ),
    CmdLineEntry( CMD_LINE_OPTION, "slice", "slice", "<X> <Y> <Z>", "Set cursor location at the given slice numbers of the first loaded volume.", 3, 3 ),
    CmdLineEntry( CMD_LINE_SWITCH, "cc", "center-at-cursor", "", "Center the view at cursor location." ),
    CmdLineEntry( CMD_LINE_SWITCH, "timecourse", "timecourse", "", "Show time course window." ),
    CmdLineEntry( CMD_LINE_SWITCH, "transform-volume", "transform-volume", "", "Show transform volume window." ),
    CmdLineEntry( CMD_LINE_SWITCH, "colorscale", "colorscale", "", "Show color scale bar on the main view." ),
    CmdLineEntry( CMD_LINE_OPTION, "write-slice-intersection", "write-slice-intersection", "<plane> <filename_format> <start> <end>", "Save polylines of the surface intersected by slice. plane can be sag/cor/hor.", 4, 4 ),
    CmdLineEntry( CMD_LINE_OPTION, "cmd", "command", "<FILE>", "Load freeview commands from a text file.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "hide", "hide", "<LAYER_TYPE>", "Hide the current layer. This is useful for loading comands by -cmd option. Valid LAYER_TYPEs are volume, surface, label, etc.", 1, 1 ),
    CmdLineEntry( CMD_LINE_OPTION, "unload", "unload", "<LAYER_TYPE>", "Unload/Close the current layer. Useful for loading comands by -cmd option. Valid LAYER_TYPEs are volume, surface, label, etc.", 1, 1 ),
    CmdLineEntry( CMD_LINE_SWITCH, "rorder", "reverse-order", "", "Load layers in reversed order." ),
    CmdLineEntry( CMD_LINE_SWITCH, "nocursor", "nocursor", "", "Hide the cursor." ),
    CmdLineEntry( CMD_LINE_SWITCH, "hide-3d-slices", "hide-3d-slices", "", "Hide slices in 3D view." ),
    CmdLineEntry( CMD_LINE_SWITCH, "hide-3d-frames", "hide-3d-frames", "", "Hide slice frames in 3D view." ),
    CmdLineEntry( CMD_LINE_SWITCH, "hide-x-slice", "hide-x-slice", "", "Hide x slice in 3D view." ),
    CmdLineEntry( CMD_LINE_SWITCH, "hide-y-slice", "hide-y-slice", "", "Hide y slice in 3D view." ),
    CmdLineEntry( CMD_LINE_SWITCH, "hide-z-slice", "hide-z-slice", "", "Hide z slice in 3D view." ),
    CmdLineEntry( CMD_LINE_OPTION, "subtitle", "subtitle", "<TEXT>", "Add subtitle to freeview window caption.", 1, 1 ),
    CmdLineEntry( CMD_LINE_SWITCH, "auto-load-surf", "auto-load-surf", "", "Do not automatically load sphere or other supplemental surface data." ),
    CmdLineEntry( CMD_LINE_SWITCH, "quit", "quit", "", "Quit freeview. Useful for scripting or loading comands by -cmd option." ),
    CmdLineEntry( CMD_LINE_SWITCH, "noquit", "noquit", "", "Do not quit freeview after screenshot command." ),
    CmdLineEntry( CMD_LINE_SWITCH, "stdin", "stdin", "", "Listening stdin for freeview command sent by other programs." ),
    CmdLineEntry( CMD_LINE_SWITCH, "verbose", "verbose", "", "Print out a lot more information, such as vertex coordinate of each click, etc." ),
    CmdLineEntry( CMD_LINE_SWITCH, "continue", "continue", "", "Continue as far as possible if an error occured" ),
    CmdLineEntry( CMD_LINE_NONE )
  };

  char progDesc[] = "Volume and surface viewer and editor for freesurfer.";

  MyCmdLineParser cmd( (const char*)"freeview", (CmdLineEntry*)cmdLineDesc );
  cmd.SetProgramDescription( progDesc );
  if ( !cmd.Parse( argc, argv ) )
  {
    return 1;
  }

#if VTK_MAJOR_VERSION > 7
  QSurfaceFormat::setDefaultFormat(QVTKOpenGLWidget::defaultFormat());
#endif

#ifdef Q_OS_LINUX
#if (QT_VERSION < QT_VERSION_CHECK(5, 0, 0))
  QApplication::setStyle( "Cleanlooks" );
#endif
#endif

  MainApplication app(argc, argv);
  app.setOrganizationName("Massachusetts General Hospital");
  app.setOrganizationDomain("nmr.mgh.harvard.edu");
#ifdef Q_OS_LINUX
  app.setApplicationName("freeview");
#else
  app.setApplicationName("FreeView");
#endif

#ifndef Q_OS_MAC
  QIcon icon(":/resource/icons/app_32.png");
  icon.addFile(":/resource/icons/app_64.png");
  app.setWindowIcon(icon);
#endif

  // global initialization
  CursorFactory::Initialize();
  qsrand(QDateTime::currentDateTime().toTime_t());
  vtkObject::GlobalWarningDisplayOff();
  // fesetround(FE_TONEAREST);
#ifdef CYGWIN
  QDir dir;
  dir.mkpath("/cygdrive/c/tmp");
#endif

  MainWindow w(NULL, &cmd);
  w.show();

  // check license
  char license_msg[2000];
  if (!chklc2(license_msg))
  {
    QMessageBox::warning(&w, "License Error", license_msg);
    w.close();
    return 1;
  }

  if (!w.ParseCommand(argc, argv, true))
  {
    w.close();
    return 1;
  }

  int ret = app.exec();

  LineProf::FinalizePetsc();
  if (w.HadError())
    ret = 1;

  return ret;
}
