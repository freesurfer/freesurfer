/**
 * @brief The internal 'program' API.
 *
 *  C_mpmProgs are overloaded classes that perform specific functions in the
 *  context of the Dijkstra system. The default contol system allows an
 *  external network-based program to drive the core dijkstra search. This
 *  process can be very slow on large-scale problems. To alleviate that,
 *  mpm_programs can be written. These map directly to external dsh scripts
 *  but without any of the network overhead.
 */
/*
 * Original Author: Rudolph Pienaar
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

#ifndef __C_MPM_PROG_H__
#define __C_MPM_PROG_H__

#include "error.h"
#include "fio.h"
#include "label.h"
#include "mri.h"
#include "mrisurf.h"

#include "C_mpmOverlay.h"
#include "env.h"
#include "general.h"

#include "oclDijkstraKernel.h"

#include <string>

const int MPMSTACKDEPTH = 64;

class C_mpmProg {

  //
  // Data members
  //

public:
  // type and var info
  std::string mstr_obj;                 // name of object class
  std::string mstr_name;                // name of object variable
  int         mid;                      // id of socket
  int         mverbosity;               // Debug related value
  int         mwarnings;                // Show warnings
  int         mstackDepth;              // Current procedure stackDepth
  std::string mstr_proc[MPMSTACKDEPTH]; // Used to track the current
                                        //+ procedure being executed
protected:
  // base class info

  s_env *mps_env; // Pointer to the main
                  //+ environment
  bool b_created;

  //
  // Method members
  //

public:
  //
  // Constructor / destructor block
  //
  void core_construct(std::string astr_name = "unnamed", int a_id = -1,
                      int a_verbosity = 0, int a_warnings = 0,
                      int a_stackDepth = 0, std::string astr_proc = "noproc");
  C_mpmProg(s_env *aps_env);
  virtual ~C_mpmProg();
  C_mpmProg(const C_mpmProg &C_mpmProg);
  C_mpmProg &operator=(const C_mpmProg &C_mpmProg);

  //
  // Error / warn /  print block
  //
  void debug_push(std::string astr_currentProc);
  void debug_pop();
  void error(std::string astr_msg = "", int code = 1);
  void warn(std::string astr_msg = "", int code = 1);
  void function_trace(std::string astr_msg       = "",
                      std::string astr_separator = "");
  void print();

  //
  // Access "housekeeping" state info
  //

  const std::string str_obj_get() const { return mstr_obj; };
  const std::string str_name_get() const { return mstr_name; };
  int               id_get() const { return mid; };
  int               verbosity_get() const { return mverbosity; };
  int               warnings_get() const { return mwarnings; };
  int               stackDepth_get() const { return mstackDepth; };

  void str_obj_set(std::string astr_val) { mstr_obj = astr_val; };
  void str_name_set(std::string astr_val) { mstr_name = astr_val; };
  void str_proc_set(int depth, std::string astr_proc) {
    mstr_proc[depth] = astr_proc;
  };
  void              id_set(int value) { mid = value; };
  void              verbosity_set(int value) { mverbosity = value; };
  void              warnings_set(int value) { mwarnings = value; };
  void              stackDepth_set(int value) { mstackDepth = value; };
  const std::string str_proc_get() const {
    return mstr_proc[stackDepth_get()];
  };
  const std::string str_proc_get(int i) { return mstr_proc[i]; };

  //
  // Core class methods
  //

  virtual int run() = 0;
};

class C_mpmProg_NOP : public C_mpmProg {

protected:
  int msleepSeconds;

public:
  C_mpmProg_NOP(s_env *aps_env);
  ~C_mpmProg_NOP();

  //
  // Access block
  //
  void sleepSeconds_set(int avalue) { msleepSeconds = avalue; };
  int  sleepSeconds_get() const { return (msleepSeconds); };

  //
  // Functional block
  //

  virtual int run();
};

///
/// \class C_mpmProg_pathFind
/// \brief This class simply finds a path between two vertices.
///
class C_mpmProg_pathFind : public C_mpmProg {

protected:
  int  mvertex_start;
  int  mvertex_end;
  int  mvertex_total;
  bool mb_surfaceRipClear;

public:
  C_mpmProg_pathFind(s_env *aps_env, int amvertex_start = 0,
                     int amvertex_end = -1);
  ~C_mpmProg_pathFind();

  //
  // Access block
  //
  void surfaceRipClear_set(bool avalue) { mb_surfaceRipClear = avalue; };
  int  surfaceRipClear_get() { return (mb_surfaceRipClear); };
  void vertexStart_set(int avalue) {
    mvertex_start        = avalue;
    mps_env->startVertex = avalue;
  };
  int  vertexStart_get() { return (mvertex_start); };
  void vertexEnd_set(int avalue) {
    mvertex_end        = avalue;
    mps_env->endVertex = avalue;
  };
  int  vertexEnd_get() { return (mvertex_end); };
  void print();

  //
  // Functional block
  //

  virtual int run();
  float       cost_compute(int start, int end);
};

///
/// \class C_mpmProg_autodijk
/// \brief This class implements the Dijkstra algorithm on the CPU.
///
class C_mpmProg_autodijk : public C_mpmProg {

protected:
  C_mpmOverlay *mpOverlayDistance; // A pointer to a 'distance'
                                   //+ mpmOverlay used to determine
                                   //+ anti-polar point.
  C_mpmOverlay *mpOverlayOrig;     // A pointer to the "original"
                                   //+ overlay created in the
                                   //+ base environment.

  int mvertex_polar; // For single sweep runs, this
                     //+ defines the "polar" vertex
                     //+ to start from
  int mvertex_start; // For multiple sweep runs, like
                     //+ worldMaps, this defines the
                     //+ start vertex
  int mvertex_step;  // the "step" vertex
  int mvertex_end;   // and the end vertex index
                     // Taken together, these
                     //+ allow for selective sweeps
                     //+ across the surface, and
                     //+ also for re-starting
                     //+ crashed processes.
  int     mvertex_total;
  int     m_costFunctionIndex;
  e_stats ms_stats; // A simple stats objects
  bool    mb_surfaceRipClear;
  bool    mb_worldMap;           // If true, generate a 'world
                                 //+ map' by looping exhaustively
                                 //+ over the entire surface and
                                 //+ only recording the cost
                                 //+ value at the vertex index
                                 //+ "farthest" from each start
                                 //+ vertex.
  bool mb_worldMapDistanceCalc;  // State calculation tracker
                                 //+ used for descriptive
                                 //+ output management
  bool mb_performExhaustive;     // If true, perform cost
                                 //+ calculations from polar
                                 //+ to every other vertex in
                                 //+ mesh -- useful only for
                                 //+ debugging/memory testing
                                 //+ since a search from
                                 //+ vertex->vertex will cover
                                 //+ whole mesh anyway in single
                                 //+ sweep.
  int mprogressIter;             // Number of iterations to
                                 //+ loop before showing
                                 //+ progress to stdout
  float *mpf_cost;               // Cost array as calculated by
                                 //+ by single call of autodijk.
                                 //+ In cases where multiple
                                 //+ calls are performed as part
                                 //+ of larger analysis, costs
                                 //+ might change. To keep
                                 //+ costs persistent, use the
                                 //+ mpf_persistent array.
  float *mpf_persistent;         // An array used to store values
                                 //+ that need to be persistent
                                 //+ across a whole autodijk run.
  float *mpf_fileSaveData;       // The save routine saves data
                                 //+ pointed to by this routine.
                                 //+ Typically pointer is managed
                                 //+ internally and should not
                                 //+ be manipulated outside this
                                 //+ class.
  bool mb_simpleStatsShow;       // If true, output some very
                                 //+ simple stats
  std::string mstr_costFileName; // Parsed from the environment
                                 //+ structure
  std::string mstr_costFullPath; // Full path to cost file

public:
  C_mpmProg_autodijk(s_env *aps_env);
  ~C_mpmProg_autodijk();

  //
  // Access block
  //
  void surfaceRipClear_set(bool avalue) { mb_surfaceRipClear = avalue; };
  int  surfaceRipClear_get() { return (mb_surfaceRipClear); };

  void costFile_set(std::string avalue) {
    mstr_costFileName = avalue;
    mstr_costFullPath = mps_env->str_workingDir + "/" + mstr_costFileName;
  };
  std::string costFile_get() { return (mstr_costFileName); };

  void vertexPolar_set(int avalue) { mvertex_polar = avalue; };
  int  vertexPolar_get() { return (mvertex_polar); };
  void vertexStart_set(int avalue) { mvertex_start = avalue; };
  int  vertexStart_get() { return (mvertex_start); };
  void vertexStep_set(int avalue) { mvertex_step = avalue; };
  int  vertexStep_get() { return (mvertex_step); };
  void vertexEnd_set(int avalue) { mvertex_end = avalue; };
  int  vertexEnd_get() { return (mvertex_end); };

  void worldMap_set(int avalue);
  bool worldMap_shouldCreate();

  void progressIter_set(int avalue) { mprogressIter = avalue; };
  int  progressIter_get() { return (mprogressIter); };
  void print();

  int vertexCosts_pack(e_stats &a_stats);
  //
  // Functional block
  //

  virtual int  run();
  float        cost_compute(int start, int end);
  e_FILEACCESS CURV_fileWrite();
};

///
/// \class C_mpmProg_autodijk_fast
/// \brief This class implements the Dijkstra algorithm using a fast
/// implementation
///        of the Single Source Shortest Path algorithm.  If compiled with
///        OpenCL, it will run a parallel accelerated version of the algorithm.
///        If not, it will just run a single threaded fast version of the
///        algorithm on the CPU
///
class C_mpmProg_autodijk_fast : public C_mpmProg_autodijk {
protected:
  ///
  /// Convert the freesurfer environment MRI surface representation into
  /// a representation that is formatted for consumption by the OpenCL
  /// algorithm
  /// \param graph Pointer to the graph data to generate
  ///
  void genOpenCLGraphRepresentation(GraphData *graph);

public:
  ///
  //  Constructor
  //
  C_mpmProg_autodijk_fast(s_env *aps_env);

  ///
  //  Destructor
  //
  virtual ~C_mpmProg_autodijk_fast();

  //
  // Functional block
  //
  virtual int run();
};

///
/// \class C_mpmProg_ROI
/// \brief This class paints an ROI about a vertex.
///
class C_mpmProg_ROI : public C_mpmProg {

protected:
  float            mf_radius;
  float            mf_plyIncrement;
  bool             mb_surfaceRipClear;
  std::vector<int> mv_vertex;
  std::string      mstr_vertexFile;
  std::string      mstr_labelFile;
  std::string      mstr_outputStem;
  bool             mb_ROIsInSeparateLabels;
  bool             mb_saveStaggered;
  bool             mb_boundaryOnly; // Toggle ROI only at border
  int              m_borderSize;

public:
  C_mpmProg_ROI(s_env *aps_env, std::string astr_vertexFile = "",
                float af_radius = 10);
  ~C_mpmProg_ROI();

  //
  // Access block
  //
  int borderSize() const { return m_borderSize; };

  bool boundaryOnly() const { return mb_boundaryOnly; };

  bool boundaryOnly(const bool &ab_val) {
    mb_boundaryOnly = ab_val;
    return mb_boundaryOnly;
  };

  void surfaceRipClear_set(bool avalue) { mb_surfaceRipClear = avalue; };
  int  surfaceRipClear_get() { return (mb_surfaceRipClear); };

  void plySaveStaggered_set(bool avalue) { mb_saveStaggered = avalue; };
  int  plySaveStaggered_get() { return (mb_saveStaggered); };

  void  radius_set(float af_value) { mf_radius = af_value; };
  float radius_get() { return mf_radius; };

  void  plyIncrement_set(float af_value) { mf_plyIncrement = af_value; };
  float plyIncrement_get() { return mf_plyIncrement; };

  std::vector<int> v_vertex_get() { return mv_vertex; };

  int border_mark();
  int vertexFile_load(std::string astr_fileName);
  int labelFile_load(std::string astr_fileName);
  int labelFile_save(std::string astr_fileName);
  int label_savePly(std::string astr_filePrefix, bool ab_staggered = false,
                    float af_plyIncrement = 1.0);

  void print();

  //
  // Functional block
  //

  virtual int run();
  float       cost_compute(int start, int end);
};

///
/// \class C_mpmProg_externalMesh
/// \brief This class uses an externalMesh with same vertices but different
/// edges
///
class C_mpmProg_externalMesh : public C_mpmProg {

protected:
  bool             mb_surfaceRipClear;
  std::vector<int> mv_vertex;
  std::string      mstr_meshFile;

public:
  C_mpmProg_externalMesh(s_env *aps_env, std::string astr_meshFile = "");
  ~C_mpmProg_externalMesh();

  //
  // Access block
  //
  void surfaceRipClear_set(bool avalue) { mb_surfaceRipClear = avalue; };
  int  surfaceRipClear_get() { return (mb_surfaceRipClear); };

  void print();

  //
  // Functional block
  //

  virtual int run();
  float       cost_compute(int start, int end);
};

#endif //__C_MPM_PROG_H__
