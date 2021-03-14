/**
 * @brief The internal 'program' API
 *
 *  The C_mpmOverlay class (and derivatives) presents an encapsulated
 *  interface to 'overlays' which are the 'curv' field component of
 *  of FreeSurfer meshes.
 *
 *  The class also maintains cost function evaluation and weight vector
 *  structures.
 *
 *  Since the class contains a pointed to the mris_pmake 'env' structure,
 *  its name is prefixed by 'mpm'.
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

#ifndef __C_MPMOVERLAY_H__
#define __C_MPMOVERLAY_H__

#include "error.h"
#include "fio.h"
#include "label.h"
#include "mri.h"
#include "mrisurf.h"

#include "env.h"
#include "general.h"
#include <string>
#include <vector>

const int MPMOVERLAYSTACKDEPTH = 64;
const int STRBUF               = 65536;

// enum typedef for the overlays relevant to FreeSurfer
typedef enum _e_overlay {
  e_K1,
  e_K2,
  e_H,
  e_K,
  e_S,
  e_BE,
  e_C,
  e_FI,
  e_thickness,
  e_curv,
  e_sulc,
  eoverlay
} EOVERLAY;

static std::string mstr_curvSuffix[] = {
    "K1.crv", "K2.crv", "K.crv",     "S.crv", "BE.crv",
    "C.crv",  "FI.crv", "thickness", "curv",  "sulc"};

static std::string mstr_curvPrefixTemplate[] = {
    "HEMI.SURF.", "HEMI.SURF.", "HEMI.SURF.", "HEMI.SURF.", "HEMI.SURF.",
    "HEMI.SURF.", "HEMI.SURF.", "HEMI.",      "HEMI.",      "HEMI."};

class C_mpmOverlay {

  //
  // Data members
  //

public:
  // type and var info
  std::string mstr_obj;    // name of object class
  std::string mstr_name;   // name of object variable
  int         mid;         // id of socket
  int         mverbosity;  // Debug related value
  int         mwarnings;   // Show warnings
  int         mstackDepth; // Current procedure stackDepth
  std::string mstr_proc[MPMOVERLAYSTACKDEPTH];
  // Used to track the current
  //+ procedure being executed
protected:
  // base class info

  s_env *mps_env; // Pointer to the main
                  //+ environment
  bool        mb_created;
  std::string mstr_surface;

  std::vector<float> mv_costWeight;       // vector of cost weights
  std::vector<float> mv_costWeightDel;    // delta weights for vector
  std::string        mstr_costWeightFile; // contains weights and del
                                          // weights

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
  C_mpmOverlay(s_env *aps_env);
  virtual ~C_mpmOverlay();
  C_mpmOverlay(const C_mpmOverlay &C_mpmOverlay);
  C_mpmOverlay &operator=(const C_mpmOverlay &C_mpmOverlay);

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

  const std::string strSurface_get() const { return mstr_surface; };

  bool          costVector_read(std::string astr_fileName = "");
  bool          costVector_write(std::string astr_fileName = "");
  virtual float costEdge_calc(int i, int j) = 0;
  std::string   curvFileName_get(EOVERLAY ae_overlay);
  e_FILEACCESS  CURV_fileRead(std::string astr_curvFileName, float *apf_curv[],
                              int *api_size);
};

//
// An "dummy" NOP subclass -- mostly for debugging.
//
class C_mpmOverlay_NOP : public C_mpmOverlay {

protected:
public:
  C_mpmOverlay_NOP(s_env *aps_env);
  ~C_mpmOverlay_NOP();

  //
  // Access block
  //

  //
  // Functional block
  //
  void          costWeightVector_init();
  virtual float costEdge_calc(int i, int j);
};

//
// An "unity" subclass -- returns a (weighted) '1' on each
// cost call.
//
class C_mpmOverlay_unity : public C_mpmOverlay {

protected:
public:
  C_mpmOverlay_unity(s_env *aps_env);
  ~C_mpmOverlay_unity();

  //
  // Access block
  //

  //
  // Functional block
  //
  void          costWeightVector_init();
  virtual float costEdge_calc(int i, int j);
};

//
// A "distance" subclass -- returns a (weighted) distance between
// two vertices by look up in the FreeSurfer mesh.
//
class C_mpmOverlay_distance : public C_mpmOverlay {

protected:
public:
  C_mpmOverlay_distance(s_env *aps_env);
  ~C_mpmOverlay_distance();

  //
  // Access block
  //

  //
  // Functional block
  //
  void          costWeightVector_init();
  virtual float costEdge_calc(int i, int j);
};

//
// A "distance" subclass -- returns a (weighted) distance between
// two vertices by look up in the FreeSurfer mesh.
//
class C_mpmOverlay_euclidean : public C_mpmOverlay {

protected:
public:
  C_mpmOverlay_euclidean(s_env *aps_env);
  ~C_mpmOverlay_euclidean();

  //
  // Access block
  //

  //
  // Functional block
  //
  void          costWeightVector_init();
  virtual float costEdge_calc(int i, int j);
};

//
// An overlay subclass that reads in the complete set of FreeSurfer
// curvature/thickness files
//
class C_mpmOverlay_FScurvs : public C_mpmOverlay {

  //
  // This class overlays several FreeSurfer "curvature" arrays
  //

protected:
  int mv_size;

public:
  C_mpmOverlay_FScurvs(s_env *aps_env, std::string astr_costWeightVector = "");
  ~C_mpmOverlay_FScurvs();

  //
  // Access block
  //
  float *overlay_getRef(EOVERLAY ae_overlay);

  //
  // Functional block
  //
  void          costWeightVector_init();
  virtual float costEdge_calc(int i, int j);
};

//
// An overlay subclass that reads in the complete set of FreeSurfer
// curvature/thickness files
//
class C_mpmOverlay_curvature : public C_mpmOverlay {

  //
  // This class reads in a single "curvature" type file and
  // uses these for edge weights.
  //

protected:
  int     mv_size;
  float * mpf_curvatureData;
  e_stats mes_curvStats;

public:
  C_mpmOverlay_curvature(s_env *aps_env, std::string astr_curvFileStem,
                         std::string astr_costWeightFile = "");
  ~C_mpmOverlay_curvature();

  //
  // Access block
  //

  //
  // Functional block
  //
  void          costWeightVector_init();
  virtual float costEdge_calc(int i, int j);
};

#endif //__C_MPMOVERLAY_H__
