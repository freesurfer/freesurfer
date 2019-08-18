
#pragma once
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
#include "mrisurf_aaa.h"
namespace SurfaceFromMRIS_MP {
    typedef MRIS_MP Representation;


    namespace Existence {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( TopologyM::Face const & src                ) : Repr_Elt(src) {}
    Face::Face ( Topology::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( XYZPositionM::Face const & src             ) : Repr_Elt(src) {}
    Face::Face ( XYZPosition::Face const & src              ) : Repr_Elt(src) {}
    Face::Face ( XYZPositionConsequencesM::Face const & src ) : Repr_Elt(src) {}
    Face::Face ( XYZPositionConsequences::Face const & src  ) : Repr_Elt(src) {}
    Face::Face ( DistortM::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( Distort::Face const & src                  ) : Repr_Elt(src) {}
    Face::Face ( AnalysisM::Face const & src                ) : Repr_Elt(src) {}
    Face::Face ( Analysis::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                              ) {}
    Vertex::Vertex ( Representation* representation, size_t idx   ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                           ) : Repr_Elt(src) {}
    Vertex::Vertex ( TopologyM::Vertex const & src                ) : Repr_Elt(src) {}
    Vertex::Vertex ( Topology::Vertex const & src                 ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPositionM::Vertex const & src             ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPosition::Vertex const & src              ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPositionConsequencesM::Vertex const & src ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPositionConsequences::Vertex const & src  ) : Repr_Elt(src) {}
    Vertex::Vertex ( DistortM::Vertex const & src                 ) : Repr_Elt(src) {}
    Vertex::Vertex ( Distort::Vertex const & src                  ) : Repr_Elt(src) {}
    Vertex::Vertex ( AnalysisM::Vertex const & src                ) : Repr_Elt(src) {}
    Vertex::Vertex ( Analysis::Vertex const & src                 ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                     ) : Repr_Elt(src) {}

    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                               ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx    ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                           ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( TopologyM::MRIS_MP const & src                ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Topology::MRIS_MP const & src                 ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPositionM::MRIS_MP const & src             ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPosition::MRIS_MP const & src              ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPositionConsequencesM::MRIS_MP const & src ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPositionConsequences::MRIS_MP const & src  ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( DistortM::MRIS_MP const & src                 ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Distort::MRIS_MP const & src                  ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AnalysisM::MRIS_MP const & src                ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Analysis::MRIS_MP const & src                 ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                     ) : Repr_Elt(src) {}



    Surface::Surface (                                               ) {}
    Surface::Surface ( Representation* representation                ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src                           ) : Repr_Elt(src) {}
    Surface::Surface ( TopologyM::Surface const & src                ) : Repr_Elt(src) {}
    Surface::Surface ( Topology::Surface const & src                 ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPositionM::Surface const & src             ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPosition::Surface const & src              ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPositionConsequencesM::Surface const & src ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPositionConsequences::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( DistortM::Surface const & src                 ) : Repr_Elt(src) {}
    Surface::Surface ( Distort::Surface const & src                  ) : Repr_Elt(src) {}
    Surface::Surface ( AnalysisM::Surface const & src                ) : Repr_Elt(src) {}
    Surface::Surface ( Analysis::Surface const & src                 ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src                     ) : Repr_Elt(src) {}

    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }


    } // namespace Existence


    namespace Topology {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( XYZPositionM::Face const & src             ) : Repr_Elt(src) {}
    Face::Face ( XYZPosition::Face const & src              ) : Repr_Elt(src) {}
    Face::Face ( XYZPositionConsequencesM::Face const & src ) : Repr_Elt(src) {}
    Face::Face ( XYZPositionConsequences::Face const & src  ) : Repr_Elt(src) {}
    Face::Face ( DistortM::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( Distort::Face const & src                  ) : Repr_Elt(src) {}
    Face::Face ( AnalysisM::Face const & src                ) : Repr_Elt(src) {}
    Face::Face ( Analysis::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                              ) {}
    Vertex::Vertex ( Representation* representation, size_t idx   ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                           ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPositionM::Vertex const & src             ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPosition::Vertex const & src              ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPositionConsequencesM::Vertex const & src ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPositionConsequences::Vertex const & src  ) : Repr_Elt(src) {}
    Vertex::Vertex ( DistortM::Vertex const & src                 ) : Repr_Elt(src) {}
    Vertex::Vertex ( Distort::Vertex const & src                  ) : Repr_Elt(src) {}
    Vertex::Vertex ( AnalysisM::Vertex const & src                ) : Repr_Elt(src) {}
    Vertex::Vertex ( Analysis::Vertex const & src                 ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                     ) : Repr_Elt(src) {}

    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                               ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx    ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                           ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPositionM::MRIS_MP const & src             ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPosition::MRIS_MP const & src              ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPositionConsequencesM::MRIS_MP const & src ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPositionConsequences::MRIS_MP const & src  ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( DistortM::MRIS_MP const & src                 ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Distort::MRIS_MP const & src                  ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AnalysisM::MRIS_MP const & src                ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Analysis::MRIS_MP const & src                 ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                     ) : Repr_Elt(src) {}



    Surface::Surface (                                               ) {}
    Surface::Surface ( Representation* representation                ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src                           ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPositionM::Surface const & src             ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPosition::Surface const & src              ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPositionConsequencesM::Surface const & src ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPositionConsequences::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( DistortM::Surface const & src                 ) : Repr_Elt(src) {}
    Surface::Surface ( Distort::Surface const & src                  ) : Repr_Elt(src) {}
    Surface::Surface ( AnalysisM::Surface const & src                ) : Repr_Elt(src) {}
    Surface::Surface ( Analysis::Surface const & src                 ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src                     ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }


    } // namespace Topology


    namespace XYZPosition {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( XYZPositionConsequencesM::Face const & src ) : Repr_Elt(src) {}
    Face::Face ( XYZPositionConsequences::Face const & src  ) : Repr_Elt(src) {}
    Face::Face ( DistortM::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( Distort::Face const & src                  ) : Repr_Elt(src) {}
    Face::Face ( AnalysisM::Face const & src                ) : Repr_Elt(src) {}
    Face::Face ( Analysis::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                              ) {}
    Vertex::Vertex ( Representation* representation, size_t idx   ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                           ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPositionConsequencesM::Vertex const & src ) : Repr_Elt(src) {}
    Vertex::Vertex ( XYZPositionConsequences::Vertex const & src  ) : Repr_Elt(src) {}
    Vertex::Vertex ( DistortM::Vertex const & src                 ) : Repr_Elt(src) {}
    Vertex::Vertex ( Distort::Vertex const & src                  ) : Repr_Elt(src) {}
    Vertex::Vertex ( AnalysisM::Vertex const & src                ) : Repr_Elt(src) {}
    Vertex::Vertex ( Analysis::Vertex const & src                 ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                     ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(CANONICAL_VERTICES,c)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                               ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx    ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                           ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPositionConsequencesM::MRIS_MP const & src ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( XYZPositionConsequences::MRIS_MP const & src  ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( DistortM::MRIS_MP const & src                 ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Distort::MRIS_MP const & src                  ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AnalysisM::MRIS_MP const & src                ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Analysis::MRIS_MP const & src                 ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                     ) : Repr_Elt(src) {}



    Surface::Surface (                                               ) {}
    Surface::Surface ( Representation* representation                ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src                           ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPositionConsequencesM::Surface const & src ) : Repr_Elt(src) {}
    Surface::Surface ( XYZPositionConsequences::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( DistortM::Surface const & src                 ) : Repr_Elt(src) {}
    Surface::Surface ( Distort::Surface const & src                  ) : Repr_Elt(src) {}
    Surface::Surface ( AnalysisM::Surface const & src                ) : Repr_Elt(src) {}
    Surface::Surface ( Analysis::Surface const & src                 ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src                     ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }


    } // namespace XYZPosition


    namespace XYZPositionConsequences {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( DistortM::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( Distort::Face const & src                  ) : Repr_Elt(src) {}
    Face::Face ( AnalysisM::Face const & src                ) : Repr_Elt(src) {}
    Face::Face ( Analysis::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    float Face::area() const {
        return repr->f_area[idx];
    }
    angles_per_triangle_t Face::angle() const {
        return repr->f_angle[idx];
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->f_orig_angle[idx];
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    FloatXYZ Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( DistortM::Vertex const & src               ) : Repr_Elt(src) {}
    Vertex::Vertex ( Distort::Vertex const & src                ) : Repr_Elt(src) {}
    Vertex::Vertex ( AnalysisM::Vertex const & src              ) : Repr_Elt(src) {}
    Vertex::Vertex ( Analysis::Vertex const & src               ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(VERTEX_NORMALS,n)
        CASE(CANONICAL_VERTICES,c)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( DistortM::MRIS_MP const & src              ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Distort::MRIS_MP const & src               ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AnalysisM::MRIS_MP const & src             ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Analysis::MRIS_MP const & src              ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( DistortM::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( Distort::Surface const & src   ) : Repr_Elt(src) {}
    Surface::Surface ( AnalysisM::Surface const & src ) : Repr_Elt(src) {}
    Surface::Surface ( Analysis::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    float Surface::xctr() const {
        return repr->xctr;
    }
    float Surface::yctr() const {
        return repr->yctr;
    }
    float Surface::zctr() const {
        return repr->zctr;
    }
    float Surface::xlo() const {
        return repr->xlo;
    }
    float Surface::ylo() const {
        return repr->ylo;
    }
    float Surface::zlo() const {
        return repr->zlo;
    }
    float Surface::xhi() const {
        return repr->xhi;
    }
    float Surface::yhi() const {
        return repr->yhi;
    }
    float Surface::zhi() const {
        return repr->zhi;
    }
    float Surface::total_area() const {
        return repr->total_area;
    }
    double Surface::avg_vertex_area() const {
        return repr->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return repr->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return repr->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return repr->orig_area;
    }
    float Surface::neg_area() const {
        return repr->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return repr->neg_orig_area;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }
    
    void Surface::set_xctr(float to) {
        repr->xctr = to;
    }
    void Surface::set_yctr(float to) {
        repr->yctr = to;
    }
    void Surface::set_zctr(float to) {
        repr->zctr = to;
    }
    void Surface::set_xlo(float to) {
        repr->xlo = to;
    }
    void Surface::set_ylo(float to) {
        repr->ylo = to;
    }
    void Surface::set_zlo(float to) {
        repr->zlo = to;
    }
    void Surface::set_xhi(float to) {
        repr->xhi = to;
    }
    void Surface::set_yhi(float to) {
        repr->yhi = to;
    }
    void Surface::set_zhi(float to) {
        repr->zhi = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }


    } // namespace XYZPositionConsequences


    namespace Distort {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( AnalysisM::Face const & src                ) : Repr_Elt(src) {}
    Face::Face ( Analysis::Face const & src                 ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    float Face::area() const {
        return repr->f_area[idx];
    }
    angles_per_triangle_t Face::angle() const {
        return repr->f_angle[idx];
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->f_orig_angle[idx];
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    FloatXYZ Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AnalysisM::Vertex const & src              ) : Repr_Elt(src) {}
    Vertex::Vertex ( Analysis::Vertex const & src               ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::wnx() const {
        return repr->v_wnx[idx];
    }
    float Vertex::wny() const {
        return repr->v_wny[idx];
    }
    float Vertex::wnz() const {  //  white normal
        return repr->v_wnz[idx];
    }
    float Vertex::dx() const {
        return repr->v_dx[idx];
    }
    float Vertex::dy() const {
        return repr->v_dy[idx];
    }
    float Vertex::dz() const {  //  current change in position
        return repr->v_dz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    float Vertex::pialx() const {
        return repr->v_pialx[idx];
    }
    float Vertex::pialy() const {
        return repr->v_pialy[idx];
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->v_pialz[idx];
    }
    float Vertex::whitex() const {
        return repr->v_whitex[idx];
    }
    float Vertex::whitey() const {
        return repr->v_whitey[idx];
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->v_whitez[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->v_neg[idx];
    }
    char Vertex::border() const {  //  flag 
        return repr->v_border[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(VERTEX_NORMALS,n)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_wnx(float to) {
        repr->v_wnx[idx] = to;
    }
    void Vertex::set_wny(float to) {
        repr->v_wny[idx] = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->v_wnz[idx] = to;
    }
    void Vertex::set_dx(float to) {
        repr->v_dx[idx] = to;
    }
    void Vertex::set_dy(float to) {
        repr->v_dy[idx] = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->v_dz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_pialx(float to) {
        repr->v_pialx[idx] = to;
    }
    void Vertex::set_pialy(float to) {
        repr->v_pialy[idx] = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->v_pialz[idx] = to;
    }
    void Vertex::set_whitex(float to) {
        repr->v_whitex[idx] = to;
    }
    void Vertex::set_whitey(float to) {
        repr->v_whitey[idx] = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->v_whitez[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->v_neg[idx] = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->v_border[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AnalysisM::MRIS_MP const & src             ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( Analysis::MRIS_MP const & src              ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AnalysisM::Surface const & src ) : Repr_Elt(src) {}
    Surface::Surface ( Analysis::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    float Surface::xctr() const {
        return repr->xctr;
    }
    float Surface::yctr() const {
        return repr->yctr;
    }
    float Surface::zctr() const {
        return repr->zctr;
    }
    float Surface::xlo() const {
        return repr->xlo;
    }
    float Surface::ylo() const {
        return repr->ylo;
    }
    float Surface::zlo() const {
        return repr->zlo;
    }
    float Surface::xhi() const {
        return repr->xhi;
    }
    float Surface::yhi() const {
        return repr->yhi;
    }
    float Surface::zhi() const {
        return repr->zhi;
    }
    float Surface::total_area() const {
        return repr->total_area;
    }
    double Surface::avg_vertex_area() const {
        return repr->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return repr->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return repr->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return repr->orig_area;
    }
    float Surface::neg_area() const {
        return repr->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return repr->neg_orig_area;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }
    
    void Surface::set_xctr(float to) {
        repr->xctr = to;
    }
    void Surface::set_yctr(float to) {
        repr->yctr = to;
    }
    void Surface::set_zctr(float to) {
        repr->zctr = to;
    }
    void Surface::set_xlo(float to) {
        repr->xlo = to;
    }
    void Surface::set_ylo(float to) {
        repr->ylo = to;
    }
    void Surface::set_zlo(float to) {
        repr->zlo = to;
    }
    void Surface::set_xhi(float to) {
        repr->xhi = to;
    }
    void Surface::set_yhi(float to) {
        repr->yhi = to;
    }
    void Surface::set_zhi(float to) {
        repr->zhi = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }


    } // namespace Distort


    namespace Analysis {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    float Face::area() const {
        return repr->f_area[idx];
    }
    angles_per_triangle_t Face::angle() const {
        return repr->f_angle[idx];
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->f_orig_angle[idx];
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    FloatXYZ Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::wnx() const {
        return repr->v_wnx[idx];
    }
    float Vertex::wny() const {
        return repr->v_wny[idx];
    }
    float Vertex::wnz() const {  //  white normal
        return repr->v_wnz[idx];
    }
    float Vertex::dx() const {
        return repr->v_dx[idx];
    }
    float Vertex::dy() const {
        return repr->v_dy[idx];
    }
    float Vertex::dz() const {  //  current change in position
        return repr->v_dz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    float Vertex::pialx() const {
        return repr->v_pialx[idx];
    }
    float Vertex::pialy() const {
        return repr->v_pialy[idx];
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->v_pialz[idx];
    }
    float Vertex::whitex() const {
        return repr->v_whitex[idx];
    }
    float Vertex::whitey() const {
        return repr->v_whitey[idx];
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->v_whitez[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->v_neg[idx];
    }
    char Vertex::border() const {  //  flag 
        return repr->v_border[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(VERTEX_NORMALS,n)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_wnx(float to) {
        repr->v_wnx[idx] = to;
    }
    void Vertex::set_wny(float to) {
        repr->v_wny[idx] = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->v_wnz[idx] = to;
    }
    void Vertex::set_dx(float to) {
        repr->v_dx[idx] = to;
    }
    void Vertex::set_dy(float to) {
        repr->v_dy[idx] = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->v_dz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_pialx(float to) {
        repr->v_pialx[idx] = to;
    }
    void Vertex::set_pialy(float to) {
        repr->v_pialy[idx] = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->v_pialz[idx] = to;
    }
    void Vertex::set_whitex(float to) {
        repr->v_whitex[idx] = to;
    }
    void Vertex::set_whitey(float to) {
        repr->v_whitey[idx] = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->v_whitez[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->v_neg[idx] = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->v_border[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    float Surface::xctr() const {
        return repr->xctr;
    }
    float Surface::yctr() const {
        return repr->yctr;
    }
    float Surface::zctr() const {
        return repr->zctr;
    }
    float Surface::xlo() const {
        return repr->xlo;
    }
    float Surface::ylo() const {
        return repr->ylo;
    }
    float Surface::zlo() const {
        return repr->zlo;
    }
    float Surface::xhi() const {
        return repr->xhi;
    }
    float Surface::yhi() const {
        return repr->yhi;
    }
    float Surface::zhi() const {
        return repr->zhi;
    }
    float Surface::total_area() const {
        return repr->total_area;
    }
    double Surface::avg_vertex_area() const {
        return repr->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return repr->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return repr->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return repr->orig_area;
    }
    float Surface::neg_area() const {
        return repr->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return repr->neg_orig_area;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }
    
    void Surface::set_xctr(float to) {
        repr->xctr = to;
    }
    void Surface::set_yctr(float to) {
        repr->yctr = to;
    }
    void Surface::set_zctr(float to) {
        repr->zctr = to;
    }
    void Surface::set_xlo(float to) {
        repr->xlo = to;
    }
    void Surface::set_ylo(float to) {
        repr->ylo = to;
    }
    void Surface::set_zlo(float to) {
        repr->zlo = to;
    }
    void Surface::set_xhi(float to) {
        repr->xhi = to;
    }
    void Surface::set_yhi(float to) {
        repr->yhi = to;
    }
    void Surface::set_zhi(float to) {
        repr->zhi = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }


    } // namespace Analysis


    namespace ExistenceM {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }
    
    void Surface::set_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane)
        repr->underlyingMRIS->status = to;
    }
    void Surface::set_origxyz_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        repr->origxyz_status = to;
    }
    void Surface::set_patch(int to) {  //  if a patch of the surface
        repr->underlyingMRIS->patch = to;
    }


    } // namespace ExistenceM


    namespace TopologyM {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    
    void Face::set_v(size_t i, Vertex to) {
        cheapAssert(repr == to.repr); repr->f_v[idx][i] = to.idx;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }


    } // namespace TopologyM


    namespace XYZPositionM {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(CANONICAL_VERTICES,c)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_x(float to) {  //  current coordinates	
        repr->v_x[idx] = to;
    }
    void Vertex::set_y(float to) {  //  use MRISsetXYZ() to set
        repr->v_y[idx] = to;
    }
    void Vertex::set_z(float to) {
        repr->v_z[idx] = to;
    }
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }


    } // namespace XYZPositionM


    namespace XYZPositionConsequencesM {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    float Face::area() const {
        return repr->f_area[idx];
    }
    angles_per_triangle_t Face::angle() const {
        return repr->f_angle[idx];
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->f_orig_angle[idx];
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    FloatXYZ Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_area(float to) {
        repr->f_area[idx] = to;
    }
    void Face::set_angle(angles_per_triangle_t to) {
        repr->f_angle[idx] = to;
    }
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_norm(FloatXYZ to) {
        repr->f_norm[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(VERTEX_NORMALS,n)
        CASE(CANONICAL_VERTICES,c)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_area(float to) {
        repr->v_area[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    float Surface::xctr() const {
        return repr->xctr;
    }
    float Surface::yctr() const {
        return repr->yctr;
    }
    float Surface::zctr() const {
        return repr->zctr;
    }
    float Surface::xlo() const {
        return repr->xlo;
    }
    float Surface::ylo() const {
        return repr->ylo;
    }
    float Surface::zlo() const {
        return repr->zlo;
    }
    float Surface::xhi() const {
        return repr->xhi;
    }
    float Surface::yhi() const {
        return repr->yhi;
    }
    float Surface::zhi() const {
        return repr->zhi;
    }
    float Surface::total_area() const {
        return repr->total_area;
    }
    double Surface::avg_vertex_area() const {
        return repr->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return repr->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return repr->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return repr->orig_area;
    }
    float Surface::neg_area() const {
        return repr->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return repr->neg_orig_area;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }
    
    void Surface::set_xctr(float to) {
        repr->xctr = to;
    }
    void Surface::set_yctr(float to) {
        repr->yctr = to;
    }
    void Surface::set_zctr(float to) {
        repr->zctr = to;
    }
    void Surface::set_xlo(float to) {
        repr->xlo = to;
    }
    void Surface::set_ylo(float to) {
        repr->ylo = to;
    }
    void Surface::set_zlo(float to) {
        repr->zlo = to;
    }
    void Surface::set_xhi(float to) {
        repr->xhi = to;
    }
    void Surface::set_yhi(float to) {
        repr->yhi = to;
    }
    void Surface::set_zhi(float to) {
        repr->zhi = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }
    void Surface::set_avg_vertex_dist(double to) {  //  set by MRIScomputeAvgInterVertexDist
        repr->avg_vertex_dist = to;
    }
    void Surface::set_std_vertex_dist(double to) {
        repr->std_vertex_dist = to;
    }
    void Surface::set_orig_area(float to) {
        repr->orig_area = to;
    }
    void Surface::set_neg_area(float to) {
        repr->neg_area = to;
    }
    void Surface::set_neg_orig_area(float to) {  //  amount of original surface in folds
        repr->neg_orig_area = to;
    }


    } // namespace XYZPositionConsequencesM


    namespace DistortM {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    float Face::area() const {
        return repr->f_area[idx];
    }
    angles_per_triangle_t Face::angle() const {
        return repr->f_angle[idx];
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->f_orig_angle[idx];
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    FloatXYZ Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::wnx() const {
        return repr->v_wnx[idx];
    }
    float Vertex::wny() const {
        return repr->v_wny[idx];
    }
    float Vertex::wnz() const {  //  white normal
        return repr->v_wnz[idx];
    }
    float Vertex::dx() const {
        return repr->v_dx[idx];
    }
    float Vertex::dy() const {
        return repr->v_dy[idx];
    }
    float Vertex::dz() const {  //  current change in position
        return repr->v_dz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    float Vertex::pialx() const {
        return repr->v_pialx[idx];
    }
    float Vertex::pialy() const {
        return repr->v_pialy[idx];
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->v_pialz[idx];
    }
    float Vertex::whitex() const {
        return repr->v_whitex[idx];
    }
    float Vertex::whitey() const {
        return repr->v_whitey[idx];
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->v_whitez[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->v_neg[idx];
    }
    char Vertex::border() const {  //  flag 
        return repr->v_border[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(VERTEX_NORMALS,n)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_wnx(float to) {
        repr->v_wnx[idx] = to;
    }
    void Vertex::set_wny(float to) {
        repr->v_wny[idx] = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->v_wnz[idx] = to;
    }
    void Vertex::set_dx(float to) {
        repr->v_dx[idx] = to;
    }
    void Vertex::set_dy(float to) {
        repr->v_dy[idx] = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->v_dz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_pialx(float to) {
        repr->v_pialx[idx] = to;
    }
    void Vertex::set_pialy(float to) {
        repr->v_pialy[idx] = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->v_pialz[idx] = to;
    }
    void Vertex::set_whitex(float to) {
        repr->v_whitex[idx] = to;
    }
    void Vertex::set_whitey(float to) {
        repr->v_whitey[idx] = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->v_whitez[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->v_neg[idx] = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->v_border[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    float Surface::xctr() const {
        return repr->xctr;
    }
    float Surface::yctr() const {
        return repr->yctr;
    }
    float Surface::zctr() const {
        return repr->zctr;
    }
    float Surface::xlo() const {
        return repr->xlo;
    }
    float Surface::ylo() const {
        return repr->ylo;
    }
    float Surface::zlo() const {
        return repr->zlo;
    }
    float Surface::xhi() const {
        return repr->xhi;
    }
    float Surface::yhi() const {
        return repr->yhi;
    }
    float Surface::zhi() const {
        return repr->zhi;
    }
    float Surface::total_area() const {
        return repr->total_area;
    }
    double Surface::avg_vertex_area() const {
        return repr->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return repr->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return repr->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return repr->orig_area;
    }
    float Surface::neg_area() const {
        return repr->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return repr->neg_orig_area;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }
    
    void Surface::set_xctr(float to) {
        repr->xctr = to;
    }
    void Surface::set_yctr(float to) {
        repr->yctr = to;
    }
    void Surface::set_zctr(float to) {
        repr->zctr = to;
    }
    void Surface::set_xlo(float to) {
        repr->xlo = to;
    }
    void Surface::set_ylo(float to) {
        repr->ylo = to;
    }
    void Surface::set_zlo(float to) {
        repr->zlo = to;
    }
    void Surface::set_xhi(float to) {
        repr->xhi = to;
    }
    void Surface::set_yhi(float to) {
        repr->yhi = to;
    }
    void Surface::set_zhi(float to) {
        repr->zhi = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }


    } // namespace DistortM


    namespace AnalysisM {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}
    Face::Face ( AllM::Face const & src                     ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    float Face::area() const {
        return repr->f_area[idx];
    }
    angles_per_triangle_t Face::angle() const {
        return repr->f_angle[idx];
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->f_orig_angle[idx];
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    FloatXYZ Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::wnx() const {
        return repr->v_wnx[idx];
    }
    float Vertex::wny() const {
        return repr->v_wny[idx];
    }
    float Vertex::wnz() const {  //  white normal
        return repr->v_wnz[idx];
    }
    float Vertex::dx() const {
        return repr->v_dx[idx];
    }
    float Vertex::dy() const {
        return repr->v_dy[idx];
    }
    float Vertex::dz() const {  //  current change in position
        return repr->v_dz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    float Vertex::pialx() const {
        return repr->v_pialx[idx];
    }
    float Vertex::pialy() const {
        return repr->v_pialy[idx];
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->v_pialz[idx];
    }
    float Vertex::whitex() const {
        return repr->v_whitex[idx];
    }
    float Vertex::whitey() const {
        return repr->v_whitey[idx];
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->v_whitez[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->v_neg[idx];
    }
    char Vertex::border() const {  //  flag 
        return repr->v_border[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(VERTEX_NORMALS,n)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_wnx(float to) {
        repr->v_wnx[idx] = to;
    }
    void Vertex::set_wny(float to) {
        repr->v_wny[idx] = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->v_wnz[idx] = to;
    }
    void Vertex::set_dx(float to) {
        repr->v_dx[idx] = to;
    }
    void Vertex::set_dy(float to) {
        repr->v_dy[idx] = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->v_dz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_pialx(float to) {
        repr->v_pialx[idx] = to;
    }
    void Vertex::set_pialy(float to) {
        repr->v_pialy[idx] = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->v_pialz[idx] = to;
    }
    void Vertex::set_whitex(float to) {
        repr->v_whitex[idx] = to;
    }
    void Vertex::set_whitey(float to) {
        repr->v_whitey[idx] = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->v_whitez[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->v_neg[idx] = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->v_border[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}
    MRIS_MP::MRIS_MP ( AllM::MRIS_MP const & src                  ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    float Surface::xctr() const {
        return repr->xctr;
    }
    float Surface::yctr() const {
        return repr->yctr;
    }
    float Surface::zctr() const {
        return repr->zctr;
    }
    float Surface::xlo() const {
        return repr->xlo;
    }
    float Surface::ylo() const {
        return repr->ylo;
    }
    float Surface::zlo() const {
        return repr->zlo;
    }
    float Surface::xhi() const {
        return repr->xhi;
    }
    float Surface::yhi() const {
        return repr->yhi;
    }
    float Surface::zhi() const {
        return repr->zhi;
    }
    float Surface::total_area() const {
        return repr->total_area;
    }
    double Surface::avg_vertex_area() const {
        return repr->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return repr->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return repr->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return repr->orig_area;
    }
    float Surface::neg_area() const {
        return repr->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return repr->neg_orig_area;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }
    
    void Surface::set_xctr(float to) {
        repr->xctr = to;
    }
    void Surface::set_yctr(float to) {
        repr->yctr = to;
    }
    void Surface::set_zctr(float to) {
        repr->zctr = to;
    }
    void Surface::set_xlo(float to) {
        repr->xlo = to;
    }
    void Surface::set_ylo(float to) {
        repr->ylo = to;
    }
    void Surface::set_zlo(float to) {
        repr->zlo = to;
    }
    void Surface::set_xhi(float to) {
        repr->xhi = to;
    }
    void Surface::set_yhi(float to) {
        repr->yhi = to;
    }
    void Surface::set_zhi(float to) {
        repr->zhi = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }


    } // namespace AnalysisM


    namespace AllM {
    Face::Face (                                            ) {}
    Face::Face ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Face::Face ( Face const & src                           ) : Repr_Elt(src) {}

    Vertex Face::v(size_t i) const {
        return Vertex(repr,repr->f_v[idx][i]);
    }
    float Face::area() const {
        return repr->f_area[idx];
    }
    angles_per_triangle_t Face::angle() const {
        return repr->f_angle[idx];
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->f_orig_angle[idx];
    }
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    FloatXYZ Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_v(size_t i, Vertex to) {
        cheapAssert(repr == to.repr); repr->f_v[idx][i] = to.idx;
    }
    void Face::set_area(float to) {
        repr->f_area[idx] = to;
    }
    void Face::set_angle(angles_per_triangle_t to) {
        repr->f_angle[idx] = to;
    }
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_norm(FloatXYZ to) {
        repr->f_norm[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}

    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->v_x[idx];
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->v_y[idx];
    }
    float Vertex::z() const {
        return repr->v_z[idx];
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->v_origx[idx];
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->v_origy[idx];
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->v_origz[idx];
    }
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::wnx() const {
        return repr->v_wnx[idx];
    }
    float Vertex::wny() const {
        return repr->v_wny[idx];
    }
    float Vertex::wnz() const {  //  white normal
        return repr->v_wnz[idx];
    }
    float Vertex::dx() const {
        return repr->v_dx[idx];
    }
    float Vertex::dy() const {
        return repr->v_dy[idx];
    }
    float Vertex::dz() const {  //  current change in position
        return repr->v_dz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::cx() const {
        return repr->v_cx[idx];
    }
    float Vertex::cy() const {
        return repr->v_cy[idx];
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->v_cz[idx];
    }
    float Vertex::pialx() const {
        return repr->v_pialx[idx];
    }
    float Vertex::pialy() const {
        return repr->v_pialy[idx];
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->v_pialz[idx];
    }
    float Vertex::whitex() const {
        return repr->v_whitex[idx];
    }
    float Vertex::whitey() const {
        return repr->v_whitey[idx];
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->v_whitez[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->v_neg[idx];
    }
    char Vertex::border() const {  //  flag 
        return repr->v_border[idx];
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->v_ripflag[idx];
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
    
    #define CASE(WHICH, FIELD) \
      case WHICH: \
        *x = this->FIELD##x();  *y = this->FIELD##y();  *z = this->FIELD##z(); \
        break;
    
      switch (which) {
        CASE(CURRENT_VERTICES,)
        CASE(ORIGINAL_VERTICES,orig)
        CASE(VERTEX_NORMALS,n)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_x(float to) {  //  current coordinates	
        repr->v_x[idx] = to;
    }
    void Vertex::set_y(float to) {  //  use MRISsetXYZ() to set
        repr->v_y[idx] = to;
    }
    void Vertex::set_z(float to) {
        repr->v_z[idx] = to;
    }
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_wnx(float to) {
        repr->v_wnx[idx] = to;
    }
    void Vertex::set_wny(float to) {
        repr->v_wny[idx] = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->v_wnz[idx] = to;
    }
    void Vertex::set_dx(float to) {
        repr->v_dx[idx] = to;
    }
    void Vertex::set_dy(float to) {
        repr->v_dy[idx] = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->v_dz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_cx(float to) {
        repr->v_cx[idx] = to;
    }
    void Vertex::set_cy(float to) {
        repr->v_cy[idx] = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->v_cz[idx] = to;
    }
    void Vertex::set_pialx(float to) {
        repr->v_pialx[idx] = to;
    }
    void Vertex::set_pialy(float to) {
        repr->v_pialy[idx] = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->v_pialz[idx] = to;
    }
    void Vertex::set_whitex(float to) {
        repr->v_whitex[idx] = to;
    }
    void Vertex::set_whitey(float to) {
        repr->v_whitey[idx] = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->v_whitez[idx] = to;
    }
    void Vertex::set_area(float to) {
        repr->v_area[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->v_neg[idx] = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->v_border[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    MRIS_MP::MRIS_MP (                                            ) {}
    MRIS_MP::MRIS_MP ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    MRIS_MP::MRIS_MP ( MRIS_MP const & src                        ) : Repr_Elt(src) {}



    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    float Surface::xctr() const {
        return repr->xctr;
    }
    float Surface::yctr() const {
        return repr->yctr;
    }
    float Surface::zctr() const {
        return repr->zctr;
    }
    float Surface::xlo() const {
        return repr->xlo;
    }
    float Surface::ylo() const {
        return repr->ylo;
    }
    float Surface::zlo() const {
        return repr->zlo;
    }
    float Surface::xhi() const {
        return repr->xhi;
    }
    float Surface::yhi() const {
        return repr->yhi;
    }
    float Surface::zhi() const {
        return repr->zhi;
    }
    float Surface::total_area() const {
        return repr->total_area;
    }
    double Surface::avg_vertex_area() const {
        return repr->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return repr->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return repr->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return repr->orig_area;
    }
    float Surface::neg_area() const {
        return repr->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return repr->neg_orig_area;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->underlyingMRIS->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->underlyingMRIS->patch;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->underlyingMRIS->noscale;
    }
    
    void Surface::set_xctr(float to) {
        repr->xctr = to;
    }
    void Surface::set_yctr(float to) {
        repr->yctr = to;
    }
    void Surface::set_zctr(float to) {
        repr->zctr = to;
    }
    void Surface::set_xlo(float to) {
        repr->xlo = to;
    }
    void Surface::set_ylo(float to) {
        repr->ylo = to;
    }
    void Surface::set_zlo(float to) {
        repr->zlo = to;
    }
    void Surface::set_xhi(float to) {
        repr->xhi = to;
    }
    void Surface::set_yhi(float to) {
        repr->yhi = to;
    }
    void Surface::set_zhi(float to) {
        repr->zhi = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }
    void Surface::set_avg_vertex_dist(double to) {  //  set by MRIScomputeAvgInterVertexDist
        repr->avg_vertex_dist = to;
    }
    void Surface::set_std_vertex_dist(double to) {
        repr->std_vertex_dist = to;
    }
    void Surface::set_orig_area(float to) {
        repr->orig_area = to;
    }
    void Surface::set_neg_area(float to) {
        repr->neg_area = to;
    }
    void Surface::set_neg_orig_area(float to) {  //  amount of original surface in folds
        repr->neg_orig_area = to;
    }
    void Surface::set_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane)
        repr->underlyingMRIS->status = to;
    }
    void Surface::set_origxyz_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        repr->origxyz_status = to;
    }
    void Surface::set_patch(int to) {  //  if a patch of the surface
        repr->underlyingMRIS->patch = to;
    }


    } // namespace AllM
} // SurfaceFromMRIS_MP
