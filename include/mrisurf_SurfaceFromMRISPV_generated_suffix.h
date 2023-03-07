
#pragma once
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
#include "mrisurf_aaa.h"
namespace SurfaceFromMRISPV {
    typedef MRISPV Representation;


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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
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

    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
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

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
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
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


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

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
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

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
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

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->f_gradNorm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( DistortM::Vertex const & src               ) : Repr_Elt(src) {}
    Vertex::Vertex ( Distort::Vertex const & src                ) : Repr_Elt(src) {}
    Vertex::Vertex ( AnalysisM::Vertex const & src              ) : Repr_Elt(src) {}
    Vertex::Vertex ( Analysis::Vertex const & src               ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( DistortM::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( Distort::Surface const & src   ) : Repr_Elt(src) {}
    Surface::Surface ( AnalysisM::Surface const & src ) : Repr_Elt(src) {}
    Surface::Surface ( Analysis::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return repr->strips[i];
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
    float Surface::x0() const {  //  center of spherical expansion
        return repr->x0;
    }
    float Surface::y0() const {
        return repr->y0;
    }
    float Surface::z0() const {
        return repr->z0;
    }
    float Surface::max_curv() const {
        return repr->max_curv;
    }
    float Surface::min_curv() const {
        return repr->min_curv;
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
    int Surface::zeros() const {
        return repr->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return repr->hemisphere;
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return repr->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return repr->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return repr->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return repr->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return repr->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::nlabels() const {
        return repr->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return repr->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        repr->strips[i] = to;
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
    void Surface::set_x0(float to) {  //  center of spherical expansion
        repr->x0 = to;
    }
    void Surface::set_y0(float to) {
        repr->y0 = to;
    }
    void Surface::set_z0(float to) {
        repr->z0 = to;
    }
    void Surface::set_max_curv(float to) {
        repr->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        repr->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        repr->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        repr->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        repr->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        repr->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        repr->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        repr->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        repr->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        repr->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        repr->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->f_gradNorm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AnalysisM::Vertex const & src              ) : Repr_Elt(src) {}
    Vertex::Vertex ( Analysis::Vertex const & src               ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    float Vertex::pnx() const {
        return repr->v_pnx[idx];
    }
    float Vertex::pny() const {
        return repr->v_pny[idx];
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->v_pnz[idx];
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
    float Vertex::onx() const {
        return repr->v_onx[idx];
    }
    float Vertex::ony() const {
        return repr->v_ony[idx];
    }
    float Vertex::onz() const {  //  original normal
        return repr->v_onz[idx];
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
    float Vertex::odx() const {
        return repr->v_odx[idx];
    }
    float Vertex::ody() const {
        return repr->v_ody[idx];
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->v_odz[idx];
    }
    float Vertex::tdx() const {
        return repr->v_tdx[idx];
    }
    float Vertex::tdy() const {
        return repr->v_tdy[idx];
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->v_tdz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::curvbak() const {
        return repr->v_curvbak[idx];
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->v_val[idx];
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->v_imag_val[idx];
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
    float Vertex::tx() const {
        return repr->v_tx[idx];
    }
    float Vertex::ty() const {
        return repr->v_ty[idx];
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->v_tz[idx];
    }
    float Vertex::t2x() const {
        return repr->v_t2x[idx];
    }
    float Vertex::t2y() const {
        return repr->v_t2y[idx];
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->v_t2z[idx];
    }
    float Vertex::targx() const {
        return repr->v_targx[idx];
    }
    float Vertex::targy() const {
        return repr->v_targy[idx];
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->v_targz[idx];
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
    float Vertex::l4x() const {
        return repr->v_l4x[idx];
    }
    float Vertex::l4y() const {
        return repr->v_l4y[idx];
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->v_l4z[idx];
    }
    float Vertex::infx() const {
        return repr->v_infx[idx];
    }
    float Vertex::infy() const {
        return repr->v_infy[idx];
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->v_infz[idx];
    }
    float Vertex::fx() const {
        return repr->v_fx[idx];
    }
    float Vertex::fy() const {
        return repr->v_fy[idx];
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->v_fz[idx];
    }
    int Vertex::px() const {
        return repr->v_px[idx];
    }
    int Vertex::qx() const {
        return repr->v_qx[idx];
    }
    int Vertex::py() const {
        return repr->v_py[idx];
    }
    int Vertex::qy() const {
        return repr->v_qy[idx];
    }
    int Vertex::pz() const {
        return repr->v_pz[idx];
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->v_qz[idx];
    }
    float Vertex::e1x() const {
        return repr->v_e1x[idx];
    }
    float Vertex::e1y() const {
        return repr->v_e1y[idx];
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_e1z[idx];
    }
    float Vertex::e2x() const {
        return repr->v_e2x[idx];
    }
    float Vertex::e2y() const {
        return repr->v_e2y[idx];
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_e2z[idx];
    }
    float Vertex::pe1x() const {
        return repr->v_pe1x[idx];
    }
    float Vertex::pe1y() const {
        return repr->v_pe1y[idx];
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_pe1z[idx];
    }
    float Vertex::pe2x() const {
        return repr->v_pe2x[idx];
    }
    float Vertex::pe2y() const {
        return repr->v_pe2y[idx];
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_pe2z[idx];
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->v_nc[idx];
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->v_val2[idx];
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->v_valbak[idx];
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->v_val2bak[idx];
    }
    float Vertex::stat() const {  //  statistic 
        return repr->v_stat[idx];
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->v_undefval[idx];
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->v_old_undefval[idx];
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->v_fixedval[idx];
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->v_fieldsign[idx];
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->v_fsmask[idx];
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->v_d[idx];
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->v_annotation[idx];
    }
    char Vertex::oripflag() const {
        return repr->v_oripflag[idx];
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->v_origripflag[idx];
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->v_vp[idx];
    }
    float Vertex::theta() const {
        return repr->v_theta[idx];
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->v_phi[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    float Vertex::group_avg_area() const {
        return repr->v_group_avg_area[idx];
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->v_K[idx];
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->v_H[idx];
    }
    float Vertex::k1() const {
        return repr->v_k1[idx];
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->v_k2[idx];
    }
    float Vertex::mean() const {
        return repr->v_mean[idx];
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->v_mean_imag[idx];
    }
    float Vertex::std_error() const {
        return repr->v_std_error[idx];
    }
    uint Vertex::flags() const {
        return repr->v_flags[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    int Vertex::cropped() const {
        return repr->v_cropped[idx];
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->v_marked[idx];
    }
    short Vertex::marked2() const {
        return repr->v_marked2[idx];
    }
    short Vertex::marked3() const {
        return repr->v_marked3[idx];
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
        CASE(PIAL_NORMALS,pn)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(TMP_VERTICES,t)
        CASE(TMP2_VERTICES,t2)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        CASE(INFLATED_VERTICES,inf)
        CASE(FLATTENED_VERTICES,f)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
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
    void Vertex::set_pnx(float to) {
        repr->v_pnx[idx] = to;
    }
    void Vertex::set_pny(float to) {
        repr->v_pny[idx] = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->v_pnz[idx] = to;
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
    void Vertex::set_onx(float to) {
        repr->v_onx[idx] = to;
    }
    void Vertex::set_ony(float to) {
        repr->v_ony[idx] = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->v_onz[idx] = to;
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
    void Vertex::set_odx(float to) {
        repr->v_odx[idx] = to;
    }
    void Vertex::set_ody(float to) {
        repr->v_ody[idx] = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->v_odz[idx] = to;
    }
    void Vertex::set_tdx(float to) {
        repr->v_tdx[idx] = to;
    }
    void Vertex::set_tdy(float to) {
        repr->v_tdy[idx] = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->v_tdz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->v_curvbak[idx] = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->v_val[idx] = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->v_imag_val[idx] = to;
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
    void Vertex::set_tx(float to) {
        repr->v_tx[idx] = to;
    }
    void Vertex::set_ty(float to) {
        repr->v_ty[idx] = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->v_tz[idx] = to;
    }
    void Vertex::set_t2x(float to) {
        repr->v_t2x[idx] = to;
    }
    void Vertex::set_t2y(float to) {
        repr->v_t2y[idx] = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->v_t2z[idx] = to;
    }
    void Vertex::set_targx(float to) {
        repr->v_targx[idx] = to;
    }
    void Vertex::set_targy(float to) {
        repr->v_targy[idx] = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->v_targz[idx] = to;
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
    void Vertex::set_l4x(float to) {
        repr->v_l4x[idx] = to;
    }
    void Vertex::set_l4y(float to) {
        repr->v_l4y[idx] = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->v_l4z[idx] = to;
    }
    void Vertex::set_infx(float to) {
        repr->v_infx[idx] = to;
    }
    void Vertex::set_infy(float to) {
        repr->v_infy[idx] = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->v_infz[idx] = to;
    }
    void Vertex::set_fx(float to) {
        repr->v_fx[idx] = to;
    }
    void Vertex::set_fy(float to) {
        repr->v_fy[idx] = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->v_fz[idx] = to;
    }
    void Vertex::set_px(int to) {
        repr->v_px[idx] = to;
    }
    void Vertex::set_qx(int to) {
        repr->v_qx[idx] = to;
    }
    void Vertex::set_py(int to) {
        repr->v_py[idx] = to;
    }
    void Vertex::set_qy(int to) {
        repr->v_qy[idx] = to;
    }
    void Vertex::set_pz(int to) {
        repr->v_pz[idx] = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->v_qz[idx] = to;
    }
    void Vertex::set_e1x(float to) {
        repr->v_e1x[idx] = to;
    }
    void Vertex::set_e1y(float to) {
        repr->v_e1y[idx] = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_e1z[idx] = to;
    }
    void Vertex::set_e2x(float to) {
        repr->v_e2x[idx] = to;
    }
    void Vertex::set_e2y(float to) {
        repr->v_e2y[idx] = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_e2z[idx] = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->v_pe1x[idx] = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->v_pe1y[idx] = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_pe1z[idx] = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->v_pe2x[idx] = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->v_pe2y[idx] = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_pe2z[idx] = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->v_nc[idx] = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->v_val2[idx] = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->v_valbak[idx] = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->v_val2bak[idx] = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->v_stat[idx] = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->v_undefval[idx] = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->v_old_undefval[idx] = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->v_fixedval[idx] = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->v_fieldsign[idx] = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->v_fsmask[idx] = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->v_d[idx] = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->v_annotation[idx] = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->v_oripflag[idx] = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->v_origripflag[idx] = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->v_vp[idx] = to;
    }
    void Vertex::set_theta(float to) {
        repr->v_theta[idx] = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->v_phi[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->v_group_avg_area[idx] = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->v_K[idx] = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->v_H[idx] = to;
    }
    void Vertex::set_k1(float to) {
        repr->v_k1[idx] = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->v_k2[idx] = to;
    }
    void Vertex::set_mean(float to) {
        repr->v_mean[idx] = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->v_mean_imag[idx] = to;
    }
    void Vertex::set_std_error(float to) {
        repr->v_std_error[idx] = to;
    }
    void Vertex::set_flags(uint to) {
        repr->v_flags[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_cropped(int to) {
        repr->v_cropped[idx] = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->v_marked[idx] = to;
    }
    void Vertex::set_marked2(short to) {
        repr->v_marked2[idx] = to;
    }
    void Vertex::set_marked3(short to) {
        repr->v_marked3[idx] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AnalysisM::Surface const & src ) : Repr_Elt(src) {}
    Surface::Surface ( Analysis::Surface const & src  ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return repr->strips[i];
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
    float Surface::x0() const {  //  center of spherical expansion
        return repr->x0;
    }
    float Surface::y0() const {
        return repr->y0;
    }
    float Surface::z0() const {
        return repr->z0;
    }
    float Surface::max_curv() const {
        return repr->max_curv;
    }
    float Surface::min_curv() const {
        return repr->min_curv;
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
    int Surface::zeros() const {
        return repr->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return repr->hemisphere;
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return repr->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return repr->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return repr->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return repr->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return repr->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::nlabels() const {
        return repr->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return repr->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        repr->strips[i] = to;
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
    void Surface::set_x0(float to) {  //  center of spherical expansion
        repr->x0 = to;
    }
    void Surface::set_y0(float to) {
        repr->y0 = to;
    }
    void Surface::set_z0(float to) {
        repr->z0 = to;
    }
    void Surface::set_max_curv(float to) {
        repr->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        repr->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        repr->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        repr->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        repr->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        repr->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        repr->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        repr->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        repr->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        repr->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        repr->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->f_gradNorm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    float Vertex::pnx() const {
        return repr->v_pnx[idx];
    }
    float Vertex::pny() const {
        return repr->v_pny[idx];
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->v_pnz[idx];
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
    float Vertex::onx() const {
        return repr->v_onx[idx];
    }
    float Vertex::ony() const {
        return repr->v_ony[idx];
    }
    float Vertex::onz() const {  //  original normal
        return repr->v_onz[idx];
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
    float Vertex::odx() const {
        return repr->v_odx[idx];
    }
    float Vertex::ody() const {
        return repr->v_ody[idx];
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->v_odz[idx];
    }
    float Vertex::tdx() const {
        return repr->v_tdx[idx];
    }
    float Vertex::tdy() const {
        return repr->v_tdy[idx];
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->v_tdz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::curvbak() const {
        return repr->v_curvbak[idx];
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->v_val[idx];
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->v_imag_val[idx];
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
    float Vertex::tx() const {
        return repr->v_tx[idx];
    }
    float Vertex::ty() const {
        return repr->v_ty[idx];
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->v_tz[idx];
    }
    float Vertex::t2x() const {
        return repr->v_t2x[idx];
    }
    float Vertex::t2y() const {
        return repr->v_t2y[idx];
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->v_t2z[idx];
    }
    float Vertex::targx() const {
        return repr->v_targx[idx];
    }
    float Vertex::targy() const {
        return repr->v_targy[idx];
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->v_targz[idx];
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
    float Vertex::l4x() const {
        return repr->v_l4x[idx];
    }
    float Vertex::l4y() const {
        return repr->v_l4y[idx];
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->v_l4z[idx];
    }
    float Vertex::infx() const {
        return repr->v_infx[idx];
    }
    float Vertex::infy() const {
        return repr->v_infy[idx];
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->v_infz[idx];
    }
    float Vertex::fx() const {
        return repr->v_fx[idx];
    }
    float Vertex::fy() const {
        return repr->v_fy[idx];
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->v_fz[idx];
    }
    int Vertex::px() const {
        return repr->v_px[idx];
    }
    int Vertex::qx() const {
        return repr->v_qx[idx];
    }
    int Vertex::py() const {
        return repr->v_py[idx];
    }
    int Vertex::qy() const {
        return repr->v_qy[idx];
    }
    int Vertex::pz() const {
        return repr->v_pz[idx];
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->v_qz[idx];
    }
    float Vertex::e1x() const {
        return repr->v_e1x[idx];
    }
    float Vertex::e1y() const {
        return repr->v_e1y[idx];
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_e1z[idx];
    }
    float Vertex::e2x() const {
        return repr->v_e2x[idx];
    }
    float Vertex::e2y() const {
        return repr->v_e2y[idx];
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_e2z[idx];
    }
    float Vertex::pe1x() const {
        return repr->v_pe1x[idx];
    }
    float Vertex::pe1y() const {
        return repr->v_pe1y[idx];
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_pe1z[idx];
    }
    float Vertex::pe2x() const {
        return repr->v_pe2x[idx];
    }
    float Vertex::pe2y() const {
        return repr->v_pe2y[idx];
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_pe2z[idx];
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->v_nc[idx];
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->v_val2[idx];
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->v_valbak[idx];
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->v_val2bak[idx];
    }
    float Vertex::stat() const {  //  statistic 
        return repr->v_stat[idx];
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->v_undefval[idx];
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->v_old_undefval[idx];
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->v_fixedval[idx];
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->v_fieldsign[idx];
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->v_fsmask[idx];
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->v_d[idx];
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->v_annotation[idx];
    }
    char Vertex::oripflag() const {
        return repr->v_oripflag[idx];
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->v_origripflag[idx];
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->v_vp[idx];
    }
    float Vertex::theta() const {
        return repr->v_theta[idx];
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->v_phi[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    float Vertex::group_avg_area() const {
        return repr->v_group_avg_area[idx];
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->v_K[idx];
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->v_H[idx];
    }
    float Vertex::k1() const {
        return repr->v_k1[idx];
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->v_k2[idx];
    }
    float Vertex::mean() const {
        return repr->v_mean[idx];
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->v_mean_imag[idx];
    }
    float Vertex::std_error() const {
        return repr->v_std_error[idx];
    }
    uint Vertex::flags() const {
        return repr->v_flags[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    int Vertex::cropped() const {
        return repr->v_cropped[idx];
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->v_marked[idx];
    }
    short Vertex::marked2() const {
        return repr->v_marked2[idx];
    }
    short Vertex::marked3() const {
        return repr->v_marked3[idx];
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
        CASE(PIAL_NORMALS,pn)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(TMP_VERTICES,t)
        CASE(TMP2_VERTICES,t2)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        CASE(INFLATED_VERTICES,inf)
        CASE(FLATTENED_VERTICES,f)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
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
    void Vertex::set_pnx(float to) {
        repr->v_pnx[idx] = to;
    }
    void Vertex::set_pny(float to) {
        repr->v_pny[idx] = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->v_pnz[idx] = to;
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
    void Vertex::set_onx(float to) {
        repr->v_onx[idx] = to;
    }
    void Vertex::set_ony(float to) {
        repr->v_ony[idx] = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->v_onz[idx] = to;
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
    void Vertex::set_odx(float to) {
        repr->v_odx[idx] = to;
    }
    void Vertex::set_ody(float to) {
        repr->v_ody[idx] = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->v_odz[idx] = to;
    }
    void Vertex::set_tdx(float to) {
        repr->v_tdx[idx] = to;
    }
    void Vertex::set_tdy(float to) {
        repr->v_tdy[idx] = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->v_tdz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->v_curvbak[idx] = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->v_val[idx] = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->v_imag_val[idx] = to;
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
    void Vertex::set_tx(float to) {
        repr->v_tx[idx] = to;
    }
    void Vertex::set_ty(float to) {
        repr->v_ty[idx] = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->v_tz[idx] = to;
    }
    void Vertex::set_t2x(float to) {
        repr->v_t2x[idx] = to;
    }
    void Vertex::set_t2y(float to) {
        repr->v_t2y[idx] = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->v_t2z[idx] = to;
    }
    void Vertex::set_targx(float to) {
        repr->v_targx[idx] = to;
    }
    void Vertex::set_targy(float to) {
        repr->v_targy[idx] = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->v_targz[idx] = to;
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
    void Vertex::set_l4x(float to) {
        repr->v_l4x[idx] = to;
    }
    void Vertex::set_l4y(float to) {
        repr->v_l4y[idx] = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->v_l4z[idx] = to;
    }
    void Vertex::set_infx(float to) {
        repr->v_infx[idx] = to;
    }
    void Vertex::set_infy(float to) {
        repr->v_infy[idx] = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->v_infz[idx] = to;
    }
    void Vertex::set_fx(float to) {
        repr->v_fx[idx] = to;
    }
    void Vertex::set_fy(float to) {
        repr->v_fy[idx] = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->v_fz[idx] = to;
    }
    void Vertex::set_px(int to) {
        repr->v_px[idx] = to;
    }
    void Vertex::set_qx(int to) {
        repr->v_qx[idx] = to;
    }
    void Vertex::set_py(int to) {
        repr->v_py[idx] = to;
    }
    void Vertex::set_qy(int to) {
        repr->v_qy[idx] = to;
    }
    void Vertex::set_pz(int to) {
        repr->v_pz[idx] = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->v_qz[idx] = to;
    }
    void Vertex::set_e1x(float to) {
        repr->v_e1x[idx] = to;
    }
    void Vertex::set_e1y(float to) {
        repr->v_e1y[idx] = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_e1z[idx] = to;
    }
    void Vertex::set_e2x(float to) {
        repr->v_e2x[idx] = to;
    }
    void Vertex::set_e2y(float to) {
        repr->v_e2y[idx] = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_e2z[idx] = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->v_pe1x[idx] = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->v_pe1y[idx] = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_pe1z[idx] = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->v_pe2x[idx] = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->v_pe2y[idx] = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_pe2z[idx] = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->v_nc[idx] = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->v_val2[idx] = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->v_valbak[idx] = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->v_val2bak[idx] = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->v_stat[idx] = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->v_undefval[idx] = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->v_old_undefval[idx] = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->v_fixedval[idx] = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->v_fieldsign[idx] = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->v_fsmask[idx] = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->v_d[idx] = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->v_annotation[idx] = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->v_oripflag[idx] = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->v_origripflag[idx] = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->v_vp[idx] = to;
    }
    void Vertex::set_theta(float to) {
        repr->v_theta[idx] = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->v_phi[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->v_group_avg_area[idx] = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->v_K[idx] = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->v_H[idx] = to;
    }
    void Vertex::set_k1(float to) {
        repr->v_k1[idx] = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->v_k2[idx] = to;
    }
    void Vertex::set_mean(float to) {
        repr->v_mean[idx] = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->v_mean_imag[idx] = to;
    }
    void Vertex::set_std_error(float to) {
        repr->v_std_error[idx] = to;
    }
    void Vertex::set_flags(uint to) {
        repr->v_flags[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_cropped(int to) {
        repr->v_cropped[idx] = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->v_marked[idx] = to;
    }
    void Vertex::set_marked2(short to) {
        repr->v_marked2[idx] = to;
    }
    void Vertex::set_marked3(short to) {
        repr->v_marked3[idx] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return repr->strips[i];
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
    float Surface::x0() const {  //  center of spherical expansion
        return repr->x0;
    }
    float Surface::y0() const {
        return repr->y0;
    }
    float Surface::z0() const {
        return repr->z0;
    }
    float Surface::max_curv() const {
        return repr->max_curv;
    }
    float Surface::min_curv() const {
        return repr->min_curv;
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
    int Surface::zeros() const {
        return repr->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return repr->hemisphere;
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return repr->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return repr->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return repr->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return repr->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return repr->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::nlabels() const {
        return repr->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return repr->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        repr->strips[i] = to;
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
    void Surface::set_x0(float to) {  //  center of spherical expansion
        repr->x0 = to;
    }
    void Surface::set_y0(float to) {
        repr->y0 = to;
    }
    void Surface::set_z0(float to) {
        repr->z0 = to;
    }
    void Surface::set_max_curv(float to) {
        repr->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        repr->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        repr->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        repr->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        repr->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        repr->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        repr->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        repr->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        repr->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        repr->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        repr->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_fname(MRIS_fname_t to) {  //  file it was originally loaded from
        repr->fname = to;
    }
    void Surface::set_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane)
        repr->status = to;
    }
    void Surface::set_origxyz_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        repr->origxyz_status = to;
    }
    void Surface::set_patch(int to) {  //  if a patch of the surface
        repr->patch = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    
    void Face::set_v(size_t i, Vertex to) {
        cheapAssert(repr == to.repr); repr->f_v[idx][i] = to.idx;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
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
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
    }
    void Vertex::set_v(size_t i, Vertex to) {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        cheapAssert(repr == to.repr); repr->v_v[idx][i] = to.idx;
    }
    void Vertex::set_vnum(short to) {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        repr->v_vnum[idx] = to;
    }
    void Vertex::set_v2num(short to) {  //  number of 1, or 2-hop neighbors                          
        repr->v_v2num[idx] = to;
    }
    void Vertex::set_v3num(short to) {  //  number of 1,2,or 3-hop neighbors                         
        repr->v_v3num[idx] = to;
    }
    void Vertex::set_vtotal(short to) {  //  total # of neighbors. copy of vnum.nsizeCur              
        repr->v_vtotal[idx] = to;
    }
    void Vertex::set_nsizeMaxClock(short to) {  //  copy of mris->nsizeMaxClock when v#num                   
        repr->v_nsizeMaxClock[idx] = to;
    }
    void Vertex::set_nsizeMax(uchar to) {  //  the max nsize that was used to fill in vnum etc          
        repr->v_nsizeMax[idx] = to;
    }
    void Vertex::set_nsizeCur(uchar to) {  //  index of the current v#num in vtotal                     
        repr->v_nsizeCur[idx] = to;
    }
    void Vertex::set_num(uchar to) {  //  number of neighboring faces                              
        repr->v_num[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->f_gradNorm[idx];
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
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }
    void Face::set_norm(PDMATRIX to) {
        repr->f_norm[idx] = to;
    }
    void Face::set_gradNorm(A3PDMATRIX to) {
        repr->f_gradNorm[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return repr->strips[i];
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
    float Surface::x0() const {  //  center of spherical expansion
        return repr->x0;
    }
    float Surface::y0() const {
        return repr->y0;
    }
    float Surface::z0() const {
        return repr->z0;
    }
    float Surface::max_curv() const {
        return repr->max_curv;
    }
    float Surface::min_curv() const {
        return repr->min_curv;
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
    int Surface::zeros() const {
        return repr->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return repr->hemisphere;
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return repr->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return repr->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return repr->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return repr->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return repr->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::nlabels() const {
        return repr->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return repr->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        repr->strips[i] = to;
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
    void Surface::set_x0(float to) {  //  center of spherical expansion
        repr->x0 = to;
    }
    void Surface::set_y0(float to) {
        repr->y0 = to;
    }
    void Surface::set_z0(float to) {
        repr->z0 = to;
    }
    void Surface::set_max_curv(float to) {
        repr->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        repr->min_curv = to;
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
    void Surface::set_zeros(int to) {
        repr->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        repr->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        repr->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        repr->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        repr->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        repr->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        repr->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        repr->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        repr->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->f_gradNorm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    float Vertex::pnx() const {
        return repr->v_pnx[idx];
    }
    float Vertex::pny() const {
        return repr->v_pny[idx];
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->v_pnz[idx];
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
    float Vertex::onx() const {
        return repr->v_onx[idx];
    }
    float Vertex::ony() const {
        return repr->v_ony[idx];
    }
    float Vertex::onz() const {  //  original normal
        return repr->v_onz[idx];
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
    float Vertex::odx() const {
        return repr->v_odx[idx];
    }
    float Vertex::ody() const {
        return repr->v_ody[idx];
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->v_odz[idx];
    }
    float Vertex::tdx() const {
        return repr->v_tdx[idx];
    }
    float Vertex::tdy() const {
        return repr->v_tdy[idx];
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->v_tdz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::curvbak() const {
        return repr->v_curvbak[idx];
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->v_val[idx];
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->v_imag_val[idx];
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
    float Vertex::tx() const {
        return repr->v_tx[idx];
    }
    float Vertex::ty() const {
        return repr->v_ty[idx];
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->v_tz[idx];
    }
    float Vertex::t2x() const {
        return repr->v_t2x[idx];
    }
    float Vertex::t2y() const {
        return repr->v_t2y[idx];
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->v_t2z[idx];
    }
    float Vertex::targx() const {
        return repr->v_targx[idx];
    }
    float Vertex::targy() const {
        return repr->v_targy[idx];
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->v_targz[idx];
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
    float Vertex::l4x() const {
        return repr->v_l4x[idx];
    }
    float Vertex::l4y() const {
        return repr->v_l4y[idx];
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->v_l4z[idx];
    }
    float Vertex::infx() const {
        return repr->v_infx[idx];
    }
    float Vertex::infy() const {
        return repr->v_infy[idx];
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->v_infz[idx];
    }
    float Vertex::fx() const {
        return repr->v_fx[idx];
    }
    float Vertex::fy() const {
        return repr->v_fy[idx];
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->v_fz[idx];
    }
    int Vertex::px() const {
        return repr->v_px[idx];
    }
    int Vertex::qx() const {
        return repr->v_qx[idx];
    }
    int Vertex::py() const {
        return repr->v_py[idx];
    }
    int Vertex::qy() const {
        return repr->v_qy[idx];
    }
    int Vertex::pz() const {
        return repr->v_pz[idx];
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->v_qz[idx];
    }
    float Vertex::e1x() const {
        return repr->v_e1x[idx];
    }
    float Vertex::e1y() const {
        return repr->v_e1y[idx];
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_e1z[idx];
    }
    float Vertex::e2x() const {
        return repr->v_e2x[idx];
    }
    float Vertex::e2y() const {
        return repr->v_e2y[idx];
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_e2z[idx];
    }
    float Vertex::pe1x() const {
        return repr->v_pe1x[idx];
    }
    float Vertex::pe1y() const {
        return repr->v_pe1y[idx];
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_pe1z[idx];
    }
    float Vertex::pe2x() const {
        return repr->v_pe2x[idx];
    }
    float Vertex::pe2y() const {
        return repr->v_pe2y[idx];
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_pe2z[idx];
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->v_nc[idx];
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->v_val2[idx];
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->v_valbak[idx];
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->v_val2bak[idx];
    }
    float Vertex::stat() const {  //  statistic 
        return repr->v_stat[idx];
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->v_undefval[idx];
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->v_old_undefval[idx];
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->v_fixedval[idx];
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->v_fieldsign[idx];
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->v_fsmask[idx];
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->v_d[idx];
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->v_annotation[idx];
    }
    char Vertex::oripflag() const {
        return repr->v_oripflag[idx];
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->v_origripflag[idx];
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->v_vp[idx];
    }
    float Vertex::theta() const {
        return repr->v_theta[idx];
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->v_phi[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    float Vertex::group_avg_area() const {
        return repr->v_group_avg_area[idx];
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->v_K[idx];
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->v_H[idx];
    }
    float Vertex::k1() const {
        return repr->v_k1[idx];
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->v_k2[idx];
    }
    float Vertex::mean() const {
        return repr->v_mean[idx];
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->v_mean_imag[idx];
    }
    float Vertex::std_error() const {
        return repr->v_std_error[idx];
    }
    uint Vertex::flags() const {
        return repr->v_flags[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    int Vertex::cropped() const {
        return repr->v_cropped[idx];
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->v_marked[idx];
    }
    short Vertex::marked2() const {
        return repr->v_marked2[idx];
    }
    short Vertex::marked3() const {
        return repr->v_marked3[idx];
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
        CASE(PIAL_NORMALS,pn)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(TMP_VERTICES,t)
        CASE(TMP2_VERTICES,t2)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        CASE(INFLATED_VERTICES,inf)
        CASE(FLATTENED_VERTICES,f)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
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
    void Vertex::set_pnx(float to) {
        repr->v_pnx[idx] = to;
    }
    void Vertex::set_pny(float to) {
        repr->v_pny[idx] = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->v_pnz[idx] = to;
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
    void Vertex::set_onx(float to) {
        repr->v_onx[idx] = to;
    }
    void Vertex::set_ony(float to) {
        repr->v_ony[idx] = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->v_onz[idx] = to;
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
    void Vertex::set_odx(float to) {
        repr->v_odx[idx] = to;
    }
    void Vertex::set_ody(float to) {
        repr->v_ody[idx] = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->v_odz[idx] = to;
    }
    void Vertex::set_tdx(float to) {
        repr->v_tdx[idx] = to;
    }
    void Vertex::set_tdy(float to) {
        repr->v_tdy[idx] = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->v_tdz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->v_curvbak[idx] = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->v_val[idx] = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->v_imag_val[idx] = to;
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
    void Vertex::set_tx(float to) {
        repr->v_tx[idx] = to;
    }
    void Vertex::set_ty(float to) {
        repr->v_ty[idx] = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->v_tz[idx] = to;
    }
    void Vertex::set_t2x(float to) {
        repr->v_t2x[idx] = to;
    }
    void Vertex::set_t2y(float to) {
        repr->v_t2y[idx] = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->v_t2z[idx] = to;
    }
    void Vertex::set_targx(float to) {
        repr->v_targx[idx] = to;
    }
    void Vertex::set_targy(float to) {
        repr->v_targy[idx] = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->v_targz[idx] = to;
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
    void Vertex::set_l4x(float to) {
        repr->v_l4x[idx] = to;
    }
    void Vertex::set_l4y(float to) {
        repr->v_l4y[idx] = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->v_l4z[idx] = to;
    }
    void Vertex::set_infx(float to) {
        repr->v_infx[idx] = to;
    }
    void Vertex::set_infy(float to) {
        repr->v_infy[idx] = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->v_infz[idx] = to;
    }
    void Vertex::set_fx(float to) {
        repr->v_fx[idx] = to;
    }
    void Vertex::set_fy(float to) {
        repr->v_fy[idx] = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->v_fz[idx] = to;
    }
    void Vertex::set_px(int to) {
        repr->v_px[idx] = to;
    }
    void Vertex::set_qx(int to) {
        repr->v_qx[idx] = to;
    }
    void Vertex::set_py(int to) {
        repr->v_py[idx] = to;
    }
    void Vertex::set_qy(int to) {
        repr->v_qy[idx] = to;
    }
    void Vertex::set_pz(int to) {
        repr->v_pz[idx] = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->v_qz[idx] = to;
    }
    void Vertex::set_e1x(float to) {
        repr->v_e1x[idx] = to;
    }
    void Vertex::set_e1y(float to) {
        repr->v_e1y[idx] = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_e1z[idx] = to;
    }
    void Vertex::set_e2x(float to) {
        repr->v_e2x[idx] = to;
    }
    void Vertex::set_e2y(float to) {
        repr->v_e2y[idx] = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_e2z[idx] = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->v_pe1x[idx] = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->v_pe1y[idx] = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_pe1z[idx] = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->v_pe2x[idx] = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->v_pe2y[idx] = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_pe2z[idx] = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->v_nc[idx] = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->v_val2[idx] = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->v_valbak[idx] = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->v_val2bak[idx] = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->v_stat[idx] = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->v_undefval[idx] = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->v_old_undefval[idx] = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->v_fixedval[idx] = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->v_fieldsign[idx] = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->v_fsmask[idx] = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->v_d[idx] = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->v_annotation[idx] = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->v_oripflag[idx] = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->v_origripflag[idx] = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->v_vp[idx] = to;
    }
    void Vertex::set_theta(float to) {
        repr->v_theta[idx] = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->v_phi[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->v_group_avg_area[idx] = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->v_K[idx] = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->v_H[idx] = to;
    }
    void Vertex::set_k1(float to) {
        repr->v_k1[idx] = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->v_k2[idx] = to;
    }
    void Vertex::set_mean(float to) {
        repr->v_mean[idx] = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->v_mean_imag[idx] = to;
    }
    void Vertex::set_std_error(float to) {
        repr->v_std_error[idx] = to;
    }
    void Vertex::set_flags(uint to) {
        repr->v_flags[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_cropped(int to) {
        repr->v_cropped[idx] = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->v_marked[idx] = to;
    }
    void Vertex::set_marked2(short to) {
        repr->v_marked2[idx] = to;
    }
    void Vertex::set_marked3(short to) {
        repr->v_marked3[idx] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return repr->strips[i];
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
    float Surface::x0() const {  //  center of spherical expansion
        return repr->x0;
    }
    float Surface::y0() const {
        return repr->y0;
    }
    float Surface::z0() const {
        return repr->z0;
    }
    float Surface::max_curv() const {
        return repr->max_curv;
    }
    float Surface::min_curv() const {
        return repr->min_curv;
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
    int Surface::zeros() const {
        return repr->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return repr->hemisphere;
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return repr->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return repr->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return repr->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return repr->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return repr->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::nlabels() const {
        return repr->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return repr->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        repr->strips[i] = to;
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
    void Surface::set_x0(float to) {  //  center of spherical expansion
        repr->x0 = to;
    }
    void Surface::set_y0(float to) {
        repr->y0 = to;
    }
    void Surface::set_z0(float to) {
        repr->z0 = to;
    }
    void Surface::set_max_curv(float to) {
        repr->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        repr->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        repr->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        repr->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        repr->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        repr->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        repr->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        repr->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        repr->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        repr->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        repr->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->f_gradNorm[idx];
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->f_orig_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    float Vertex::pnx() const {
        return repr->v_pnx[idx];
    }
    float Vertex::pny() const {
        return repr->v_pny[idx];
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->v_pnz[idx];
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
    float Vertex::onx() const {
        return repr->v_onx[idx];
    }
    float Vertex::ony() const {
        return repr->v_ony[idx];
    }
    float Vertex::onz() const {  //  original normal
        return repr->v_onz[idx];
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
    float Vertex::odx() const {
        return repr->v_odx[idx];
    }
    float Vertex::ody() const {
        return repr->v_ody[idx];
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->v_odz[idx];
    }
    float Vertex::tdx() const {
        return repr->v_tdx[idx];
    }
    float Vertex::tdy() const {
        return repr->v_tdy[idx];
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->v_tdz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::curvbak() const {
        return repr->v_curvbak[idx];
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->v_val[idx];
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->v_imag_val[idx];
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
    float Vertex::tx() const {
        return repr->v_tx[idx];
    }
    float Vertex::ty() const {
        return repr->v_ty[idx];
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->v_tz[idx];
    }
    float Vertex::t2x() const {
        return repr->v_t2x[idx];
    }
    float Vertex::t2y() const {
        return repr->v_t2y[idx];
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->v_t2z[idx];
    }
    float Vertex::targx() const {
        return repr->v_targx[idx];
    }
    float Vertex::targy() const {
        return repr->v_targy[idx];
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->v_targz[idx];
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
    float Vertex::l4x() const {
        return repr->v_l4x[idx];
    }
    float Vertex::l4y() const {
        return repr->v_l4y[idx];
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->v_l4z[idx];
    }
    float Vertex::infx() const {
        return repr->v_infx[idx];
    }
    float Vertex::infy() const {
        return repr->v_infy[idx];
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->v_infz[idx];
    }
    float Vertex::fx() const {
        return repr->v_fx[idx];
    }
    float Vertex::fy() const {
        return repr->v_fy[idx];
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->v_fz[idx];
    }
    int Vertex::px() const {
        return repr->v_px[idx];
    }
    int Vertex::qx() const {
        return repr->v_qx[idx];
    }
    int Vertex::py() const {
        return repr->v_py[idx];
    }
    int Vertex::qy() const {
        return repr->v_qy[idx];
    }
    int Vertex::pz() const {
        return repr->v_pz[idx];
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->v_qz[idx];
    }
    float Vertex::e1x() const {
        return repr->v_e1x[idx];
    }
    float Vertex::e1y() const {
        return repr->v_e1y[idx];
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_e1z[idx];
    }
    float Vertex::e2x() const {
        return repr->v_e2x[idx];
    }
    float Vertex::e2y() const {
        return repr->v_e2y[idx];
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_e2z[idx];
    }
    float Vertex::pe1x() const {
        return repr->v_pe1x[idx];
    }
    float Vertex::pe1y() const {
        return repr->v_pe1y[idx];
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_pe1z[idx];
    }
    float Vertex::pe2x() const {
        return repr->v_pe2x[idx];
    }
    float Vertex::pe2y() const {
        return repr->v_pe2y[idx];
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_pe2z[idx];
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->v_nc[idx];
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->v_val2[idx];
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->v_valbak[idx];
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->v_val2bak[idx];
    }
    float Vertex::stat() const {  //  statistic 
        return repr->v_stat[idx];
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->v_undefval[idx];
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->v_old_undefval[idx];
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->v_fixedval[idx];
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->v_fieldsign[idx];
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->v_fsmask[idx];
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->v_d[idx];
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->v_annotation[idx];
    }
    char Vertex::oripflag() const {
        return repr->v_oripflag[idx];
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->v_origripflag[idx];
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->v_vp[idx];
    }
    float Vertex::theta() const {
        return repr->v_theta[idx];
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->v_phi[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    float Vertex::group_avg_area() const {
        return repr->v_group_avg_area[idx];
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->v_K[idx];
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->v_H[idx];
    }
    float Vertex::k1() const {
        return repr->v_k1[idx];
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->v_k2[idx];
    }
    float Vertex::mean() const {
        return repr->v_mean[idx];
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->v_mean_imag[idx];
    }
    float Vertex::std_error() const {
        return repr->v_std_error[idx];
    }
    uint Vertex::flags() const {
        return repr->v_flags[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    int Vertex::cropped() const {
        return repr->v_cropped[idx];
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->v_marked[idx];
    }
    short Vertex::marked2() const {
        return repr->v_marked2[idx];
    }
    short Vertex::marked3() const {
        return repr->v_marked3[idx];
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
        CASE(PIAL_NORMALS,pn)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(TMP_VERTICES,t)
        CASE(TMP2_VERTICES,t2)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        CASE(INFLATED_VERTICES,inf)
        CASE(FLATTENED_VERTICES,f)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
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
    void Vertex::set_pnx(float to) {
        repr->v_pnx[idx] = to;
    }
    void Vertex::set_pny(float to) {
        repr->v_pny[idx] = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->v_pnz[idx] = to;
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
    void Vertex::set_onx(float to) {
        repr->v_onx[idx] = to;
    }
    void Vertex::set_ony(float to) {
        repr->v_ony[idx] = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->v_onz[idx] = to;
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
    void Vertex::set_odx(float to) {
        repr->v_odx[idx] = to;
    }
    void Vertex::set_ody(float to) {
        repr->v_ody[idx] = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->v_odz[idx] = to;
    }
    void Vertex::set_tdx(float to) {
        repr->v_tdx[idx] = to;
    }
    void Vertex::set_tdy(float to) {
        repr->v_tdy[idx] = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->v_tdz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->v_curvbak[idx] = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->v_val[idx] = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->v_imag_val[idx] = to;
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
    void Vertex::set_tx(float to) {
        repr->v_tx[idx] = to;
    }
    void Vertex::set_ty(float to) {
        repr->v_ty[idx] = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->v_tz[idx] = to;
    }
    void Vertex::set_t2x(float to) {
        repr->v_t2x[idx] = to;
    }
    void Vertex::set_t2y(float to) {
        repr->v_t2y[idx] = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->v_t2z[idx] = to;
    }
    void Vertex::set_targx(float to) {
        repr->v_targx[idx] = to;
    }
    void Vertex::set_targy(float to) {
        repr->v_targy[idx] = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->v_targz[idx] = to;
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
    void Vertex::set_l4x(float to) {
        repr->v_l4x[idx] = to;
    }
    void Vertex::set_l4y(float to) {
        repr->v_l4y[idx] = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->v_l4z[idx] = to;
    }
    void Vertex::set_infx(float to) {
        repr->v_infx[idx] = to;
    }
    void Vertex::set_infy(float to) {
        repr->v_infy[idx] = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->v_infz[idx] = to;
    }
    void Vertex::set_fx(float to) {
        repr->v_fx[idx] = to;
    }
    void Vertex::set_fy(float to) {
        repr->v_fy[idx] = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->v_fz[idx] = to;
    }
    void Vertex::set_px(int to) {
        repr->v_px[idx] = to;
    }
    void Vertex::set_qx(int to) {
        repr->v_qx[idx] = to;
    }
    void Vertex::set_py(int to) {
        repr->v_py[idx] = to;
    }
    void Vertex::set_qy(int to) {
        repr->v_qy[idx] = to;
    }
    void Vertex::set_pz(int to) {
        repr->v_pz[idx] = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->v_qz[idx] = to;
    }
    void Vertex::set_e1x(float to) {
        repr->v_e1x[idx] = to;
    }
    void Vertex::set_e1y(float to) {
        repr->v_e1y[idx] = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_e1z[idx] = to;
    }
    void Vertex::set_e2x(float to) {
        repr->v_e2x[idx] = to;
    }
    void Vertex::set_e2y(float to) {
        repr->v_e2y[idx] = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_e2z[idx] = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->v_pe1x[idx] = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->v_pe1y[idx] = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_pe1z[idx] = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->v_pe2x[idx] = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->v_pe2y[idx] = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_pe2z[idx] = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->v_nc[idx] = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->v_val2[idx] = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->v_valbak[idx] = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->v_val2bak[idx] = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->v_stat[idx] = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->v_undefval[idx] = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->v_old_undefval[idx] = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->v_fixedval[idx] = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->v_fieldsign[idx] = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->v_fsmask[idx] = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->v_d[idx] = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->v_annotation[idx] = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->v_oripflag[idx] = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->v_origripflag[idx] = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->v_vp[idx] = to;
    }
    void Vertex::set_theta(float to) {
        repr->v_theta[idx] = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->v_phi[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->v_group_avg_area[idx] = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->v_K[idx] = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->v_H[idx] = to;
    }
    void Vertex::set_k1(float to) {
        repr->v_k1[idx] = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->v_k2[idx] = to;
    }
    void Vertex::set_mean(float to) {
        repr->v_mean[idx] = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->v_mean_imag[idx] = to;
    }
    void Vertex::set_std_error(float to) {
        repr->v_std_error[idx] = to;
    }
    void Vertex::set_flags(uint to) {
        repr->v_flags[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_cropped(int to) {
        repr->v_cropped[idx] = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->v_marked[idx] = to;
    }
    void Vertex::set_marked2(short to) {
        repr->v_marked2[idx] = to;
    }
    void Vertex::set_marked3(short to) {
        repr->v_marked3[idx] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return repr->strips[i];
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
    float Surface::x0() const {  //  center of spherical expansion
        return repr->x0;
    }
    float Surface::y0() const {
        return repr->y0;
    }
    float Surface::z0() const {
        return repr->z0;
    }
    float Surface::max_curv() const {
        return repr->max_curv;
    }
    float Surface::min_curv() const {
        return repr->min_curv;
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
    int Surface::zeros() const {
        return repr->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return repr->hemisphere;
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return repr->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return repr->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return repr->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return repr->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return repr->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::nlabels() const {
        return repr->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return repr->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        repr->strips[i] = to;
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
    void Surface::set_x0(float to) {  //  center of spherical expansion
        repr->x0 = to;
    }
    void Surface::set_y0(float to) {
        repr->y0 = to;
    }
    void Surface::set_z0(float to) {
        repr->z0 = to;
    }
    void Surface::set_max_curv(float to) {
        repr->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        repr->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        repr->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        repr->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        repr->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        repr->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        repr->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        repr->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        repr->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        repr->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        repr->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        repr->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        repr->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
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
    char Face::oripflag() const {
        return repr->f_oripflag[idx];
    }
    int Face::marked() const {
        return repr->f_marked[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->f_gradNorm[idx];
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
    void Face::set_oripflag(char to) {
        repr->f_oripflag[idx] = to;
    }
    void Face::set_marked(int to) {
        repr->f_marked[idx] = to;
    }
    void Face::set_norm(PDMATRIX to) {
        repr->f_norm[idx] = to;
    }
    void Face::set_gradNorm(A3PDMATRIX to) {
        repr->f_gradNorm[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr, repr->v_f[idx][i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->v_n[idx][i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->v_e[idx][i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->v_v[idx][i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->v_vnum[idx];
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->v_v2num[idx];
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->v_v3num[idx];
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->v_vtotal[idx];
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->v_nsizeMaxClock[idx];
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->v_nsizeMax[idx];
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->v_nsizeCur[idx];
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->v_dist_orig[idx][i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_capacity[idx];
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->v_dist_orig_capacity[idx];
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
    float Vertex::pnx() const {
        return repr->v_pnx[idx];
    }
    float Vertex::pny() const {
        return repr->v_pny[idx];
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->v_pnz[idx];
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
    float Vertex::onx() const {
        return repr->v_onx[idx];
    }
    float Vertex::ony() const {
        return repr->v_ony[idx];
    }
    float Vertex::onz() const {  //  original normal
        return repr->v_onz[idx];
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
    float Vertex::odx() const {
        return repr->v_odx[idx];
    }
    float Vertex::ody() const {
        return repr->v_ody[idx];
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->v_odz[idx];
    }
    float Vertex::tdx() const {
        return repr->v_tdx[idx];
    }
    float Vertex::tdy() const {
        return repr->v_tdy[idx];
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->v_tdz[idx];
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->v_curv[idx];
    }
    float Vertex::curvbak() const {
        return repr->v_curvbak[idx];
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->v_val[idx];
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->v_imag_val[idx];
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
    float Vertex::tx() const {
        return repr->v_tx[idx];
    }
    float Vertex::ty() const {
        return repr->v_ty[idx];
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->v_tz[idx];
    }
    float Vertex::t2x() const {
        return repr->v_t2x[idx];
    }
    float Vertex::t2y() const {
        return repr->v_t2y[idx];
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->v_t2z[idx];
    }
    float Vertex::targx() const {
        return repr->v_targx[idx];
    }
    float Vertex::targy() const {
        return repr->v_targy[idx];
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->v_targz[idx];
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
    float Vertex::l4x() const {
        return repr->v_l4x[idx];
    }
    float Vertex::l4y() const {
        return repr->v_l4y[idx];
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->v_l4z[idx];
    }
    float Vertex::infx() const {
        return repr->v_infx[idx];
    }
    float Vertex::infy() const {
        return repr->v_infy[idx];
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->v_infz[idx];
    }
    float Vertex::fx() const {
        return repr->v_fx[idx];
    }
    float Vertex::fy() const {
        return repr->v_fy[idx];
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->v_fz[idx];
    }
    int Vertex::px() const {
        return repr->v_px[idx];
    }
    int Vertex::qx() const {
        return repr->v_qx[idx];
    }
    int Vertex::py() const {
        return repr->v_py[idx];
    }
    int Vertex::qy() const {
        return repr->v_qy[idx];
    }
    int Vertex::pz() const {
        return repr->v_pz[idx];
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->v_qz[idx];
    }
    float Vertex::e1x() const {
        return repr->v_e1x[idx];
    }
    float Vertex::e1y() const {
        return repr->v_e1y[idx];
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_e1z[idx];
    }
    float Vertex::e2x() const {
        return repr->v_e2x[idx];
    }
    float Vertex::e2y() const {
        return repr->v_e2y[idx];
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_e2z[idx];
    }
    float Vertex::pe1x() const {
        return repr->v_pe1x[idx];
    }
    float Vertex::pe1y() const {
        return repr->v_pe1y[idx];
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->v_pe1z[idx];
    }
    float Vertex::pe2x() const {
        return repr->v_pe2x[idx];
    }
    float Vertex::pe2y() const {
        return repr->v_pe2y[idx];
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->v_pe2z[idx];
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->v_nc[idx];
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->v_val2[idx];
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->v_valbak[idx];
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->v_val2bak[idx];
    }
    float Vertex::stat() const {  //  statistic 
        return repr->v_stat[idx];
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->v_undefval[idx];
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->v_old_undefval[idx];
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->v_fixedval[idx];
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->v_fieldsign[idx];
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->v_fsmask[idx];
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->v_d[idx];
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->v_annotation[idx];
    }
    char Vertex::oripflag() const {
        return repr->v_oripflag[idx];
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->v_origripflag[idx];
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->v_vp[idx];
    }
    float Vertex::theta() const {
        return repr->v_theta[idx];
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->v_phi[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
    }
    float Vertex::group_avg_area() const {
        return repr->v_group_avg_area[idx];
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->v_K[idx];
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->v_H[idx];
    }
    float Vertex::k1() const {
        return repr->v_k1[idx];
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->v_k2[idx];
    }
    float Vertex::mean() const {
        return repr->v_mean[idx];
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->v_mean_imag[idx];
    }
    float Vertex::std_error() const {
        return repr->v_std_error[idx];
    }
    uint Vertex::flags() const {
        return repr->v_flags[idx];
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->v_fno[idx];
    }
    int Vertex::cropped() const {
        return repr->v_cropped[idx];
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->v_marked[idx];
    }
    short Vertex::marked2() const {
        return repr->v_marked2[idx];
    }
    short Vertex::marked3() const {
        return repr->v_marked3[idx];
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
        CASE(PIAL_NORMALS,pn)
        CASE(WHITE_NORMALS,wn)
        CASE(CANONICAL_VERTICES,c)
        CASE(TMP_VERTICES,t)
        CASE(TMP2_VERTICES,t2)
        CASE(PIAL_VERTICES,pial)
        CASE(WHITE_VERTICES,white)
        CASE(INFLATED_VERTICES,inf)
        CASE(FLATTENED_VERTICES,f)
        default:
          *x = *y = *z = 0.0;
          ErrorExit(ERROR_UNSUPPORTED, "which_coords: unsupported which %d", which);
          break;
      }
    
    #undef CASE
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->v_f[idx][i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->v_n[idx][i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->v_e[idx][i] = to;
    }
    void Vertex::set_v(size_t i, Vertex to) {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        cheapAssert(repr == to.repr); repr->v_v[idx][i] = to.idx;
    }
    void Vertex::set_vnum(short to) {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        repr->v_vnum[idx] = to;
    }
    void Vertex::set_v2num(short to) {  //  number of 1, or 2-hop neighbors                          
        repr->v_v2num[idx] = to;
    }
    void Vertex::set_v3num(short to) {  //  number of 1,2,or 3-hop neighbors                         
        repr->v_v3num[idx] = to;
    }
    void Vertex::set_vtotal(short to) {  //  total # of neighbors. copy of vnum.nsizeCur              
        repr->v_vtotal[idx] = to;
    }
    void Vertex::set_nsizeMaxClock(short to) {  //  copy of mris->nsizeMaxClock when v#num                   
        repr->v_nsizeMaxClock[idx] = to;
    }
    void Vertex::set_nsizeMax(uchar to) {  //  the max nsize that was used to fill in vnum etc          
        repr->v_nsizeMax[idx] = to;
    }
    void Vertex::set_nsizeCur(uchar to) {  //  index of the current v#num in vtotal                     
        repr->v_nsizeCur[idx] = to;
    }
    void Vertex::set_num(uchar to) {  //  number of neighboring faces                              
        repr->v_num[idx] = to;
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
    void Vertex::set_pnx(float to) {
        repr->v_pnx[idx] = to;
    }
    void Vertex::set_pny(float to) {
        repr->v_pny[idx] = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->v_pnz[idx] = to;
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
    void Vertex::set_onx(float to) {
        repr->v_onx[idx] = to;
    }
    void Vertex::set_ony(float to) {
        repr->v_ony[idx] = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->v_onz[idx] = to;
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
    void Vertex::set_odx(float to) {
        repr->v_odx[idx] = to;
    }
    void Vertex::set_ody(float to) {
        repr->v_ody[idx] = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->v_odz[idx] = to;
    }
    void Vertex::set_tdx(float to) {
        repr->v_tdx[idx] = to;
    }
    void Vertex::set_tdy(float to) {
        repr->v_tdy[idx] = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->v_tdz[idx] = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->v_curv[idx] = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->v_curvbak[idx] = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->v_val[idx] = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->v_imag_val[idx] = to;
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
    void Vertex::set_tx(float to) {
        repr->v_tx[idx] = to;
    }
    void Vertex::set_ty(float to) {
        repr->v_ty[idx] = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->v_tz[idx] = to;
    }
    void Vertex::set_t2x(float to) {
        repr->v_t2x[idx] = to;
    }
    void Vertex::set_t2y(float to) {
        repr->v_t2y[idx] = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->v_t2z[idx] = to;
    }
    void Vertex::set_targx(float to) {
        repr->v_targx[idx] = to;
    }
    void Vertex::set_targy(float to) {
        repr->v_targy[idx] = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->v_targz[idx] = to;
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
    void Vertex::set_l4x(float to) {
        repr->v_l4x[idx] = to;
    }
    void Vertex::set_l4y(float to) {
        repr->v_l4y[idx] = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->v_l4z[idx] = to;
    }
    void Vertex::set_infx(float to) {
        repr->v_infx[idx] = to;
    }
    void Vertex::set_infy(float to) {
        repr->v_infy[idx] = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->v_infz[idx] = to;
    }
    void Vertex::set_fx(float to) {
        repr->v_fx[idx] = to;
    }
    void Vertex::set_fy(float to) {
        repr->v_fy[idx] = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->v_fz[idx] = to;
    }
    void Vertex::set_px(int to) {
        repr->v_px[idx] = to;
    }
    void Vertex::set_qx(int to) {
        repr->v_qx[idx] = to;
    }
    void Vertex::set_py(int to) {
        repr->v_py[idx] = to;
    }
    void Vertex::set_qy(int to) {
        repr->v_qy[idx] = to;
    }
    void Vertex::set_pz(int to) {
        repr->v_pz[idx] = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->v_qz[idx] = to;
    }
    void Vertex::set_e1x(float to) {
        repr->v_e1x[idx] = to;
    }
    void Vertex::set_e1y(float to) {
        repr->v_e1y[idx] = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_e1z[idx] = to;
    }
    void Vertex::set_e2x(float to) {
        repr->v_e2x[idx] = to;
    }
    void Vertex::set_e2y(float to) {
        repr->v_e2y[idx] = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_e2z[idx] = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->v_pe1x[idx] = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->v_pe1y[idx] = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->v_pe1z[idx] = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->v_pe2x[idx] = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->v_pe2y[idx] = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->v_pe2z[idx] = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->v_nc[idx] = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->v_val2[idx] = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->v_valbak[idx] = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->v_val2bak[idx] = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->v_stat[idx] = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->v_undefval[idx] = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->v_old_undefval[idx] = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->v_fixedval[idx] = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->v_fieldsign[idx] = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->v_fsmask[idx] = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->v_d[idx] = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->v_annotation[idx] = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->v_oripflag[idx] = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->v_origripflag[idx] = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->v_vp[idx] = to;
    }
    void Vertex::set_theta(float to) {
        repr->v_theta[idx] = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->v_phi[idx] = to;
    }
    void Vertex::set_area(float to) {
        repr->v_area[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->v_group_avg_area[idx] = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->v_K[idx] = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->v_H[idx] = to;
    }
    void Vertex::set_k1(float to) {
        repr->v_k1[idx] = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->v_k2[idx] = to;
    }
    void Vertex::set_mean(float to) {
        repr->v_mean[idx] = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->v_mean_imag[idx] = to;
    }
    void Vertex::set_std_error(float to) {
        repr->v_std_error[idx] = to;
    }
    void Vertex::set_flags(uint to) {
        repr->v_flags[idx] = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->v_fno[idx] = to;
    }
    void Vertex::set_cropped(int to) {
        repr->v_cropped[idx] = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->v_marked[idx] = to;
    }
    void Vertex::set_marked2(short to) {
        repr->v_marked2[idx] = to;
    }
    void Vertex::set_marked3(short to) {
        repr->v_marked3[idx] = to;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return repr->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return repr->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return repr->nedges;
    }
    int Surface::ncorners() const {  //  # of triangle corners
        return repr->ncorners;
    }
    int Surface::nstrips() const {
        return repr->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr, i);
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return repr->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return repr->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return repr->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(repr, i);
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return repr->edges[i];
    }
    MRI_CORNER Surface::corners(size_t i) const {
        return repr->corners[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return repr->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return repr->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return repr->strips[i];
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
    float Surface::x0() const {  //  center of spherical expansion
        return repr->x0;
    }
    float Surface::y0() const {
        return repr->y0;
    }
    float Surface::z0() const {
        return repr->z0;
    }
    float Surface::max_curv() const {
        return repr->max_curv;
    }
    float Surface::min_curv() const {
        return repr->min_curv;
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
    int Surface::zeros() const {
        return repr->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return repr->hemisphere;
    }
    int Surface::initialized() const {
        return repr->initialized;
    }
    PLTA Surface::lta() const {
        return repr->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return repr->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return repr->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return repr->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    float Surface::a() const {
        return repr->a;
    }
    float Surface::b() const {
        return repr->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return repr->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return repr->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return repr->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return repr->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return repr->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return repr->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return repr->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return repr->patch;
    }
    int Surface::nlabels() const {
        return repr->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return repr->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return repr->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return repr->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return repr->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return repr->gamma;
    }
    float Surface::da() const {
        return repr->da;
    }
    float Surface::db() const {
        return repr->db;
    }
    float Surface::dg() const {  //  old deltas
        return repr->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return repr->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return repr->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return repr->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return repr->subject_name;
    }
    float Surface::canon_area() const {
        return repr->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return repr->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return repr->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return repr->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return repr->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return repr->ct;
    }
    int Surface::orig_xyzspace() const {  //  xyz coordinate space of surface read by MRISread() before any conversion, 0=tkregister space, 1=scanner space
        return repr->orig_xyzspace;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return repr->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return repr->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return repr->cmdlines;
    }
    int Surface::ncmds() const {
        return repr->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return repr->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return repr->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return repr->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return repr->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return repr->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return repr->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return repr->mht;
    }
    p_void Surface::temps() const {
        return repr->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        repr->strips[i] = to;
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
    void Surface::set_x0(float to) {  //  center of spherical expansion
        repr->x0 = to;
    }
    void Surface::set_y0(float to) {
        repr->y0 = to;
    }
    void Surface::set_z0(float to) {
        repr->z0 = to;
    }
    void Surface::set_max_curv(float to) {
        repr->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        repr->min_curv = to;
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
    void Surface::set_zeros(int to) {
        repr->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        repr->hemisphere = to;
    }
    void Surface::set_fname(MRIS_fname_t to) {  //  file it was originally loaded from
        repr->fname = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        repr->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        repr->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        repr->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        repr->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        repr->Ktotal = to;
    }
    void Surface::set_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane)
        repr->status = to;
    }
    void Surface::set_origxyz_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        repr->origxyz_status = to;
    }
    void Surface::set_patch(int to) {  //  if a patch of the surface
        repr->patch = to;
    }
    void Surface::set_nlabels(int to) {
        repr->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        repr->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        repr->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        repr->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        repr->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        repr->gamma = to;
    }
    void Surface::set_da(float to) {
        repr->da = to;
    }
    void Surface::set_db(float to) {
        repr->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        repr->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        repr->type = to;
    }


    } // namespace AllM
} // SurfaceFromMRISPV
