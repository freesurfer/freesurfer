
#pragma once
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
#include "mrisurf_aaa.h"
namespace SurfaceFromMRIS {
    typedef MRIS Representation;


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
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
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
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
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
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
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
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    float Face::area() const {
        return repr->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return repr->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return repr->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
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
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return repr->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return repr->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return repr->vertices[idx].nz;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    float Vertex::area() const {
        return repr->vertices[idx].area;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        repr->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        repr->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->vertices[idx].nz = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
    Vertex Surface::v_temporal_pole() const {
        return Vertex(repr, repr->v_temporal_pole - repr->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(repr, repr->v_frontal_pole - repr->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(repr, repr->v_occipital_pole - repr->vertices);
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
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_temporal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_frontal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_occipital_pole = repr->vertices + to.idx;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    float Face::area() const {
        return repr->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return repr->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return repr->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AnalysisM::Vertex const & src              ) : Repr_Elt(src) {}
    Vertex::Vertex ( Analysis::Vertex const & src               ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return repr->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return repr->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return repr->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return repr->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return repr->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return repr->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return repr->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return repr->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return repr->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return repr->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return repr->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return repr->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return repr->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return repr->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return repr->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return repr->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return repr->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return repr->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return repr->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return repr->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return repr->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->vertices[idx].tz;
    }
    float Vertex::t2x() const {
        return repr->vertices[idx].t2x;
    }
    float Vertex::t2y() const {
        return repr->vertices[idx].t2y;
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->vertices[idx].t2z;
    }
    float Vertex::targx() const {
        return repr->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return repr->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return repr->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return repr->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return repr->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return repr->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return repr->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return repr->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return repr->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return repr->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return repr->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return repr->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->vertices[idx].fz;
    }
    int Vertex::px() const {
        return repr->vertices[idx].px;
    }
    int Vertex::qx() const {
        return repr->vertices[idx].qx;
    }
    int Vertex::py() const {
        return repr->vertices[idx].py;
    }
    int Vertex::qy() const {
        return repr->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return repr->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return repr->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return repr->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return repr->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return repr->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return repr->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return repr->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return repr->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return repr->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return repr->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return repr->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return repr->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->vertices[idx].phi;
    }
    float Vertex::area() const {
        return repr->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return repr->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return repr->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->vertices[idx].H;
    }
    float Vertex::k1() const {
        return repr->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return repr->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return repr->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return repr->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return repr->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return repr->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return repr->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return repr->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        repr->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        repr->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        repr->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        repr->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        repr->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        repr->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        repr->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        repr->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        repr->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        repr->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        repr->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        repr->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        repr->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        repr->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        repr->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        repr->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->vertices[idx].tz = to;
    }
    void Vertex::set_t2x(float to) {
        repr->vertices[idx].t2x = to;
    }
    void Vertex::set_t2y(float to) {
        repr->vertices[idx].t2y = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->vertices[idx].t2z = to;
    }
    void Vertex::set_targx(float to) {
        repr->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        repr->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        repr->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        repr->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        repr->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        repr->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        repr->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        repr->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        repr->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        repr->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        repr->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        repr->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        repr->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        repr->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        repr->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        repr->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        repr->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        repr->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        repr->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        repr->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        repr->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        repr->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->vertices[idx].phi = to;
    }
    void Vertex::set_origarea(float to) {
        repr->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        repr->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        repr->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        repr->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        repr->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        repr->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        repr->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        repr->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
    Vertex Surface::v_temporal_pole() const {
        return Vertex(repr, repr->v_temporal_pole - repr->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(repr, repr->v_frontal_pole - repr->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(repr, repr->v_occipital_pole - repr->vertices);
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
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_temporal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_frontal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_occipital_pole = repr->vertices + to.idx;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    float Face::area() const {
        return repr->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return repr->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return repr->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return repr->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return repr->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return repr->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return repr->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return repr->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return repr->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return repr->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return repr->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return repr->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return repr->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return repr->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return repr->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return repr->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return repr->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return repr->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return repr->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return repr->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return repr->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return repr->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return repr->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return repr->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->vertices[idx].tz;
    }
    float Vertex::t2x() const {
        return repr->vertices[idx].t2x;
    }
    float Vertex::t2y() const {
        return repr->vertices[idx].t2y;
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->vertices[idx].t2z;
    }
    float Vertex::targx() const {
        return repr->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return repr->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return repr->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return repr->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return repr->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return repr->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return repr->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return repr->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return repr->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return repr->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return repr->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return repr->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->vertices[idx].fz;
    }
    int Vertex::px() const {
        return repr->vertices[idx].px;
    }
    int Vertex::qx() const {
        return repr->vertices[idx].qx;
    }
    int Vertex::py() const {
        return repr->vertices[idx].py;
    }
    int Vertex::qy() const {
        return repr->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return repr->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return repr->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return repr->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return repr->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return repr->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return repr->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return repr->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return repr->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return repr->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return repr->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return repr->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return repr->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->vertices[idx].phi;
    }
    float Vertex::area() const {
        return repr->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return repr->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return repr->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->vertices[idx].H;
    }
    float Vertex::k1() const {
        return repr->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return repr->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return repr->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return repr->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return repr->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return repr->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return repr->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return repr->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        repr->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        repr->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        repr->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        repr->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        repr->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        repr->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        repr->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        repr->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        repr->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        repr->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        repr->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        repr->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        repr->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        repr->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        repr->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        repr->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->vertices[idx].tz = to;
    }
    void Vertex::set_t2x(float to) {
        repr->vertices[idx].t2x = to;
    }
    void Vertex::set_t2y(float to) {
        repr->vertices[idx].t2y = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->vertices[idx].t2z = to;
    }
    void Vertex::set_targx(float to) {
        repr->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        repr->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        repr->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        repr->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        repr->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        repr->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        repr->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        repr->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        repr->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        repr->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        repr->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        repr->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        repr->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        repr->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        repr->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        repr->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        repr->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        repr->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        repr->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        repr->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        repr->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        repr->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->vertices[idx].phi = to;
    }
    void Vertex::set_origarea(float to) {
        repr->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        repr->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        repr->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        repr->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        repr->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        repr->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        repr->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        repr->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
    Vertex Surface::v_temporal_pole() const {
        return Vertex(repr, repr->v_temporal_pole - repr->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(repr, repr->v_frontal_pole - repr->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(repr, repr->v_occipital_pole - repr->vertices);
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
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_temporal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_frontal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_occipital_pole = repr->vertices + to.idx;
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
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    
    void Face::set_v(size_t i, Vertex to) {
        cheapAssert(repr == to.repr); repr->faces[idx].v[i] = to.idx;
    }
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_v(size_t i, Vertex to) {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].v[i] = to.idx;
    }
    void Vertex::set_vnum(short to) {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        repr->vertices_topology[idx].vnum = to;
    }
    void Vertex::set_v2num(short to) {  //  number of 1, or 2-hop neighbors                          
        repr->vertices_topology[idx].v2num = to;
    }
    void Vertex::set_v3num(short to) {  //  number of 1,2,or 3-hop neighbors                         
        repr->vertices_topology[idx].v3num = to;
    }
    void Vertex::set_vtotal(short to) {  //  total # of neighbors. copy of vnum.nsizeCur              
        repr->vertices_topology[idx].vtotal = to;
    }
    void Vertex::set_nsizeMaxClock(short to) {  //  copy of mris->nsizeMaxClock when v#num                   
        repr->vertices_topology[idx].nsizeMaxClock = to;
    }
    void Vertex::set_nsizeMax(uchar to) {  //  the max nsize that was used to fill in vnum etc          
        repr->vertices_topology[idx].nsizeMax = to;
    }
    void Vertex::set_nsizeCur(uchar to) {  //  index of the current v#num in vtotal                     
        repr->vertices_topology[idx].nsizeCur = to;
    }
    void Vertex::set_num(uchar to) {  //  number of neighboring faces                              
        repr->vertices_topology[idx].num = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_x(float to) {  //  current coordinates	
        repr->vertices[idx].x = to;
    }
    void Vertex::set_y(float to) {  //  use MRISsetXYZ() to set
        repr->vertices[idx].y = to;
    }
    void Vertex::set_z(float to) {
        repr->vertices[idx].z = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    float Face::area() const {
        return repr->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return repr->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return repr->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->faces[idx].gradNorm;
    }
    
    void Face::set_area(float to) {
        repr->faces[idx].area = to;
    }
    void Face::set_angle(angles_per_triangle_t to) {
        repr->faces[idx].angle = to;
    }
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }
    void Face::set_norm(PDMATRIX to) {
        repr->faces[idx].norm = to;
    }
    void Face::set_gradNorm(A3PDMATRIX to) {
        repr->faces[idx].gradNorm = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return repr->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return repr->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return repr->vertices[idx].nz;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    float Vertex::area() const {
        return repr->vertices[idx].area;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        repr->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        repr->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->vertices[idx].nz = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_area(float to) {
        repr->vertices[idx].area = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
    Vertex Surface::v_temporal_pole() const {
        return Vertex(repr, repr->v_temporal_pole - repr->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(repr, repr->v_frontal_pole - repr->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(repr, repr->v_occipital_pole - repr->vertices);
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
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_temporal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_frontal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_occipital_pole = repr->vertices + to.idx;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    float Face::area() const {
        return repr->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return repr->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return repr->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return repr->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return repr->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return repr->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return repr->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return repr->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return repr->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return repr->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return repr->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return repr->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return repr->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return repr->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return repr->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return repr->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return repr->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return repr->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return repr->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return repr->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return repr->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return repr->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return repr->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return repr->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->vertices[idx].tz;
    }
    float Vertex::t2x() const {
        return repr->vertices[idx].t2x;
    }
    float Vertex::t2y() const {
        return repr->vertices[idx].t2y;
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->vertices[idx].t2z;
    }
    float Vertex::targx() const {
        return repr->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return repr->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return repr->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return repr->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return repr->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return repr->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return repr->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return repr->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return repr->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return repr->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return repr->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return repr->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->vertices[idx].fz;
    }
    int Vertex::px() const {
        return repr->vertices[idx].px;
    }
    int Vertex::qx() const {
        return repr->vertices[idx].qx;
    }
    int Vertex::py() const {
        return repr->vertices[idx].py;
    }
    int Vertex::qy() const {
        return repr->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return repr->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return repr->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return repr->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return repr->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return repr->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return repr->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return repr->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return repr->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return repr->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return repr->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return repr->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return repr->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->vertices[idx].phi;
    }
    float Vertex::area() const {
        return repr->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return repr->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return repr->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->vertices[idx].H;
    }
    float Vertex::k1() const {
        return repr->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return repr->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return repr->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return repr->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return repr->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return repr->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return repr->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return repr->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        repr->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        repr->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        repr->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        repr->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        repr->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        repr->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        repr->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        repr->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        repr->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        repr->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        repr->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        repr->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        repr->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        repr->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        repr->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        repr->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->vertices[idx].tz = to;
    }
    void Vertex::set_t2x(float to) {
        repr->vertices[idx].t2x = to;
    }
    void Vertex::set_t2y(float to) {
        repr->vertices[idx].t2y = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->vertices[idx].t2z = to;
    }
    void Vertex::set_targx(float to) {
        repr->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        repr->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        repr->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        repr->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        repr->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        repr->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        repr->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        repr->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        repr->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        repr->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        repr->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        repr->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        repr->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        repr->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        repr->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        repr->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        repr->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        repr->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        repr->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        repr->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        repr->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        repr->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->vertices[idx].phi = to;
    }
    void Vertex::set_origarea(float to) {
        repr->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        repr->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        repr->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        repr->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        repr->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        repr->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        repr->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        repr->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
    Vertex Surface::v_temporal_pole() const {
        return Vertex(repr, repr->v_temporal_pole - repr->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(repr, repr->v_frontal_pole - repr->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(repr, repr->v_occipital_pole - repr->vertices);
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
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_temporal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_frontal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_occipital_pole = repr->vertices + to.idx;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    float Face::area() const {
        return repr->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return repr->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return repr->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return repr->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return repr->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return repr->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return repr->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return repr->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return repr->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return repr->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return repr->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return repr->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return repr->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return repr->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return repr->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return repr->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return repr->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return repr->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return repr->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return repr->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return repr->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return repr->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return repr->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return repr->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->vertices[idx].tz;
    }
    float Vertex::t2x() const {
        return repr->vertices[idx].t2x;
    }
    float Vertex::t2y() const {
        return repr->vertices[idx].t2y;
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->vertices[idx].t2z;
    }
    float Vertex::targx() const {
        return repr->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return repr->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return repr->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return repr->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return repr->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return repr->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return repr->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return repr->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return repr->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return repr->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return repr->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return repr->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->vertices[idx].fz;
    }
    int Vertex::px() const {
        return repr->vertices[idx].px;
    }
    int Vertex::qx() const {
        return repr->vertices[idx].qx;
    }
    int Vertex::py() const {
        return repr->vertices[idx].py;
    }
    int Vertex::qy() const {
        return repr->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return repr->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return repr->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return repr->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return repr->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return repr->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return repr->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return repr->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return repr->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return repr->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return repr->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return repr->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return repr->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->vertices[idx].phi;
    }
    float Vertex::area() const {
        return repr->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return repr->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return repr->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->vertices[idx].H;
    }
    float Vertex::k1() const {
        return repr->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return repr->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return repr->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return repr->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return repr->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return repr->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return repr->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return repr->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        repr->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        repr->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        repr->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        repr->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        repr->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        repr->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        repr->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        repr->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        repr->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        repr->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        repr->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        repr->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        repr->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        repr->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        repr->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        repr->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->vertices[idx].tz = to;
    }
    void Vertex::set_t2x(float to) {
        repr->vertices[idx].t2x = to;
    }
    void Vertex::set_t2y(float to) {
        repr->vertices[idx].t2y = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->vertices[idx].t2z = to;
    }
    void Vertex::set_targx(float to) {
        repr->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        repr->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        repr->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        repr->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        repr->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        repr->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        repr->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        repr->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        repr->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        repr->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        repr->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        repr->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        repr->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        repr->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        repr->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        repr->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        repr->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        repr->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        repr->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        repr->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        repr->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        repr->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->vertices[idx].phi = to;
    }
    void Vertex::set_origarea(float to) {
        repr->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        repr->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        repr->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        repr->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        repr->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        repr->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        repr->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        repr->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
    Vertex Surface::v_temporal_pole() const {
        return Vertex(repr, repr->v_temporal_pole - repr->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(repr, repr->v_frontal_pole - repr->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(repr, repr->v_occipital_pole - repr->vertices);
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
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_temporal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_frontal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_occipital_pole = repr->vertices + to.idx;
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
        return Vertex(repr,repr->faces[idx].v[i]);
    }
    float Face::area() const {
        return repr->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return repr->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return repr->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return repr->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return repr->faces[idx].oripflag;
    }
    int Face::marked() const {
        return repr->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return repr->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return repr->faces[idx].gradNorm;
    }
    
    void Face::set_v(size_t i, Vertex to) {
        cheapAssert(repr == to.repr); repr->faces[idx].v[i] = to.idx;
    }
    void Face::set_area(float to) {
        repr->faces[idx].area = to;
    }
    void Face::set_angle(angles_per_triangle_t to) {
        repr->faces[idx].angle = to;
    }
    void Face::set_orig_angle(angles_per_triangle_t to) {
        repr->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        repr->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        repr->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        repr->faces[idx].marked = to;
    }
    void Face::set_norm(PDMATRIX to) {
        repr->faces[idx].norm = to;
    }
    void Face::set_gradNorm(A3PDMATRIX to) {
        repr->faces[idx].gradNorm = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return repr->vertices_topology[idx].n[i];
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return repr->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(repr,repr->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        return repr->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return repr->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return repr->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return repr->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return repr->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return repr->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return repr->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return repr->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return repr->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return repr->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return repr->vertices[idx].y;
    }
    float Vertex::z() const {
        return repr->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return repr->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return repr->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return repr->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return repr->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return repr->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return repr->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return repr->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return repr->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return repr->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return repr->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return repr->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return repr->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return repr->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return repr->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return repr->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return repr->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return repr->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return repr->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return repr->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return repr->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return repr->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return repr->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return repr->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return repr->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return repr->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return repr->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return repr->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return repr->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return repr->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return repr->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return repr->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return repr->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return repr->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return repr->vertices[idx].tz;
    }
    float Vertex::t2x() const {
        return repr->vertices[idx].t2x;
    }
    float Vertex::t2y() const {
        return repr->vertices[idx].t2y;
    }
    float Vertex::t2z() const {  //  another tmp coordinate storage
        return repr->vertices[idx].t2z;
    }
    float Vertex::targx() const {
        return repr->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return repr->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return repr->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return repr->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return repr->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return repr->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return repr->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return repr->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return repr->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return repr->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return repr->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return repr->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return repr->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return repr->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return repr->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return repr->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return repr->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return repr->vertices[idx].fz;
    }
    int Vertex::px() const {
        return repr->vertices[idx].px;
    }
    int Vertex::qx() const {
        return repr->vertices[idx].qx;
    }
    int Vertex::py() const {
        return repr->vertices[idx].py;
    }
    int Vertex::qy() const {
        return repr->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return repr->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return repr->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return repr->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return repr->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return repr->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return repr->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return repr->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return repr->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return repr->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return repr->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return repr->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return repr->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return repr->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return repr->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return repr->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return repr->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return repr->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return repr->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return repr->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return repr->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return repr->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return repr->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return repr->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return repr->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return repr->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return repr->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return repr->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return repr->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return repr->vertices[idx].phi;
    }
    float Vertex::area() const {
        return repr->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return repr->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return repr->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return repr->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return repr->vertices[idx].H;
    }
    float Vertex::k1() const {
        return repr->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return repr->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return repr->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return repr->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return repr->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return repr->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return repr->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return repr->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return repr->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return repr->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return repr->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return repr->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return repr->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return repr->vertices[idx].ripflag;
    }
    void Vertex::which_coords(int which, float *x, float *y, float *z) const {
         MRISvertexCoord2XYZ_float(&repr->vertices[idx], which, x, y, z);
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        repr->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        repr->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_v(size_t i, Vertex to) {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        cheapAssert(repr == to.repr); repr->vertices_topology[idx].v[i] = to.idx;
    }
    void Vertex::set_vnum(short to) {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i) 
        repr->vertices_topology[idx].vnum = to;
    }
    void Vertex::set_v2num(short to) {  //  number of 1, or 2-hop neighbors                          
        repr->vertices_topology[idx].v2num = to;
    }
    void Vertex::set_v3num(short to) {  //  number of 1,2,or 3-hop neighbors                         
        repr->vertices_topology[idx].v3num = to;
    }
    void Vertex::set_vtotal(short to) {  //  total # of neighbors. copy of vnum.nsizeCur              
        repr->vertices_topology[idx].vtotal = to;
    }
    void Vertex::set_nsizeMaxClock(short to) {  //  copy of mris->nsizeMaxClock when v#num                   
        repr->vertices_topology[idx].nsizeMaxClock = to;
    }
    void Vertex::set_nsizeMax(uchar to) {  //  the max nsize that was used to fill in vnum etc          
        repr->vertices_topology[idx].nsizeMax = to;
    }
    void Vertex::set_nsizeCur(uchar to) {  //  index of the current v#num in vtotal                     
        repr->vertices_topology[idx].nsizeCur = to;
    }
    void Vertex::set_num(uchar to) {  //  number of neighboring faces                              
        repr->vertices_topology[idx].num = to;
    }
    void Vertex::set_x(float to) {  //  current coordinates	
        repr->vertices[idx].x = to;
    }
    void Vertex::set_y(float to) {  //  use MRISsetXYZ() to set
        repr->vertices[idx].y = to;
    }
    void Vertex::set_z(float to) {
        repr->vertices[idx].z = to;
    }
    void Vertex::set_nx(float to) {
        repr->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        repr->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        repr->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        repr->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        repr->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        repr->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        repr->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        repr->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        repr->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        repr->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        repr->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        repr->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        repr->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        repr->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        repr->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        repr->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        repr->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        repr->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        repr->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        repr->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        repr->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        repr->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        repr->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        repr->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        repr->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        repr->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        repr->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        repr->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        repr->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        repr->vertices[idx].tz = to;
    }
    void Vertex::set_t2x(float to) {
        repr->vertices[idx].t2x = to;
    }
    void Vertex::set_t2y(float to) {
        repr->vertices[idx].t2y = to;
    }
    void Vertex::set_t2z(float to) {  //  another tmp coordinate storage
        repr->vertices[idx].t2z = to;
    }
    void Vertex::set_targx(float to) {
        repr->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        repr->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        repr->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        repr->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        repr->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        repr->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        repr->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        repr->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        repr->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        repr->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        repr->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        repr->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        repr->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        repr->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        repr->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        repr->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        repr->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        repr->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        repr->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        repr->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        repr->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        repr->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        repr->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        repr->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        repr->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        repr->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        repr->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        repr->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        repr->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        repr->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        repr->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        repr->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        repr->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        repr->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        repr->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        repr->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        repr->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        repr->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        repr->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        repr->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        repr->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        repr->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        repr->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        repr->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        repr->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        repr->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        repr->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        repr->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        repr->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        repr->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        repr->vertices[idx].phi = to;
    }
    void Vertex::set_area(float to) {
        repr->vertices[idx].area = to;
    }
    void Vertex::set_origarea(float to) {
        repr->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        repr->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        repr->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        repr->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        repr->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        repr->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        repr->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        repr->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        repr->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        repr->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        repr->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        repr->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        repr->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        repr->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        repr->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        repr->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        repr->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->vertices[idx].ripflag = to;
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
    Vertex Surface::v_temporal_pole() const {
        return Vertex(repr, repr->v_temporal_pole - repr->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(repr, repr->v_frontal_pole - repr->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(repr, repr->v_occipital_pole - repr->vertices);
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
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_temporal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_frontal_pole = repr->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(repr == to.repr); repr->v_occipital_pole = repr->vertices + to.idx;
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
} // SurfaceFromMRIS
