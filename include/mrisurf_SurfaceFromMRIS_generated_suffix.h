
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
#define SEPARATE_VERTEX_TOPOLOGY
namespace SurfaceFromMRIS {


    namespace Existence {
    Face::Face (                                                                 ) {}
    Face::Face ( MRIS* mris, size_t idx                     ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src                           ) : MRIS_Elt(src) {}     
    Face::Face ( TopologyM::Face const & src                ) : MRIS_Elt(src) {}     
    Face::Face ( Topology::Face const & src                 ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPositionM::Face const & src             ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPosition::Face const & src              ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPositionConsequencesM::Face const & src ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPositionConsequences::Face const & src  ) : MRIS_Elt(src) {}     
    Face::Face ( DistortM::Face const & src                 ) : MRIS_Elt(src) {}     
    Face::Face ( Distort::Face const & src                  ) : MRIS_Elt(src) {}     
    Face::Face ( AnalysisM::Face const & src                ) : MRIS_Elt(src) {}     
    Face::Face ( Analysis::Face const & src                 ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src                     ) : MRIS_Elt(src) {}     

    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                                              ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx                       ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src                           ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( TopologyM::Vertex const & src                ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Topology::Vertex const & src                 ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPositionM::Vertex const & src             ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPosition::Vertex const & src              ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPositionConsequencesM::Vertex const & src ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPositionConsequences::Vertex const & src  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( DistortM::Vertex const & src                 ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Distort::Vertex const & src                  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AnalysisM::Vertex const & src                ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Analysis::Vertex const & src                 ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src                     ) : MRIS_Elt(src) {}     

    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                                               ) {}                   
    Surface::Surface ( MRIS* mris                                    ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src                           ) : MRIS_Elt(src) {}   
    Surface::Surface ( TopologyM::Surface const & src                ) : MRIS_Elt(src) {}   
    Surface::Surface ( Topology::Surface const & src                 ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPositionM::Surface const & src             ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPosition::Surface const & src              ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPositionConsequencesM::Surface const & src ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPositionConsequences::Surface const & src  ) : MRIS_Elt(src) {}   
    Surface::Surface ( DistortM::Surface const & src                 ) : MRIS_Elt(src) {}   
    Surface::Surface ( Distort::Surface const & src                  ) : MRIS_Elt(src) {}   
    Surface::Surface ( AnalysisM::Surface const & src                ) : MRIS_Elt(src) {}   
    Surface::Surface ( Analysis::Surface const & src                 ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src                     ) : MRIS_Elt(src) {}   

    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }


    } // namespace Existence


    namespace Topology {
    Face::Face (                                            ) {}                     
    Face::Face ( MRIS* mris, size_t idx                     ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src                           ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPositionM::Face const & src             ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPosition::Face const & src              ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPositionConsequencesM::Face const & src ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPositionConsequences::Face const & src  ) : MRIS_Elt(src) {}     
    Face::Face ( DistortM::Face const & src                 ) : MRIS_Elt(src) {}     
    Face::Face ( Distort::Face const & src                  ) : MRIS_Elt(src) {}     
    Face::Face ( AnalysisM::Face const & src                ) : MRIS_Elt(src) {}     
    Face::Face ( Analysis::Face const & src                 ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src                     ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                                              ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx                       ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src                           ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPositionM::Vertex const & src             ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPosition::Vertex const & src              ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPositionConsequencesM::Vertex const & src ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPositionConsequences::Vertex const & src  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( DistortM::Vertex const & src                 ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Distort::Vertex const & src                  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AnalysisM::Vertex const & src                ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Analysis::Vertex const & src                 ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src                     ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                                               ) {}                   
    Surface::Surface ( MRIS* mris                                    ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src                           ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPositionM::Surface const & src             ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPosition::Surface const & src              ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPositionConsequencesM::Surface const & src ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPositionConsequences::Surface const & src  ) : MRIS_Elt(src) {}   
    Surface::Surface ( DistortM::Surface const & src                 ) : MRIS_Elt(src) {}   
    Surface::Surface ( Distort::Surface const & src                  ) : MRIS_Elt(src) {}   
    Surface::Surface ( AnalysisM::Surface const & src                ) : MRIS_Elt(src) {}   
    Surface::Surface ( Analysis::Surface const & src                 ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src                     ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }


    } // namespace Topology


    namespace XYZPosition {
    Face::Face (                                            ) {}                     
    Face::Face ( MRIS* mris, size_t idx                     ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src                           ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPositionConsequencesM::Face const & src ) : MRIS_Elt(src) {}     
    Face::Face ( XYZPositionConsequences::Face const & src  ) : MRIS_Elt(src) {}     
    Face::Face ( DistortM::Face const & src                 ) : MRIS_Elt(src) {}     
    Face::Face ( Distort::Face const & src                  ) : MRIS_Elt(src) {}     
    Face::Face ( AnalysisM::Face const & src                ) : MRIS_Elt(src) {}     
    Face::Face ( Analysis::Face const & src                 ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src                     ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                                              ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx                       ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src                           ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPositionConsequencesM::Vertex const & src ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( XYZPositionConsequences::Vertex const & src  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( DistortM::Vertex const & src                 ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Distort::Vertex const & src                  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AnalysisM::Vertex const & src                ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Analysis::Vertex const & src                 ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src                     ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                                               ) {}                   
    Surface::Surface ( MRIS* mris                                    ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src                           ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPositionConsequencesM::Surface const & src ) : MRIS_Elt(src) {}   
    Surface::Surface ( XYZPositionConsequences::Surface const & src  ) : MRIS_Elt(src) {}   
    Surface::Surface ( DistortM::Surface const & src                 ) : MRIS_Elt(src) {}   
    Surface::Surface ( Distort::Surface const & src                  ) : MRIS_Elt(src) {}   
    Surface::Surface ( AnalysisM::Surface const & src                ) : MRIS_Elt(src) {}   
    Surface::Surface ( Analysis::Surface const & src                 ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src                     ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace XYZPosition


    namespace XYZPositionConsequences {
    Face::Face (                             ) {}                     
    Face::Face ( MRIS* mris, size_t idx      ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src            ) : MRIS_Elt(src) {}     
    Face::Face ( DistortM::Face const & src  ) : MRIS_Elt(src) {}     
    Face::Face ( Distort::Face const & src   ) : MRIS_Elt(src) {}     
    Face::Face ( AnalysisM::Face const & src ) : MRIS_Elt(src) {}     
    Face::Face ( Analysis::Face const & src  ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src      ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    float Face::area() const {
        return mris->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return mris->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return mris->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return mris->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return mris->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        mris->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                               ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx        ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src            ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( DistortM::Vertex const & src  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Distort::Vertex const & src   ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AnalysisM::Vertex const & src ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Analysis::Vertex const & src  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src      ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return mris->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return mris->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return mris->vertices[idx].nz;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    float Vertex::area() const {
        return mris->vertices[idx].area;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        mris->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        mris->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        mris->vertices[idx].nz = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                                ) {}                   
    Surface::Surface ( MRIS* mris                     ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src            ) : MRIS_Elt(src) {}   
    Surface::Surface ( DistortM::Surface const & src  ) : MRIS_Elt(src) {}   
    Surface::Surface ( Distort::Surface const & src   ) : MRIS_Elt(src) {}   
    Surface::Surface ( AnalysisM::Surface const & src ) : MRIS_Elt(src) {}   
    Surface::Surface ( Analysis::Surface const & src  ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src      ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return mris->strips[i];
    }
    float Surface::xctr() const {
        return mris->xctr;
    }
    float Surface::yctr() const {
        return mris->yctr;
    }
    float Surface::zctr() const {
        return mris->zctr;
    }
    float Surface::xlo() const {
        return mris->xlo;
    }
    float Surface::ylo() const {
        return mris->ylo;
    }
    float Surface::zlo() const {
        return mris->zlo;
    }
    float Surface::xhi() const {
        return mris->xhi;
    }
    float Surface::yhi() const {
        return mris->yhi;
    }
    float Surface::zhi() const {
        return mris->zhi;
    }
    float Surface::x0() const {  //  center of spherical expansion
        return mris->x0;
    }
    float Surface::y0() const {
        return mris->y0;
    }
    float Surface::z0() const {
        return mris->z0;
    }
    Vertex Surface::v_temporal_pole() const {
        return Vertex(mris,mris->v_temporal_pole - mris->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(mris,mris->v_frontal_pole - mris->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(mris,mris->v_occipital_pole - mris->vertices);
    }
    float Surface::max_curv() const {
        return mris->max_curv;
    }
    float Surface::min_curv() const {
        return mris->min_curv;
    }
    float Surface::total_area() const {
        return mris->total_area;
    }
    double Surface::avg_vertex_area() const {
        return mris->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return mris->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return mris->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return mris->orig_area;
    }
    float Surface::neg_area() const {
        return mris->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return mris->neg_orig_area;
    }
    int Surface::zeros() const {
        return mris->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return mris->hemisphere;
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return mris->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return mris->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return mris->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return mris->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return mris->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::nlabels() const {
        return mris->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return mris->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        mris->strips[i] = to;
    }
    void Surface::set_xctr(float to) {
        mris->xctr = to;
    }
    void Surface::set_yctr(float to) {
        mris->yctr = to;
    }
    void Surface::set_zctr(float to) {
        mris->zctr = to;
    }
    void Surface::set_xlo(float to) {
        mris->xlo = to;
    }
    void Surface::set_ylo(float to) {
        mris->ylo = to;
    }
    void Surface::set_zlo(float to) {
        mris->zlo = to;
    }
    void Surface::set_xhi(float to) {
        mris->xhi = to;
    }
    void Surface::set_yhi(float to) {
        mris->yhi = to;
    }
    void Surface::set_zhi(float to) {
        mris->zhi = to;
    }
    void Surface::set_x0(float to) {  //  center of spherical expansion
        mris->x0 = to;
    }
    void Surface::set_y0(float to) {
        mris->y0 = to;
    }
    void Surface::set_z0(float to) {
        mris->z0 = to;
    }
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_temporal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_frontal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_occipital_pole = mris->vertices + to.idx;
    }
    void Surface::set_max_curv(float to) {
        mris->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        mris->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        mris->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        mris->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        mris->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        mris->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        mris->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        mris->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        mris->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        mris->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        mris->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        mris->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        mris->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace XYZPositionConsequences


    namespace Distort {
    Face::Face (                             ) {}                     
    Face::Face ( MRIS* mris, size_t idx      ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src            ) : MRIS_Elt(src) {}     
    Face::Face ( AnalysisM::Face const & src ) : MRIS_Elt(src) {}     
    Face::Face ( Analysis::Face const & src  ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src      ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    float Face::area() const {
        return mris->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return mris->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return mris->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return mris->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return mris->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        mris->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                               ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx        ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src            ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AnalysisM::Vertex const & src ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( Analysis::Vertex const & src  ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src      ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return mris->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return mris->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return mris->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return mris->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return mris->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return mris->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return mris->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return mris->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return mris->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return mris->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return mris->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return mris->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return mris->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return mris->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return mris->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return mris->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return mris->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return mris->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return mris->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return mris->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return mris->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return mris->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return mris->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return mris->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return mris->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return mris->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return mris->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz;
    }
    float Vertex::tx2() const {
        return mris->vertices[idx].tx2;
    }
    float Vertex::ty2() const {
        return mris->vertices[idx].ty2;
    }
    float Vertex::tz2() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz2;
    }
    float Vertex::targx() const {
        return mris->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return mris->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return mris->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return mris->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return mris->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return mris->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return mris->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return mris->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return mris->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return mris->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return mris->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return mris->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return mris->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return mris->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return mris->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return mris->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return mris->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return mris->vertices[idx].fz;
    }
    int Vertex::px() const {
        return mris->vertices[idx].px;
    }
    int Vertex::qx() const {
        return mris->vertices[idx].qx;
    }
    int Vertex::py() const {
        return mris->vertices[idx].py;
    }
    int Vertex::qy() const {
        return mris->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return mris->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return mris->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return mris->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return mris->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return mris->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return mris->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return mris->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return mris->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return mris->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return mris->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return mris->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return mris->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return mris->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return mris->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return mris->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return mris->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return mris->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return mris->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return mris->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return mris->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return mris->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return mris->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return mris->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return mris->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return mris->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return mris->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return mris->vertices[idx].phi;
    }
    float Vertex::area() const {
        return mris->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return mris->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return mris->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return mris->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return mris->vertices[idx].H;
    }
    float Vertex::k1() const {
        return mris->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return mris->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return mris->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return mris->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return mris->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return mris->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return mris->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return mris->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return mris->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return mris->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return mris->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return mris->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return mris->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        mris->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        mris->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        mris->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        mris->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        mris->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        mris->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        mris->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        mris->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        mris->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        mris->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        mris->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        mris->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        mris->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        mris->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        mris->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        mris->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        mris->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        mris->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        mris->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        mris->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        mris->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        mris->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        mris->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        mris->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        mris->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        mris->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        mris->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz = to;
    }
    void Vertex::set_tx2(float to) {
        mris->vertices[idx].tx2 = to;
    }
    void Vertex::set_ty2(float to) {
        mris->vertices[idx].ty2 = to;
    }
    void Vertex::set_tz2(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz2 = to;
    }
    void Vertex::set_targx(float to) {
        mris->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        mris->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        mris->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        mris->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        mris->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        mris->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        mris->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        mris->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        mris->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        mris->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        mris->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        mris->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        mris->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        mris->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        mris->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        mris->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        mris->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        mris->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        mris->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        mris->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        mris->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        mris->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        mris->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        mris->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        mris->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        mris->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        mris->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        mris->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        mris->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        mris->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        mris->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        mris->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        mris->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        mris->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        mris->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        mris->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        mris->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        mris->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        mris->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        mris->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        mris->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        mris->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        mris->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        mris->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        mris->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        mris->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        mris->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        mris->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        mris->vertices[idx].phi = to;
    }
    void Vertex::set_origarea(float to) {
        mris->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        mris->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        mris->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        mris->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        mris->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        mris->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        mris->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        mris->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        mris->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        mris->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        mris->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        mris->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        mris->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        mris->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        mris->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        mris->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        mris->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                                ) {}                   
    Surface::Surface ( MRIS* mris                     ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src            ) : MRIS_Elt(src) {}   
    Surface::Surface ( AnalysisM::Surface const & src ) : MRIS_Elt(src) {}   
    Surface::Surface ( Analysis::Surface const & src  ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src      ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return mris->strips[i];
    }
    float Surface::xctr() const {
        return mris->xctr;
    }
    float Surface::yctr() const {
        return mris->yctr;
    }
    float Surface::zctr() const {
        return mris->zctr;
    }
    float Surface::xlo() const {
        return mris->xlo;
    }
    float Surface::ylo() const {
        return mris->ylo;
    }
    float Surface::zlo() const {
        return mris->zlo;
    }
    float Surface::xhi() const {
        return mris->xhi;
    }
    float Surface::yhi() const {
        return mris->yhi;
    }
    float Surface::zhi() const {
        return mris->zhi;
    }
    float Surface::x0() const {  //  center of spherical expansion
        return mris->x0;
    }
    float Surface::y0() const {
        return mris->y0;
    }
    float Surface::z0() const {
        return mris->z0;
    }
    Vertex Surface::v_temporal_pole() const {
        return Vertex(mris,mris->v_temporal_pole - mris->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(mris,mris->v_frontal_pole - mris->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(mris,mris->v_occipital_pole - mris->vertices);
    }
    float Surface::max_curv() const {
        return mris->max_curv;
    }
    float Surface::min_curv() const {
        return mris->min_curv;
    }
    float Surface::total_area() const {
        return mris->total_area;
    }
    double Surface::avg_vertex_area() const {
        return mris->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return mris->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return mris->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return mris->orig_area;
    }
    float Surface::neg_area() const {
        return mris->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return mris->neg_orig_area;
    }
    int Surface::zeros() const {
        return mris->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return mris->hemisphere;
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return mris->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return mris->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return mris->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return mris->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return mris->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::nlabels() const {
        return mris->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return mris->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        mris->strips[i] = to;
    }
    void Surface::set_xctr(float to) {
        mris->xctr = to;
    }
    void Surface::set_yctr(float to) {
        mris->yctr = to;
    }
    void Surface::set_zctr(float to) {
        mris->zctr = to;
    }
    void Surface::set_xlo(float to) {
        mris->xlo = to;
    }
    void Surface::set_ylo(float to) {
        mris->ylo = to;
    }
    void Surface::set_zlo(float to) {
        mris->zlo = to;
    }
    void Surface::set_xhi(float to) {
        mris->xhi = to;
    }
    void Surface::set_yhi(float to) {
        mris->yhi = to;
    }
    void Surface::set_zhi(float to) {
        mris->zhi = to;
    }
    void Surface::set_x0(float to) {  //  center of spherical expansion
        mris->x0 = to;
    }
    void Surface::set_y0(float to) {
        mris->y0 = to;
    }
    void Surface::set_z0(float to) {
        mris->z0 = to;
    }
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_temporal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_frontal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_occipital_pole = mris->vertices + to.idx;
    }
    void Surface::set_max_curv(float to) {
        mris->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        mris->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        mris->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        mris->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        mris->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        mris->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        mris->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        mris->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        mris->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        mris->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        mris->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        mris->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        mris->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace Distort


    namespace Analysis {
    Face::Face (                        ) {}                     
    Face::Face ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src       ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    float Face::area() const {
        return mris->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return mris->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return mris->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return mris->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return mris->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        mris->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                          ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx   ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src       ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return mris->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return mris->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return mris->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return mris->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return mris->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return mris->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return mris->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return mris->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return mris->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return mris->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return mris->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return mris->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return mris->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return mris->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return mris->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return mris->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return mris->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return mris->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return mris->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return mris->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return mris->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return mris->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return mris->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return mris->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return mris->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return mris->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return mris->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz;
    }
    float Vertex::tx2() const {
        return mris->vertices[idx].tx2;
    }
    float Vertex::ty2() const {
        return mris->vertices[idx].ty2;
    }
    float Vertex::tz2() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz2;
    }
    float Vertex::targx() const {
        return mris->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return mris->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return mris->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return mris->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return mris->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return mris->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return mris->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return mris->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return mris->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return mris->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return mris->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return mris->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return mris->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return mris->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return mris->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return mris->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return mris->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return mris->vertices[idx].fz;
    }
    int Vertex::px() const {
        return mris->vertices[idx].px;
    }
    int Vertex::qx() const {
        return mris->vertices[idx].qx;
    }
    int Vertex::py() const {
        return mris->vertices[idx].py;
    }
    int Vertex::qy() const {
        return mris->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return mris->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return mris->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return mris->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return mris->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return mris->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return mris->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return mris->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return mris->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return mris->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return mris->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return mris->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return mris->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return mris->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return mris->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return mris->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return mris->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return mris->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return mris->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return mris->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return mris->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return mris->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return mris->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return mris->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return mris->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return mris->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return mris->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return mris->vertices[idx].phi;
    }
    float Vertex::area() const {
        return mris->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return mris->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return mris->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return mris->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return mris->vertices[idx].H;
    }
    float Vertex::k1() const {
        return mris->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return mris->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return mris->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return mris->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return mris->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return mris->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return mris->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return mris->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return mris->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return mris->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return mris->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return mris->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return mris->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        mris->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        mris->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        mris->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        mris->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        mris->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        mris->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        mris->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        mris->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        mris->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        mris->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        mris->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        mris->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        mris->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        mris->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        mris->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        mris->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        mris->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        mris->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        mris->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        mris->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        mris->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        mris->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        mris->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        mris->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        mris->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        mris->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        mris->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz = to;
    }
    void Vertex::set_tx2(float to) {
        mris->vertices[idx].tx2 = to;
    }
    void Vertex::set_ty2(float to) {
        mris->vertices[idx].ty2 = to;
    }
    void Vertex::set_tz2(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz2 = to;
    }
    void Vertex::set_targx(float to) {
        mris->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        mris->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        mris->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        mris->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        mris->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        mris->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        mris->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        mris->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        mris->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        mris->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        mris->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        mris->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        mris->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        mris->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        mris->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        mris->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        mris->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        mris->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        mris->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        mris->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        mris->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        mris->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        mris->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        mris->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        mris->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        mris->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        mris->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        mris->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        mris->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        mris->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        mris->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        mris->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        mris->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        mris->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        mris->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        mris->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        mris->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        mris->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        mris->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        mris->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        mris->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        mris->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        mris->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        mris->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        mris->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        mris->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        mris->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        mris->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        mris->vertices[idx].phi = to;
    }
    void Vertex::set_origarea(float to) {
        mris->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        mris->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        mris->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        mris->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        mris->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        mris->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        mris->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        mris->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        mris->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        mris->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        mris->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        mris->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        mris->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        mris->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        mris->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        mris->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        mris->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                           ) {}                   
    Surface::Surface ( MRIS* mris                ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src       ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return mris->strips[i];
    }
    float Surface::xctr() const {
        return mris->xctr;
    }
    float Surface::yctr() const {
        return mris->yctr;
    }
    float Surface::zctr() const {
        return mris->zctr;
    }
    float Surface::xlo() const {
        return mris->xlo;
    }
    float Surface::ylo() const {
        return mris->ylo;
    }
    float Surface::zlo() const {
        return mris->zlo;
    }
    float Surface::xhi() const {
        return mris->xhi;
    }
    float Surface::yhi() const {
        return mris->yhi;
    }
    float Surface::zhi() const {
        return mris->zhi;
    }
    float Surface::x0() const {  //  center of spherical expansion
        return mris->x0;
    }
    float Surface::y0() const {
        return mris->y0;
    }
    float Surface::z0() const {
        return mris->z0;
    }
    Vertex Surface::v_temporal_pole() const {
        return Vertex(mris,mris->v_temporal_pole - mris->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(mris,mris->v_frontal_pole - mris->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(mris,mris->v_occipital_pole - mris->vertices);
    }
    float Surface::max_curv() const {
        return mris->max_curv;
    }
    float Surface::min_curv() const {
        return mris->min_curv;
    }
    float Surface::total_area() const {
        return mris->total_area;
    }
    double Surface::avg_vertex_area() const {
        return mris->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return mris->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return mris->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return mris->orig_area;
    }
    float Surface::neg_area() const {
        return mris->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return mris->neg_orig_area;
    }
    int Surface::zeros() const {
        return mris->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return mris->hemisphere;
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return mris->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return mris->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return mris->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return mris->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return mris->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::nlabels() const {
        return mris->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return mris->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        mris->strips[i] = to;
    }
    void Surface::set_xctr(float to) {
        mris->xctr = to;
    }
    void Surface::set_yctr(float to) {
        mris->yctr = to;
    }
    void Surface::set_zctr(float to) {
        mris->zctr = to;
    }
    void Surface::set_xlo(float to) {
        mris->xlo = to;
    }
    void Surface::set_ylo(float to) {
        mris->ylo = to;
    }
    void Surface::set_zlo(float to) {
        mris->zlo = to;
    }
    void Surface::set_xhi(float to) {
        mris->xhi = to;
    }
    void Surface::set_yhi(float to) {
        mris->yhi = to;
    }
    void Surface::set_zhi(float to) {
        mris->zhi = to;
    }
    void Surface::set_x0(float to) {  //  center of spherical expansion
        mris->x0 = to;
    }
    void Surface::set_y0(float to) {
        mris->y0 = to;
    }
    void Surface::set_z0(float to) {
        mris->z0 = to;
    }
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_temporal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_frontal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_occipital_pole = mris->vertices + to.idx;
    }
    void Surface::set_max_curv(float to) {
        mris->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        mris->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        mris->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        mris->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        mris->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        mris->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        mris->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        mris->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        mris->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        mris->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        mris->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        mris->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        mris->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace Analysis


    namespace ExistenceM {
    Face::Face (                        ) {}                     
    Face::Face ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src       ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src ) : MRIS_Elt(src) {}     

    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                          ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx   ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src       ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src ) : MRIS_Elt(src) {}     

    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                           ) {}                   
    Surface::Surface ( MRIS* mris                ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src       ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src ) : MRIS_Elt(src) {}   

    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_fname(MRIS_fname_t to) {  //  file it was originally loaded from
        mris->fname = to;
    }
    void Surface::set_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane)
        mris->status = to;
    }
    void Surface::set_origxyz_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        mris->origxyz_status = to;
    }
    void Surface::set_patch(int to) {  //  if a patch of the surface
        mris->patch = to;
    }


    } // namespace ExistenceM


    namespace TopologyM {
    Face::Face (                        ) {}                     
    Face::Face ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src       ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    
    void Face::set_v(vertices_per_face_t to) {
        mris->faces[idx].v = to;
    }
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                          ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx   ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src       ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_v(size_t i, Vertex to) {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].v[i] = to.idx;
    }
    void Vertex::set_vnum(short to) {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        mris->vertices_topology[idx].vnum = to;
    }
    void Vertex::set_v2num(short to) {  //  number of 1, or 2-hop neighbors                          
        mris->vertices_topology[idx].v2num = to;
    }
    void Vertex::set_v3num(short to) {  //  number of 1,2,or 3-hop neighbors                         
        mris->vertices_topology[idx].v3num = to;
    }
    void Vertex::set_vtotal(short to) {  //  total # of neighbors. copy of vnum.nsizeCur              
        mris->vertices_topology[idx].vtotal = to;
    }
    void Vertex::set_nsizeMaxClock(short to) {  //  copy of mris->nsizeMaxClock when v#num                   
        mris->vertices_topology[idx].nsizeMaxClock = to;
    }
    void Vertex::set_nsizeMax(uchar to) {  //  the max nsize that was used to fill in vnum etc          
        mris->vertices_topology[idx].nsizeMax = to;
    }
    void Vertex::set_nsizeCur(uchar to) {  //  index of the current v#num in vtotal                     
        mris->vertices_topology[idx].nsizeCur = to;
    }
    void Vertex::set_num(uchar to) {  //  number of neighboring faces                              
        mris->vertices_topology[idx].num = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                           ) {}                   
    Surface::Surface ( MRIS* mris                ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src       ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }


    } // namespace TopologyM


    namespace XYZPositionM {
    Face::Face (                        ) {}                     
    Face::Face ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src       ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                          ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx   ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src       ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_x(float to) {  //  current coordinates	
        mris->vertices[idx].x = to;
    }
    void Vertex::set_y(float to) {  //  use MRISsetXYZ() to set
        mris->vertices[idx].y = to;
    }
    void Vertex::set_z(float to) {
        mris->vertices[idx].z = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                           ) {}                   
    Surface::Surface ( MRIS* mris                ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src       ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace XYZPositionM


    namespace XYZPositionConsequencesM {
    Face::Face (                        ) {}                     
    Face::Face ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src       ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    float Face::area() const {
        return mris->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return mris->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return mris->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return mris->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return mris->faces[idx].gradNorm;
    }
    
    void Face::set_area(float to) {
        mris->faces[idx].area = to;
    }
    void Face::set_angle(angles_per_triangle_t to) {
        mris->faces[idx].angle = to;
    }
    void Face::set_orig_angle(angles_per_triangle_t to) {
        mris->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }
    void Face::set_norm(PDMATRIX to) {
        mris->faces[idx].norm = to;
    }
    void Face::set_gradNorm(A3PDMATRIX to) {
        mris->faces[idx].gradNorm = to;
    }


    Vertex::Vertex (                          ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx   ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src       ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return mris->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return mris->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return mris->vertices[idx].nz;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    float Vertex::area() const {
        return mris->vertices[idx].area;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        mris->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        mris->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        mris->vertices[idx].nz = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_area(float to) {
        mris->vertices[idx].area = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                           ) {}                   
    Surface::Surface ( MRIS* mris                ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src       ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return mris->strips[i];
    }
    float Surface::xctr() const {
        return mris->xctr;
    }
    float Surface::yctr() const {
        return mris->yctr;
    }
    float Surface::zctr() const {
        return mris->zctr;
    }
    float Surface::xlo() const {
        return mris->xlo;
    }
    float Surface::ylo() const {
        return mris->ylo;
    }
    float Surface::zlo() const {
        return mris->zlo;
    }
    float Surface::xhi() const {
        return mris->xhi;
    }
    float Surface::yhi() const {
        return mris->yhi;
    }
    float Surface::zhi() const {
        return mris->zhi;
    }
    float Surface::x0() const {  //  center of spherical expansion
        return mris->x0;
    }
    float Surface::y0() const {
        return mris->y0;
    }
    float Surface::z0() const {
        return mris->z0;
    }
    Vertex Surface::v_temporal_pole() const {
        return Vertex(mris,mris->v_temporal_pole - mris->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(mris,mris->v_frontal_pole - mris->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(mris,mris->v_occipital_pole - mris->vertices);
    }
    float Surface::max_curv() const {
        return mris->max_curv;
    }
    float Surface::min_curv() const {
        return mris->min_curv;
    }
    float Surface::total_area() const {
        return mris->total_area;
    }
    double Surface::avg_vertex_area() const {
        return mris->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return mris->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return mris->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return mris->orig_area;
    }
    float Surface::neg_area() const {
        return mris->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return mris->neg_orig_area;
    }
    int Surface::zeros() const {
        return mris->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return mris->hemisphere;
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return mris->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return mris->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return mris->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return mris->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return mris->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::nlabels() const {
        return mris->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return mris->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        mris->strips[i] = to;
    }
    void Surface::set_xctr(float to) {
        mris->xctr = to;
    }
    void Surface::set_yctr(float to) {
        mris->yctr = to;
    }
    void Surface::set_zctr(float to) {
        mris->zctr = to;
    }
    void Surface::set_xlo(float to) {
        mris->xlo = to;
    }
    void Surface::set_ylo(float to) {
        mris->ylo = to;
    }
    void Surface::set_zlo(float to) {
        mris->zlo = to;
    }
    void Surface::set_xhi(float to) {
        mris->xhi = to;
    }
    void Surface::set_yhi(float to) {
        mris->yhi = to;
    }
    void Surface::set_zhi(float to) {
        mris->zhi = to;
    }
    void Surface::set_x0(float to) {  //  center of spherical expansion
        mris->x0 = to;
    }
    void Surface::set_y0(float to) {
        mris->y0 = to;
    }
    void Surface::set_z0(float to) {
        mris->z0 = to;
    }
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_temporal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_frontal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_occipital_pole = mris->vertices + to.idx;
    }
    void Surface::set_max_curv(float to) {
        mris->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        mris->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        mris->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        mris->avg_vertex_area = to;
    }
    void Surface::set_avg_vertex_dist(double to) {  //  set by MRIScomputeAvgInterVertexDist
        mris->avg_vertex_dist = to;
    }
    void Surface::set_std_vertex_dist(double to) {
        mris->std_vertex_dist = to;
    }
    void Surface::set_orig_area(float to) {
        mris->orig_area = to;
    }
    void Surface::set_neg_area(float to) {
        mris->neg_area = to;
    }
    void Surface::set_neg_orig_area(float to) {  //  amount of original surface in folds
        mris->neg_orig_area = to;
    }
    void Surface::set_zeros(int to) {
        mris->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        mris->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        mris->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        mris->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        mris->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        mris->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        mris->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        mris->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        mris->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace XYZPositionConsequencesM


    namespace DistortM {
    Face::Face (                        ) {}                     
    Face::Face ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src       ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    float Face::area() const {
        return mris->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return mris->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return mris->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return mris->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return mris->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        mris->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                          ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx   ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src       ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return mris->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return mris->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return mris->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return mris->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return mris->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return mris->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return mris->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return mris->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return mris->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return mris->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return mris->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return mris->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return mris->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return mris->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return mris->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return mris->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return mris->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return mris->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return mris->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return mris->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return mris->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return mris->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return mris->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return mris->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return mris->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return mris->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return mris->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz;
    }
    float Vertex::tx2() const {
        return mris->vertices[idx].tx2;
    }
    float Vertex::ty2() const {
        return mris->vertices[idx].ty2;
    }
    float Vertex::tz2() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz2;
    }
    float Vertex::targx() const {
        return mris->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return mris->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return mris->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return mris->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return mris->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return mris->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return mris->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return mris->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return mris->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return mris->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return mris->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return mris->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return mris->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return mris->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return mris->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return mris->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return mris->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return mris->vertices[idx].fz;
    }
    int Vertex::px() const {
        return mris->vertices[idx].px;
    }
    int Vertex::qx() const {
        return mris->vertices[idx].qx;
    }
    int Vertex::py() const {
        return mris->vertices[idx].py;
    }
    int Vertex::qy() const {
        return mris->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return mris->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return mris->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return mris->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return mris->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return mris->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return mris->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return mris->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return mris->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return mris->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return mris->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return mris->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return mris->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return mris->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return mris->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return mris->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return mris->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return mris->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return mris->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return mris->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return mris->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return mris->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return mris->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return mris->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return mris->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return mris->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return mris->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return mris->vertices[idx].phi;
    }
    float Vertex::area() const {
        return mris->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return mris->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return mris->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return mris->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return mris->vertices[idx].H;
    }
    float Vertex::k1() const {
        return mris->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return mris->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return mris->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return mris->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return mris->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return mris->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return mris->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return mris->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return mris->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return mris->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return mris->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return mris->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return mris->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        mris->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        mris->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        mris->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        mris->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        mris->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        mris->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        mris->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        mris->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        mris->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        mris->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        mris->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        mris->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        mris->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        mris->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        mris->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        mris->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        mris->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        mris->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        mris->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        mris->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        mris->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        mris->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        mris->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        mris->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        mris->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        mris->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        mris->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz = to;
    }
    void Vertex::set_tx2(float to) {
        mris->vertices[idx].tx2 = to;
    }
    void Vertex::set_ty2(float to) {
        mris->vertices[idx].ty2 = to;
    }
    void Vertex::set_tz2(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz2 = to;
    }
    void Vertex::set_targx(float to) {
        mris->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        mris->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        mris->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        mris->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        mris->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        mris->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        mris->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        mris->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        mris->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        mris->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        mris->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        mris->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        mris->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        mris->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        mris->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        mris->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        mris->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        mris->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        mris->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        mris->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        mris->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        mris->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        mris->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        mris->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        mris->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        mris->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        mris->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        mris->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        mris->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        mris->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        mris->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        mris->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        mris->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        mris->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        mris->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        mris->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        mris->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        mris->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        mris->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        mris->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        mris->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        mris->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        mris->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        mris->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        mris->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        mris->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        mris->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        mris->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        mris->vertices[idx].phi = to;
    }
    void Vertex::set_origarea(float to) {
        mris->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        mris->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        mris->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        mris->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        mris->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        mris->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        mris->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        mris->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        mris->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        mris->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        mris->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        mris->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        mris->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        mris->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        mris->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        mris->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        mris->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                           ) {}                   
    Surface::Surface ( MRIS* mris                ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src       ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return mris->strips[i];
    }
    float Surface::xctr() const {
        return mris->xctr;
    }
    float Surface::yctr() const {
        return mris->yctr;
    }
    float Surface::zctr() const {
        return mris->zctr;
    }
    float Surface::xlo() const {
        return mris->xlo;
    }
    float Surface::ylo() const {
        return mris->ylo;
    }
    float Surface::zlo() const {
        return mris->zlo;
    }
    float Surface::xhi() const {
        return mris->xhi;
    }
    float Surface::yhi() const {
        return mris->yhi;
    }
    float Surface::zhi() const {
        return mris->zhi;
    }
    float Surface::x0() const {  //  center of spherical expansion
        return mris->x0;
    }
    float Surface::y0() const {
        return mris->y0;
    }
    float Surface::z0() const {
        return mris->z0;
    }
    Vertex Surface::v_temporal_pole() const {
        return Vertex(mris,mris->v_temporal_pole - mris->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(mris,mris->v_frontal_pole - mris->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(mris,mris->v_occipital_pole - mris->vertices);
    }
    float Surface::max_curv() const {
        return mris->max_curv;
    }
    float Surface::min_curv() const {
        return mris->min_curv;
    }
    float Surface::total_area() const {
        return mris->total_area;
    }
    double Surface::avg_vertex_area() const {
        return mris->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return mris->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return mris->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return mris->orig_area;
    }
    float Surface::neg_area() const {
        return mris->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return mris->neg_orig_area;
    }
    int Surface::zeros() const {
        return mris->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return mris->hemisphere;
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return mris->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return mris->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return mris->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return mris->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return mris->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::nlabels() const {
        return mris->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return mris->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        mris->strips[i] = to;
    }
    void Surface::set_xctr(float to) {
        mris->xctr = to;
    }
    void Surface::set_yctr(float to) {
        mris->yctr = to;
    }
    void Surface::set_zctr(float to) {
        mris->zctr = to;
    }
    void Surface::set_xlo(float to) {
        mris->xlo = to;
    }
    void Surface::set_ylo(float to) {
        mris->ylo = to;
    }
    void Surface::set_zlo(float to) {
        mris->zlo = to;
    }
    void Surface::set_xhi(float to) {
        mris->xhi = to;
    }
    void Surface::set_yhi(float to) {
        mris->yhi = to;
    }
    void Surface::set_zhi(float to) {
        mris->zhi = to;
    }
    void Surface::set_x0(float to) {  //  center of spherical expansion
        mris->x0 = to;
    }
    void Surface::set_y0(float to) {
        mris->y0 = to;
    }
    void Surface::set_z0(float to) {
        mris->z0 = to;
    }
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_temporal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_frontal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_occipital_pole = mris->vertices + to.idx;
    }
    void Surface::set_max_curv(float to) {
        mris->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        mris->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        mris->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        mris->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        mris->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        mris->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        mris->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        mris->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        mris->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        mris->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        mris->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        mris->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        mris->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace DistortM


    namespace AnalysisM {
    Face::Face (                        ) {}                     
    Face::Face ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src       ) : MRIS_Elt(src) {}     
    Face::Face ( AllM::Face const & src ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    float Face::area() const {
        return mris->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return mris->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return mris->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return mris->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return mris->faces[idx].gradNorm;
    }
    
    void Face::set_orig_angle(angles_per_triangle_t to) {
        mris->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }


    Vertex::Vertex (                          ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx   ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src       ) : MRIS_Elt(src) {}     
    Vertex::Vertex ( AllM::Vertex const & src ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return mris->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return mris->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return mris->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return mris->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return mris->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return mris->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return mris->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return mris->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return mris->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return mris->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return mris->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return mris->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return mris->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return mris->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return mris->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return mris->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return mris->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return mris->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return mris->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return mris->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return mris->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return mris->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return mris->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return mris->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return mris->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return mris->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return mris->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz;
    }
    float Vertex::tx2() const {
        return mris->vertices[idx].tx2;
    }
    float Vertex::ty2() const {
        return mris->vertices[idx].ty2;
    }
    float Vertex::tz2() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz2;
    }
    float Vertex::targx() const {
        return mris->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return mris->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return mris->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return mris->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return mris->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return mris->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return mris->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return mris->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return mris->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return mris->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return mris->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return mris->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return mris->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return mris->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return mris->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return mris->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return mris->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return mris->vertices[idx].fz;
    }
    int Vertex::px() const {
        return mris->vertices[idx].px;
    }
    int Vertex::qx() const {
        return mris->vertices[idx].qx;
    }
    int Vertex::py() const {
        return mris->vertices[idx].py;
    }
    int Vertex::qy() const {
        return mris->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return mris->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return mris->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return mris->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return mris->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return mris->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return mris->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return mris->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return mris->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return mris->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return mris->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return mris->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return mris->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return mris->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return mris->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return mris->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return mris->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return mris->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return mris->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return mris->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return mris->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return mris->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return mris->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return mris->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return mris->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return mris->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return mris->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return mris->vertices[idx].phi;
    }
    float Vertex::area() const {
        return mris->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return mris->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return mris->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return mris->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return mris->vertices[idx].H;
    }
    float Vertex::k1() const {
        return mris->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return mris->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return mris->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return mris->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return mris->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return mris->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return mris->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return mris->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return mris->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return mris->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return mris->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return mris->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return mris->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_nx(float to) {
        mris->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        mris->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        mris->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        mris->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        mris->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        mris->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        mris->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        mris->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        mris->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        mris->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        mris->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        mris->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        mris->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        mris->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        mris->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        mris->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        mris->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        mris->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        mris->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        mris->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        mris->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        mris->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        mris->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        mris->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        mris->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        mris->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        mris->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz = to;
    }
    void Vertex::set_tx2(float to) {
        mris->vertices[idx].tx2 = to;
    }
    void Vertex::set_ty2(float to) {
        mris->vertices[idx].ty2 = to;
    }
    void Vertex::set_tz2(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz2 = to;
    }
    void Vertex::set_targx(float to) {
        mris->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        mris->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        mris->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        mris->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        mris->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        mris->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        mris->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        mris->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        mris->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        mris->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        mris->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        mris->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        mris->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        mris->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        mris->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        mris->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        mris->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        mris->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        mris->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        mris->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        mris->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        mris->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        mris->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        mris->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        mris->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        mris->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        mris->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        mris->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        mris->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        mris->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        mris->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        mris->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        mris->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        mris->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        mris->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        mris->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        mris->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        mris->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        mris->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        mris->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        mris->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        mris->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        mris->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        mris->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        mris->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        mris->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        mris->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        mris->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        mris->vertices[idx].phi = to;
    }
    void Vertex::set_origarea(float to) {
        mris->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        mris->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        mris->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        mris->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        mris->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        mris->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        mris->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        mris->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        mris->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        mris->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        mris->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        mris->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        mris->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        mris->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        mris->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        mris->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        mris->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                           ) {}                   
    Surface::Surface ( MRIS* mris                ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src       ) : MRIS_Elt(src) {}   
    Surface::Surface ( AllM::Surface const & src ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return mris->strips[i];
    }
    float Surface::xctr() const {
        return mris->xctr;
    }
    float Surface::yctr() const {
        return mris->yctr;
    }
    float Surface::zctr() const {
        return mris->zctr;
    }
    float Surface::xlo() const {
        return mris->xlo;
    }
    float Surface::ylo() const {
        return mris->ylo;
    }
    float Surface::zlo() const {
        return mris->zlo;
    }
    float Surface::xhi() const {
        return mris->xhi;
    }
    float Surface::yhi() const {
        return mris->yhi;
    }
    float Surface::zhi() const {
        return mris->zhi;
    }
    float Surface::x0() const {  //  center of spherical expansion
        return mris->x0;
    }
    float Surface::y0() const {
        return mris->y0;
    }
    float Surface::z0() const {
        return mris->z0;
    }
    Vertex Surface::v_temporal_pole() const {
        return Vertex(mris,mris->v_temporal_pole - mris->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(mris,mris->v_frontal_pole - mris->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(mris,mris->v_occipital_pole - mris->vertices);
    }
    float Surface::max_curv() const {
        return mris->max_curv;
    }
    float Surface::min_curv() const {
        return mris->min_curv;
    }
    float Surface::total_area() const {
        return mris->total_area;
    }
    double Surface::avg_vertex_area() const {
        return mris->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return mris->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return mris->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return mris->orig_area;
    }
    float Surface::neg_area() const {
        return mris->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return mris->neg_orig_area;
    }
    int Surface::zeros() const {
        return mris->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return mris->hemisphere;
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return mris->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return mris->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return mris->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return mris->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return mris->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::nlabels() const {
        return mris->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return mris->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        mris->strips[i] = to;
    }
    void Surface::set_xctr(float to) {
        mris->xctr = to;
    }
    void Surface::set_yctr(float to) {
        mris->yctr = to;
    }
    void Surface::set_zctr(float to) {
        mris->zctr = to;
    }
    void Surface::set_xlo(float to) {
        mris->xlo = to;
    }
    void Surface::set_ylo(float to) {
        mris->ylo = to;
    }
    void Surface::set_zlo(float to) {
        mris->zlo = to;
    }
    void Surface::set_xhi(float to) {
        mris->xhi = to;
    }
    void Surface::set_yhi(float to) {
        mris->yhi = to;
    }
    void Surface::set_zhi(float to) {
        mris->zhi = to;
    }
    void Surface::set_x0(float to) {  //  center of spherical expansion
        mris->x0 = to;
    }
    void Surface::set_y0(float to) {
        mris->y0 = to;
    }
    void Surface::set_z0(float to) {
        mris->z0 = to;
    }
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_temporal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_frontal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_occipital_pole = mris->vertices + to.idx;
    }
    void Surface::set_max_curv(float to) {
        mris->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        mris->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        mris->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        mris->avg_vertex_area = to;
    }
    void Surface::set_zeros(int to) {
        mris->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        mris->hemisphere = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        mris->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        mris->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        mris->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        mris->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        mris->Ktotal = to;
    }
    void Surface::set_nlabels(int to) {
        mris->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        mris->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace AnalysisM


    namespace AllM {
    Face::Face (                        ) {}                     
    Face::Face ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Face::Face ( Face const & src       ) : MRIS_Elt(src) {}     

    vertices_per_face_t Face::v() const {
        return mris->faces[idx].v;
    }
    float Face::area() const {
        return mris->faces[idx].area;
    }
    angles_per_triangle_t Face::angle() const {
        return mris->faces[idx].angle;
    }
    angles_per_triangle_t Face::orig_angle() const {
        return mris->faces[idx].orig_angle;
    }
    char Face::ripflag() const {
        return mris->faces[idx].ripflag;
    }
    char Face::oripflag() const {
        return mris->faces[idx].oripflag;
    }
    int Face::marked() const {
        return mris->faces[idx].marked;
    }
    PDMATRIX Face::norm() const {
        return mris->faces[idx].norm;
    }
    A3PDMATRIX Face::gradNorm() const {
        return mris->faces[idx].gradNorm;
    }
    
    void Face::set_v(vertices_per_face_t to) {
        mris->faces[idx].v = to;
    }
    void Face::set_area(float to) {
        mris->faces[idx].area = to;
    }
    void Face::set_angle(angles_per_triangle_t to) {
        mris->faces[idx].angle = to;
    }
    void Face::set_orig_angle(angles_per_triangle_t to) {
        mris->faces[idx].orig_angle = to;
    }
    void Face::set_ripflag(char to) {
        mris->faces[idx].ripflag = to;
    }
    void Face::set_oripflag(char to) {
        mris->faces[idx].oripflag = to;
    }
    void Face::set_marked(int to) {
        mris->faces[idx].marked = to;
    }
    void Face::set_norm(PDMATRIX to) {
        mris->faces[idx].norm = to;
    }
    void Face::set_gradNorm(A3PDMATRIX to) {
        mris->faces[idx].gradNorm = to;
    }


    Vertex::Vertex (                        ) {}                     
    Vertex::Vertex ( MRIS* mris, size_t idx ) : MRIS_Elt(mris,idx) {}
    Vertex::Vertex ( Vertex const & src     ) : MRIS_Elt(src) {}     

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(mris,mris->vertices_topology[idx].f[i]);
    }
    size_t Vertex::n(size_t i) const {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        return size_t(mris->vertices_topology[idx].n[i]);
    }
    int Vertex::e(size_t i) const {  //  edge state for neighboring vertices                      
        return mris->vertices_topology[idx].e[i];
    }
    Vertex Vertex::v(size_t i) const {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        return Vertex(mris,mris->vertices_topology[idx].v[i]);
    }
    short Vertex::vnum() const {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        return mris->vertices_topology[idx].vnum;
    }
    short Vertex::v2num() const {  //  number of 1, or 2-hop neighbors                          
        return mris->vertices_topology[idx].v2num;
    }
    short Vertex::v3num() const {  //  number of 1,2,or 3-hop neighbors                         
        return mris->vertices_topology[idx].v3num;
    }
    short Vertex::vtotal() const {  //  total # of neighbors. copy of vnum.nsizeCur              
        return mris->vertices_topology[idx].vtotal;
    }
    short Vertex::nsizeMaxClock() const {  //  copy of mris->nsizeMaxClock when v#num                   
        return mris->vertices_topology[idx].nsizeMaxClock;
    }
    uchar Vertex::nsizeMax() const {  //  the max nsize that was used to fill in vnum etc          
        return mris->vertices_topology[idx].nsizeMax;
    }
    uchar Vertex::nsizeCur() const {  //  index of the current v#num in vtotal                     
        return mris->vertices_topology[idx].nsizeCur;
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return mris->vertices_topology[idx].num;
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return mris->vertices[idx].dist[i];
    }
    float Vertex::dist_orig(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on origxyz
        return mris->vertices[idx].dist_orig[i];
    }
    int Vertex::dist_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_capacity;
    }
    int Vertex::dist_orig_capacity() const {  //  -- should contain at least vtx_vtotal elements   
        return mris->vertices[idx].dist_orig_capacity;
    }
    float Vertex::x() const {  //  current coordinates	
        return mris->vertices[idx].x;
    }
    float Vertex::y() const {  //  use MRISsetXYZ() to set
        return mris->vertices[idx].y;
    }
    float Vertex::z() const {
        return mris->vertices[idx].z;
    }
    float Vertex::origx() const {  //  original coordinates, see also MRIS::origxyz_status
        return mris->vertices[idx].origx;
    }
    float Vertex::origy() const {  //  use MRISsetOriginalXYZ(, 
        return mris->vertices[idx].origy;
    }
    float Vertex::origz() const {  //  or MRISsetOriginalXYZfromXYZ to set
        return mris->vertices[idx].origz;
    }
    float Vertex::nx() const {
        return mris->vertices[idx].nx;
    }
    float Vertex::ny() const {
        return mris->vertices[idx].ny;
    }
    float Vertex::nz() const {  //  curr normal
        return mris->vertices[idx].nz;
    }
    float Vertex::pnx() const {
        return mris->vertices[idx].pnx;
    }
    float Vertex::pny() const {
        return mris->vertices[idx].pny;
    }
    float Vertex::pnz() const {  //  pial normal
        return mris->vertices[idx].pnz;
    }
    float Vertex::wnx() const {
        return mris->vertices[idx].wnx;
    }
    float Vertex::wny() const {
        return mris->vertices[idx].wny;
    }
    float Vertex::wnz() const {  //  white normal
        return mris->vertices[idx].wnz;
    }
    float Vertex::onx() const {
        return mris->vertices[idx].onx;
    }
    float Vertex::ony() const {
        return mris->vertices[idx].ony;
    }
    float Vertex::onz() const {  //  original normal
        return mris->vertices[idx].onz;
    }
    float Vertex::dx() const {
        return mris->vertices[idx].dx;
    }
    float Vertex::dy() const {
        return mris->vertices[idx].dy;
    }
    float Vertex::dz() const {  //  current change in position
        return mris->vertices[idx].dz;
    }
    float Vertex::odx() const {
        return mris->vertices[idx].odx;
    }
    float Vertex::ody() const {
        return mris->vertices[idx].ody;
    }
    float Vertex::odz() const {  //  last change of position (for momentum, 
        return mris->vertices[idx].odz;
    }
    float Vertex::tdx() const {
        return mris->vertices[idx].tdx;
    }
    float Vertex::tdy() const {
        return mris->vertices[idx].tdy;
    }
    float Vertex::tdz() const {  //  temporary storage for averaging gradient
        return mris->vertices[idx].tdz;
    }
    float Vertex::curv() const {  //  curr curvature
        return mris->vertices[idx].curv;
    }
    float Vertex::curvbak() const {
        return mris->vertices[idx].curvbak;
    }
    float Vertex::val() const {  //  scalar data value (file: rh.val, sig2-rh.w)
        return mris->vertices[idx].val;
    }
    float Vertex::imag_val() const {  //  imaginary part of complex data value
        return mris->vertices[idx].imag_val;
    }
    float Vertex::cx() const {
        return mris->vertices[idx].cx;
    }
    float Vertex::cy() const {
        return mris->vertices[idx].cy;
    }
    float Vertex::cz() const {  //  coordinates in canonical coordinate system
        return mris->vertices[idx].cz;
    }
    float Vertex::tx() const {
        return mris->vertices[idx].tx;
    }
    float Vertex::ty() const {
        return mris->vertices[idx].ty;
    }
    float Vertex::tz() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz;
    }
    float Vertex::tx2() const {
        return mris->vertices[idx].tx2;
    }
    float Vertex::ty2() const {
        return mris->vertices[idx].ty2;
    }
    float Vertex::tz2() const {  //  tmp coordinate storage
        return mris->vertices[idx].tz2;
    }
    float Vertex::targx() const {
        return mris->vertices[idx].targx;
    }
    float Vertex::targy() const {
        return mris->vertices[idx].targy;
    }
    float Vertex::targz() const {  //  target coordinates
        return mris->vertices[idx].targz;
    }
    float Vertex::pialx() const {
        return mris->vertices[idx].pialx;
    }
    float Vertex::pialy() const {
        return mris->vertices[idx].pialy;
    }
    float Vertex::pialz() const {  //  pial surface coordinates
        return mris->vertices[idx].pialz;
    }
    float Vertex::whitex() const {
        return mris->vertices[idx].whitex;
    }
    float Vertex::whitey() const {
        return mris->vertices[idx].whitey;
    }
    float Vertex::whitez() const {  //  white surface coordinates
        return mris->vertices[idx].whitez;
    }
    float Vertex::l4x() const {
        return mris->vertices[idx].l4x;
    }
    float Vertex::l4y() const {
        return mris->vertices[idx].l4y;
    }
    float Vertex::l4z() const {  //  layerIV surface coordinates
        return mris->vertices[idx].l4z;
    }
    float Vertex::infx() const {
        return mris->vertices[idx].infx;
    }
    float Vertex::infy() const {
        return mris->vertices[idx].infy;
    }
    float Vertex::infz() const {  //  inflated coordinates
        return mris->vertices[idx].infz;
    }
    float Vertex::fx() const {
        return mris->vertices[idx].fx;
    }
    float Vertex::fy() const {
        return mris->vertices[idx].fy;
    }
    float Vertex::fz() const {  //  flattened coordinates
        return mris->vertices[idx].fz;
    }
    int Vertex::px() const {
        return mris->vertices[idx].px;
    }
    int Vertex::qx() const {
        return mris->vertices[idx].qx;
    }
    int Vertex::py() const {
        return mris->vertices[idx].py;
    }
    int Vertex::qy() const {
        return mris->vertices[idx].qy;
    }
    int Vertex::pz() const {
        return mris->vertices[idx].pz;
    }
    int Vertex::qz() const {  //  rational coordinates for exact calculations
        return mris->vertices[idx].qz;
    }
    float Vertex::e1x() const {
        return mris->vertices[idx].e1x;
    }
    float Vertex::e1y() const {
        return mris->vertices[idx].e1y;
    }
    float Vertex::e1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].e1z;
    }
    float Vertex::e2x() const {
        return mris->vertices[idx].e2x;
    }
    float Vertex::e2y() const {
        return mris->vertices[idx].e2y;
    }
    float Vertex::e2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].e2z;
    }
    float Vertex::pe1x() const {
        return mris->vertices[idx].pe1x;
    }
    float Vertex::pe1y() const {
        return mris->vertices[idx].pe1y;
    }
    float Vertex::pe1z() const {  //  1st basis vector for the local tangent plane
        return mris->vertices[idx].pe1z;
    }
    float Vertex::pe2x() const {
        return mris->vertices[idx].pe2x;
    }
    float Vertex::pe2y() const {
        return mris->vertices[idx].pe2y;
    }
    float Vertex::pe2z() const {  //  2nd basis vector for the local tangent plane
        return mris->vertices[idx].pe2z;
    }
    float Vertex::nc() const {  //  curr length normal comp 
        return mris->vertices[idx].nc;
    }
    float Vertex::val2() const {  //  complex comp data value (file: sig3-rh.w) 
        return mris->vertices[idx].val2;
    }
    float Vertex::valbak() const {  //  scalar data stack 
        return mris->vertices[idx].valbak;
    }
    float Vertex::val2bak() const {  //  complex comp data stack 
        return mris->vertices[idx].val2bak;
    }
    float Vertex::stat() const {  //  statistic 
        return mris->vertices[idx].stat;
    }
    int Vertex::undefval() const {  //  [previously dist=0] 
        return mris->vertices[idx].undefval;
    }
    int Vertex::old_undefval() const {  //  for smooth_val_sparse 
        return mris->vertices[idx].old_undefval;
    }
    int Vertex::fixedval() const {  //  [previously val=0] 
        return mris->vertices[idx].fixedval;
    }
    float Vertex::fieldsign() const {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        return mris->vertices[idx].fieldsign;
    }
    float Vertex::fsmask() const {  //  significance mask (file: rh.fm) 
        return mris->vertices[idx].fsmask;
    }
    float Vertex::d() const {  //  for distance calculations 
        return mris->vertices[idx].d;
    }
    int Vertex::annotation() const {  //  area label (defunct--now from label file name!) 
        return mris->vertices[idx].annotation;
    }
    char Vertex::oripflag() const {
        return mris->vertices[idx].oripflag;
    }
    char Vertex::origripflag() const {  //  cuts flags 
        return mris->vertices[idx].origripflag;
    }
    p_void Vertex::vp() const {  //  to store user's information 
        return mris->vertices[idx].vp;
    }
    float Vertex::theta() const {
        return mris->vertices[idx].theta;
    }
    float Vertex::phi() const {  //  parameterization 
        return mris->vertices[idx].phi;
    }
    float Vertex::area() const {
        return mris->vertices[idx].area;
    }
    float Vertex::origarea() const {
        return mris->vertices[idx].origarea;
    }
    float Vertex::group_avg_area() const {
        return mris->vertices[idx].group_avg_area;
    }
    float Vertex::K() const {  //  Gaussian curvature 
        return mris->vertices[idx].K;
    }
    float Vertex::H() const {  //  mean curvature 
        return mris->vertices[idx].H;
    }
    float Vertex::k1() const {
        return mris->vertices[idx].k1;
    }
    float Vertex::k2() const {  //  the principal curvatures 
        return mris->vertices[idx].k2;
    }
    float Vertex::mean() const {
        return mris->vertices[idx].mean;
    }
    float Vertex::mean_imag() const {  //  imaginary part of complex statistic 
        return mris->vertices[idx].mean_imag;
    }
    float Vertex::std_error() const {
        return mris->vertices[idx].std_error;
    }
    uint Vertex::flags() const {
        return mris->vertices[idx].flags;
    }
    int Vertex::fno() const {  //  face that this vertex is in 
        return mris->vertices[idx].fno;
    }
    int Vertex::cropped() const {
        return mris->vertices[idx].cropped;
    }
    short Vertex::marked() const {  //  for a variety of uses 
        return mris->vertices[idx].marked;
    }
    short Vertex::marked2() const {
        return mris->vertices[idx].marked2;
    }
    short Vertex::marked3() const {
        return mris->vertices[idx].marked3;
    }
    char Vertex::neg() const {  //  1 if the normal vector is inverted 
        return mris->vertices[idx].neg;
    }
    char Vertex::border() const {  //  flag 
        return mris->vertices[idx].border;
    }
    char Vertex::ripflag() const {  //  vertex no longer exists - placed last to load the next vertex into cache
        return mris->vertices[idx].ripflag;
    }
    
    void Vertex::set_f(size_t i, Face to) {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].f[i] = to.idx;
    }
    void Vertex::set_n(size_t i, size_t to) {  // size() is num.    array[v->num] the face.v[*] index for this vertex        
        mris->vertices_topology[idx].n[i] = uchar(to);
    }
    void Vertex::set_e(size_t i, int to) {  //  edge state for neighboring vertices                      
        mris->vertices_topology[idx].e[i] = to;
    }
    void Vertex::set_v(size_t i, Vertex to) {  // size() is vtotal.    array[v->vtotal or more] of vno, head sorted by hops     
        cheapAssert(mris == to.mris); mris->vertices_topology[idx].v[i] = to.idx;
    }
    void Vertex::set_vnum(short to) {  //  number of 1-hop neighbors    should use [p]VERTEXvnum(i, 
        mris->vertices_topology[idx].vnum = to;
    }
    void Vertex::set_v2num(short to) {  //  number of 1, or 2-hop neighbors                          
        mris->vertices_topology[idx].v2num = to;
    }
    void Vertex::set_v3num(short to) {  //  number of 1,2,or 3-hop neighbors                         
        mris->vertices_topology[idx].v3num = to;
    }
    void Vertex::set_vtotal(short to) {  //  total # of neighbors. copy of vnum.nsizeCur              
        mris->vertices_topology[idx].vtotal = to;
    }
    void Vertex::set_nsizeMaxClock(short to) {  //  copy of mris->nsizeMaxClock when v#num                   
        mris->vertices_topology[idx].nsizeMaxClock = to;
    }
    void Vertex::set_nsizeMax(uchar to) {  //  the max nsize that was used to fill in vnum etc          
        mris->vertices_topology[idx].nsizeMax = to;
    }
    void Vertex::set_nsizeCur(uchar to) {  //  index of the current v#num in vtotal                     
        mris->vertices_topology[idx].nsizeCur = to;
    }
    void Vertex::set_num(uchar to) {  //  number of neighboring faces                              
        mris->vertices_topology[idx].num = to;
    }
    void Vertex::set_x(float to) {  //  current coordinates	
        mris->vertices[idx].x = to;
    }
    void Vertex::set_y(float to) {  //  use MRISsetXYZ() to set
        mris->vertices[idx].y = to;
    }
    void Vertex::set_z(float to) {
        mris->vertices[idx].z = to;
    }
    void Vertex::set_nx(float to) {
        mris->vertices[idx].nx = to;
    }
    void Vertex::set_ny(float to) {
        mris->vertices[idx].ny = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        mris->vertices[idx].nz = to;
    }
    void Vertex::set_pnx(float to) {
        mris->vertices[idx].pnx = to;
    }
    void Vertex::set_pny(float to) {
        mris->vertices[idx].pny = to;
    }
    void Vertex::set_pnz(float to) {  //  pial normal
        mris->vertices[idx].pnz = to;
    }
    void Vertex::set_wnx(float to) {
        mris->vertices[idx].wnx = to;
    }
    void Vertex::set_wny(float to) {
        mris->vertices[idx].wny = to;
    }
    void Vertex::set_wnz(float to) {  //  white normal
        mris->vertices[idx].wnz = to;
    }
    void Vertex::set_onx(float to) {
        mris->vertices[idx].onx = to;
    }
    void Vertex::set_ony(float to) {
        mris->vertices[idx].ony = to;
    }
    void Vertex::set_onz(float to) {  //  original normal
        mris->vertices[idx].onz = to;
    }
    void Vertex::set_dx(float to) {
        mris->vertices[idx].dx = to;
    }
    void Vertex::set_dy(float to) {
        mris->vertices[idx].dy = to;
    }
    void Vertex::set_dz(float to) {  //  current change in position
        mris->vertices[idx].dz = to;
    }
    void Vertex::set_odx(float to) {
        mris->vertices[idx].odx = to;
    }
    void Vertex::set_ody(float to) {
        mris->vertices[idx].ody = to;
    }
    void Vertex::set_odz(float to) {  //  last change of position (for momentum, 
        mris->vertices[idx].odz = to;
    }
    void Vertex::set_tdx(float to) {
        mris->vertices[idx].tdx = to;
    }
    void Vertex::set_tdy(float to) {
        mris->vertices[idx].tdy = to;
    }
    void Vertex::set_tdz(float to) {  //  temporary storage for averaging gradient
        mris->vertices[idx].tdz = to;
    }
    void Vertex::set_curv(float to) {  //  curr curvature
        mris->vertices[idx].curv = to;
    }
    void Vertex::set_curvbak(float to) {
        mris->vertices[idx].curvbak = to;
    }
    void Vertex::set_val(float to) {  //  scalar data value (file: rh.val, sig2-rh.w)
        mris->vertices[idx].val = to;
    }
    void Vertex::set_imag_val(float to) {  //  imaginary part of complex data value
        mris->vertices[idx].imag_val = to;
    }
    void Vertex::set_cx(float to) {
        mris->vertices[idx].cx = to;
    }
    void Vertex::set_cy(float to) {
        mris->vertices[idx].cy = to;
    }
    void Vertex::set_cz(float to) {  //  coordinates in canonical coordinate system
        mris->vertices[idx].cz = to;
    }
    void Vertex::set_tx(float to) {
        mris->vertices[idx].tx = to;
    }
    void Vertex::set_ty(float to) {
        mris->vertices[idx].ty = to;
    }
    void Vertex::set_tz(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz = to;
    }
    void Vertex::set_tx2(float to) {
        mris->vertices[idx].tx2 = to;
    }
    void Vertex::set_ty2(float to) {
        mris->vertices[idx].ty2 = to;
    }
    void Vertex::set_tz2(float to) {  //  tmp coordinate storage
        mris->vertices[idx].tz2 = to;
    }
    void Vertex::set_targx(float to) {
        mris->vertices[idx].targx = to;
    }
    void Vertex::set_targy(float to) {
        mris->vertices[idx].targy = to;
    }
    void Vertex::set_targz(float to) {  //  target coordinates
        mris->vertices[idx].targz = to;
    }
    void Vertex::set_pialx(float to) {
        mris->vertices[idx].pialx = to;
    }
    void Vertex::set_pialy(float to) {
        mris->vertices[idx].pialy = to;
    }
    void Vertex::set_pialz(float to) {  //  pial surface coordinates
        mris->vertices[idx].pialz = to;
    }
    void Vertex::set_whitex(float to) {
        mris->vertices[idx].whitex = to;
    }
    void Vertex::set_whitey(float to) {
        mris->vertices[idx].whitey = to;
    }
    void Vertex::set_whitez(float to) {  //  white surface coordinates
        mris->vertices[idx].whitez = to;
    }
    void Vertex::set_l4x(float to) {
        mris->vertices[idx].l4x = to;
    }
    void Vertex::set_l4y(float to) {
        mris->vertices[idx].l4y = to;
    }
    void Vertex::set_l4z(float to) {  //  layerIV surface coordinates
        mris->vertices[idx].l4z = to;
    }
    void Vertex::set_infx(float to) {
        mris->vertices[idx].infx = to;
    }
    void Vertex::set_infy(float to) {
        mris->vertices[idx].infy = to;
    }
    void Vertex::set_infz(float to) {  //  inflated coordinates
        mris->vertices[idx].infz = to;
    }
    void Vertex::set_fx(float to) {
        mris->vertices[idx].fx = to;
    }
    void Vertex::set_fy(float to) {
        mris->vertices[idx].fy = to;
    }
    void Vertex::set_fz(float to) {  //  flattened coordinates
        mris->vertices[idx].fz = to;
    }
    void Vertex::set_px(int to) {
        mris->vertices[idx].px = to;
    }
    void Vertex::set_qx(int to) {
        mris->vertices[idx].qx = to;
    }
    void Vertex::set_py(int to) {
        mris->vertices[idx].py = to;
    }
    void Vertex::set_qy(int to) {
        mris->vertices[idx].qy = to;
    }
    void Vertex::set_pz(int to) {
        mris->vertices[idx].pz = to;
    }
    void Vertex::set_qz(int to) {  //  rational coordinates for exact calculations
        mris->vertices[idx].qz = to;
    }
    void Vertex::set_e1x(float to) {
        mris->vertices[idx].e1x = to;
    }
    void Vertex::set_e1y(float to) {
        mris->vertices[idx].e1y = to;
    }
    void Vertex::set_e1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].e1z = to;
    }
    void Vertex::set_e2x(float to) {
        mris->vertices[idx].e2x = to;
    }
    void Vertex::set_e2y(float to) {
        mris->vertices[idx].e2y = to;
    }
    void Vertex::set_e2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].e2z = to;
    }
    void Vertex::set_pe1x(float to) {
        mris->vertices[idx].pe1x = to;
    }
    void Vertex::set_pe1y(float to) {
        mris->vertices[idx].pe1y = to;
    }
    void Vertex::set_pe1z(float to) {  //  1st basis vector for the local tangent plane
        mris->vertices[idx].pe1z = to;
    }
    void Vertex::set_pe2x(float to) {
        mris->vertices[idx].pe2x = to;
    }
    void Vertex::set_pe2y(float to) {
        mris->vertices[idx].pe2y = to;
    }
    void Vertex::set_pe2z(float to) {  //  2nd basis vector for the local tangent plane
        mris->vertices[idx].pe2z = to;
    }
    void Vertex::set_nc(float to) {  //  curr length normal comp 
        mris->vertices[idx].nc = to;
    }
    void Vertex::set_val2(float to) {  //  complex comp data value (file: sig3-rh.w) 
        mris->vertices[idx].val2 = to;
    }
    void Vertex::set_valbak(float to) {  //  scalar data stack 
        mris->vertices[idx].valbak = to;
    }
    void Vertex::set_val2bak(float to) {  //  complex comp data stack 
        mris->vertices[idx].val2bak = to;
    }
    void Vertex::set_stat(float to) {  //  statistic 
        mris->vertices[idx].stat = to;
    }
    void Vertex::set_undefval(int to) {  //  [previously dist=0] 
        mris->vertices[idx].undefval = to;
    }
    void Vertex::set_old_undefval(int to) {  //  for smooth_val_sparse 
        mris->vertices[idx].old_undefval = to;
    }
    void Vertex::set_fixedval(int to) {  //  [previously val=0] 
        mris->vertices[idx].fixedval = to;
    }
    void Vertex::set_fieldsign(float to) {  //  fieldsign--final: -1, "0", "1" (file: rh.fs) 
        mris->vertices[idx].fieldsign = to;
    }
    void Vertex::set_fsmask(float to) {  //  significance mask (file: rh.fm) 
        mris->vertices[idx].fsmask = to;
    }
    void Vertex::set_d(float to) {  //  for distance calculations 
        mris->vertices[idx].d = to;
    }
    void Vertex::set_annotation(int to) {  //  area label (defunct--now from label file name!) 
        mris->vertices[idx].annotation = to;
    }
    void Vertex::set_oripflag(char to) {
        mris->vertices[idx].oripflag = to;
    }
    void Vertex::set_origripflag(char to) {  //  cuts flags 
        mris->vertices[idx].origripflag = to;
    }
    void Vertex::set_vp(p_void to) {  //  to store user's information 
        mris->vertices[idx].vp = to;
    }
    void Vertex::set_theta(float to) {
        mris->vertices[idx].theta = to;
    }
    void Vertex::set_phi(float to) {  //  parameterization 
        mris->vertices[idx].phi = to;
    }
    void Vertex::set_area(float to) {
        mris->vertices[idx].area = to;
    }
    void Vertex::set_origarea(float to) {
        mris->vertices[idx].origarea = to;
    }
    void Vertex::set_group_avg_area(float to) {
        mris->vertices[idx].group_avg_area = to;
    }
    void Vertex::set_K(float to) {  //  Gaussian curvature 
        mris->vertices[idx].K = to;
    }
    void Vertex::set_H(float to) {  //  mean curvature 
        mris->vertices[idx].H = to;
    }
    void Vertex::set_k1(float to) {
        mris->vertices[idx].k1 = to;
    }
    void Vertex::set_k2(float to) {  //  the principal curvatures 
        mris->vertices[idx].k2 = to;
    }
    void Vertex::set_mean(float to) {
        mris->vertices[idx].mean = to;
    }
    void Vertex::set_mean_imag(float to) {  //  imaginary part of complex statistic 
        mris->vertices[idx].mean_imag = to;
    }
    void Vertex::set_std_error(float to) {
        mris->vertices[idx].std_error = to;
    }
    void Vertex::set_flags(uint to) {
        mris->vertices[idx].flags = to;
    }
    void Vertex::set_fno(int to) {  //  face that this vertex is in 
        mris->vertices[idx].fno = to;
    }
    void Vertex::set_cropped(int to) {
        mris->vertices[idx].cropped = to;
    }
    void Vertex::set_marked(short to) {  //  for a variety of uses 
        mris->vertices[idx].marked = to;
    }
    void Vertex::set_marked2(short to) {
        mris->vertices[idx].marked2 = to;
    }
    void Vertex::set_marked3(short to) {
        mris->vertices[idx].marked3 = to;
    }
    void Vertex::set_neg(char to) {  //  1 if the normal vector is inverted 
        mris->vertices[idx].neg = to;
    }
    void Vertex::set_border(char to) {  //  flag 
        mris->vertices[idx].border = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        mris->vertices[idx].ripflag = to;
    }


    Surface::Surface (                     ) {}                   
    Surface::Surface ( MRIS* mris          ) : MRIS_Elt(mris,0) {}
    Surface::Surface ( Surface const & src ) : MRIS_Elt(src) {}   

    int Surface::nverticesFrozen() const {  //  # of vertices on surface is frozen
        return mris->nverticesFrozen;
    }
    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return mris->nfaces;
    }
    bool Surface::faceAttachmentDeferred() const {  //  defer connecting faces to vertices for performance reasons
        return mris->faceAttachmentDeferred;
    }
    int Surface::nedges() const {  //  # of edges on surface
        return mris->nedges;
    }
    int Surface::nstrips() const {
        return mris->nstrips;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(mris,idx);;
    }
    p_p_void Surface::dist_storage() const {  //  the malloced/realloced vertex dist fields, so those fields can be quickly nulled and restored
        return mris->dist_storage;
    }
    p_p_void Surface::dist_orig_storage() const {  //  the malloced/realloced vertex dist_orig fields, so those fields can be quickly nulled and restored
        return mris->dist_orig_storage;
    }
    int Surface::tempsAssigned() const {  //  State of various temp fields that can be borrowed if not already in use
        return mris->tempsAssigned;
    }
    Face Surface::faces(size_t i) const {
        return Face(mris,idx);;
    }
    MRI_EDGE Surface::edges(size_t i) const {
        return mris->edges[i];
    }
    FaceNormCacheEntry Surface::faceNormCacheEntries(size_t i) const {
        return mris->faceNormCacheEntries[i];
    }
    FaceNormDeferredEntry Surface::faceNormDeferredEntries(size_t i) const {
        return mris->faceNormDeferredEntries[i];
    }
    STRIP Surface::strips(size_t i) const {
        return mris->strips[i];
    }
    float Surface::xctr() const {
        return mris->xctr;
    }
    float Surface::yctr() const {
        return mris->yctr;
    }
    float Surface::zctr() const {
        return mris->zctr;
    }
    float Surface::xlo() const {
        return mris->xlo;
    }
    float Surface::ylo() const {
        return mris->ylo;
    }
    float Surface::zlo() const {
        return mris->zlo;
    }
    float Surface::xhi() const {
        return mris->xhi;
    }
    float Surface::yhi() const {
        return mris->yhi;
    }
    float Surface::zhi() const {
        return mris->zhi;
    }
    float Surface::x0() const {  //  center of spherical expansion
        return mris->x0;
    }
    float Surface::y0() const {
        return mris->y0;
    }
    float Surface::z0() const {
        return mris->z0;
    }
    Vertex Surface::v_temporal_pole() const {
        return Vertex(mris,mris->v_temporal_pole - mris->vertices);
    }
    Vertex Surface::v_frontal_pole() const {
        return Vertex(mris,mris->v_frontal_pole - mris->vertices);
    }
    Vertex Surface::v_occipital_pole() const {
        return Vertex(mris,mris->v_occipital_pole - mris->vertices);
    }
    float Surface::max_curv() const {
        return mris->max_curv;
    }
    float Surface::min_curv() const {
        return mris->min_curv;
    }
    float Surface::total_area() const {
        return mris->total_area;
    }
    double Surface::avg_vertex_area() const {
        return mris->avg_vertex_area;
    }
    double Surface::avg_vertex_dist() const {  //  set by MRIScomputeAvgInterVertexDist
        return mris->avg_vertex_dist;
    }
    double Surface::std_vertex_dist() const {
        return mris->std_vertex_dist;
    }
    float Surface::orig_area() const {
        return mris->orig_area;
    }
    float Surface::neg_area() const {
        return mris->neg_area;
    }
    float Surface::neg_orig_area() const {  //  amount of original surface in folds
        return mris->neg_orig_area;
    }
    int Surface::zeros() const {
        return mris->zeros;
    }
    int Surface::hemisphere() const {  //  which hemisphere
        return mris->hemisphere;
    }
    int Surface::initialized() const {
        return mris->initialized;
    }
    PLTA Surface::lta() const {
        return mris->lta;
    }
    PMATRIX Surface::SRASToTalSRAS_() const {
        return mris->SRASToTalSRAS_;
    }
    PMATRIX Surface::TalSRASToSRAS_() const {
        return mris->TalSRASToSRAS_;
    }
    int Surface::free_transform() const {
        return mris->free_transform;
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return mris->radius;
    }
    float Surface::a() const {
        return mris->a;
    }
    float Surface::b() const {
        return mris->b;
    }
    float Surface::c() const {  //  ellipsoid parameters
        return mris->c;
    }
    MRIS_fname_t Surface::fname() const {  //  file it was originally loaded from
        return mris->fname;
    }
    float Surface::Hmin() const {  //  min mean curvature
        return mris->Hmin;
    }
    float Surface::Hmax() const {  //  max mean curvature
        return mris->Hmax;
    }
    float Surface::Kmin() const {  //  min Gaussian curvature
        return mris->Kmin;
    }
    float Surface::Kmax() const {  //  max Gaussian curvature
        return mris->Kmax;
    }
    double Surface::Ktotal() const {  //  total Gaussian curvature
        return mris->Ktotal;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return mris->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return mris->origxyz_status;
    }
    int Surface::patch() const {  //  if a patch of the surface
        return mris->patch;
    }
    int Surface::nlabels() const {
        return mris->nlabels;
    }
    PMRIS_AREA_LABEL Surface::labels() const {  //  nlabels of these (may be null)
        return mris->labels;
    }
    p_void Surface::vp() const {  //  for misc. use
        return mris->vp;
    }
    float Surface::alpha() const {  //  rotation around z-axis
        return mris->alpha;
    }
    float Surface::beta() const {  //  rotation around y-axis
        return mris->beta;
    }
    float Surface::gamma() const {  //  rotation around x-axis
        return mris->gamma;
    }
    float Surface::da() const {
        return mris->da;
    }
    float Surface::db() const {
        return mris->db;
    }
    float Surface::dg() const {  //  old deltas
        return mris->dg;
    }
    int Surface::type() const {  //  what type of surface was this initially
        return mris->type;
    }
    int Surface::max_vertices() const {  //  may be bigger than nvertices, set by calling MRISreallocVerticesAndFaces
        return mris->max_vertices;
    }
    int Surface::max_faces() const {  //  may be bigger than nfaces,    set by calling MRISreallocVerticesAndFaces
        return mris->max_faces;
    }
    MRIS_subject_name_t Surface::subject_name() const {  //  name of the subject
        return mris->subject_name;
    }
    float Surface::canon_area() const {
        return mris->canon_area;
    }
    int Surface::noscale() const {  //  don't scale by surface area if true
        return mris->noscale;
    }
    float Surface::dx2(size_t i) const {  //  an extra set of gradient (not always alloced)
        return mris->dx2[i];
    }
    float Surface::dy2(size_t i) const {
        return mris->dy2[i];
    }
    float Surface::dz2(size_t i) const {
        return mris->dz2[i];
    }
    PCOLOR_TABLE Surface::ct() const {
        return mris->ct;
    }
    int Surface::useRealRAS() const {  //  if 0 (default), vertex position is a conformed volume RAS with c_(r,"a","s")=0.  else is a real RAS (volume stored RAS)
        return mris->useRealRAS;
    }
    VOL_GEOM Surface::vg() const {  //  volume info from which this surface is created. valid iff vg.valid = 1
        return mris->vg;
    }
    MRIS_cmdlines_t Surface::cmdlines() const {
        return mris->cmdlines;
    }
    int Surface::ncmds() const {
        return mris->ncmds;
    }
    float Surface::group_avg_surface_area() const {  //  average of total surface area for group
        return mris->group_avg_surface_area;
    }
    int Surface::group_avg_vtxarea_loaded() const {  //  average vertex area for group at each vertex
        return mris->group_avg_vtxarea_loaded;
    }
    int Surface::triangle_links_removed() const {  //  for quad surfaces
        return mris->triangle_links_removed;
    }
    p_void Surface::user_parms() const {  //  for whatever the user wants to hang here
        return mris->user_parms;
    }
    PMATRIX Surface::m_sras2vox() const {  //  for converting surface ras to voxel
        return mris->m_sras2vox;
    }
    PMRI Surface::mri_sras2vox() const {  //  volume that the above matrix is for
        return mris->mri_sras2vox;
    }
    p_void Surface::mht() const {
        return mris->mht;
    }
    p_void Surface::temps() const {
        return mris->temps;
    }
    
    void Surface::set_strips(size_t i, STRIP to) {
        mris->strips[i] = to;
    }
    void Surface::set_xctr(float to) {
        mris->xctr = to;
    }
    void Surface::set_yctr(float to) {
        mris->yctr = to;
    }
    void Surface::set_zctr(float to) {
        mris->zctr = to;
    }
    void Surface::set_xlo(float to) {
        mris->xlo = to;
    }
    void Surface::set_ylo(float to) {
        mris->ylo = to;
    }
    void Surface::set_zlo(float to) {
        mris->zlo = to;
    }
    void Surface::set_xhi(float to) {
        mris->xhi = to;
    }
    void Surface::set_yhi(float to) {
        mris->yhi = to;
    }
    void Surface::set_zhi(float to) {
        mris->zhi = to;
    }
    void Surface::set_x0(float to) {  //  center of spherical expansion
        mris->x0 = to;
    }
    void Surface::set_y0(float to) {
        mris->y0 = to;
    }
    void Surface::set_z0(float to) {
        mris->z0 = to;
    }
    void Surface::set_v_temporal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_temporal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_frontal_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_frontal_pole = mris->vertices + to.idx;
    }
    void Surface::set_v_occipital_pole(Vertex to) {
        cheapAssert(mris == to.mris); mris->v_occipital_pole = mris->vertices + to.idx;
    }
    void Surface::set_max_curv(float to) {
        mris->max_curv = to;
    }
    void Surface::set_min_curv(float to) {
        mris->min_curv = to;
    }
    void Surface::set_total_area(float to) {
        mris->total_area = to;
    }
    void Surface::set_avg_vertex_area(double to) {
        mris->avg_vertex_area = to;
    }
    void Surface::set_avg_vertex_dist(double to) {  //  set by MRIScomputeAvgInterVertexDist
        mris->avg_vertex_dist = to;
    }
    void Surface::set_std_vertex_dist(double to) {
        mris->std_vertex_dist = to;
    }
    void Surface::set_orig_area(float to) {
        mris->orig_area = to;
    }
    void Surface::set_neg_area(float to) {
        mris->neg_area = to;
    }
    void Surface::set_neg_orig_area(float to) {  //  amount of original surface in folds
        mris->neg_orig_area = to;
    }
    void Surface::set_zeros(int to) {
        mris->zeros = to;
    }
    void Surface::set_hemisphere(int to) {  //  which hemisphere
        mris->hemisphere = to;
    }
    void Surface::set_fname(MRIS_fname_t to) {  //  file it was originally loaded from
        mris->fname = to;
    }
    void Surface::set_Hmin(float to) {  //  min mean curvature
        mris->Hmin = to;
    }
    void Surface::set_Hmax(float to) {  //  max mean curvature
        mris->Hmax = to;
    }
    void Surface::set_Kmin(float to) {  //  min Gaussian curvature
        mris->Kmin = to;
    }
    void Surface::set_Kmax(float to) {  //  max Gaussian curvature
        mris->Kmax = to;
    }
    void Surface::set_Ktotal(double to) {  //  total Gaussian curvature
        mris->Ktotal = to;
    }
    void Surface::set_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane)
        mris->status = to;
    }
    void Surface::set_origxyz_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        mris->origxyz_status = to;
    }
    void Surface::set_patch(int to) {  //  if a patch of the surface
        mris->patch = to;
    }
    void Surface::set_nlabels(int to) {
        mris->nlabels = to;
    }
    void Surface::set_labels(PMRIS_AREA_LABEL to) {  //  nlabels of these (may be null)
        mris->labels = to;
    }
    void Surface::set_vp(p_void to) {  //  for misc. use
        mris->vp = to;
    }
    void Surface::set_alpha(float to) {  //  rotation around z-axis
        mris->alpha = to;
    }
    void Surface::set_beta(float to) {  //  rotation around y-axis
        mris->beta = to;
    }
    void Surface::set_gamma(float to) {  //  rotation around x-axis
        mris->gamma = to;
    }
    void Surface::set_da(float to) {
        mris->da = to;
    }
    void Surface::set_db(float to) {
        mris->db = to;
    }
    void Surface::set_dg(float to) {  //  old deltas
        mris->dg = to;
    }
    void Surface::set_type(int to) {  //  what type of surface was this initially
        mris->type = to;
    }


    } // namespace AllM
} // namespace SurfaceFromMRIS
