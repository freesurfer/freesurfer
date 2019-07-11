
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
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
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
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

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
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
        CASE(VERTEX_NORMALS,n)
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
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
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
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
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

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
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
        CASE(VERTEX_NORMALS,n)
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
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
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
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
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
        CASE(VERTEX_NORMALS,n)
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
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
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
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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


    Surface::Surface (                                ) {}
    Surface::Surface ( Representation* representation ) : Repr_Elt(representation,0) {}
    Surface::Surface ( Surface const & src            ) : Repr_Elt(src) {}
    Surface::Surface ( AllM::Surface const & src      ) : Repr_Elt(src) {}

    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
    }
    
    void Surface::set_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane)
        repr->status = to;
    }
    void Surface::set_origxyz_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        repr->origxyz_status = to;
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

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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
    void Vertex::set_x(float to) {  //  current coordinates	
        repr->v_x[idx] = to;
    }
    void Vertex::set_y(float to) {  //  use MRISsetXYZ() to set
        repr->v_y[idx] = to;
    }
    void Vertex::set_z(float to) {
        repr->v_z[idx] = to;
    }
    void Vertex::set_ripflag(char to) {  //  vertex no longer exists - placed last to load the next vertex into cache
        repr->v_ripflag[idx] = to;
    }


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
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
    }
    double Surface::radius() const {  //  radius (if status==MRIS_SPHERE)
        return repr->radius;
    }
    MRIS_Status Surface::status() const {  //  type of surface (e.g. sphere, plane)
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_area(float to) {
        repr->f_area[idx] = to;
    }
    void Face::set_angle(angles_per_triangle_t to) {
        repr->f_angle[idx] = to;
    }
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_norm(PDMATRIX to) {
        repr->f_norm[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
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
        CASE(VERTEX_NORMALS,n)
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
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
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
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
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
        CASE(VERTEX_NORMALS,n)
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
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
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
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    PDMATRIX Face::norm() const {
        return repr->f_norm[idx];
    }
    
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}
    Vertex::Vertex ( AllM::Vertex const & src                   ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
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
        CASE(VERTEX_NORMALS,n)
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
    void Vertex::set_nx(float to) {
        repr->v_nx[idx] = to;
    }
    void Vertex::set_ny(float to) {
        repr->v_ny[idx] = to;
    }
    void Vertex::set_nz(float to) {  //  curr normal
        repr->v_nz[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
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
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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
    char Face::ripflag() const {
        return repr->f_ripflag[idx];
    }
    PDMATRIX Face::norm() const {
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
    void Face::set_ripflag(char to) {
        repr->f_ripflag[idx] = to;
    }
    void Face::set_norm(PDMATRIX to) {
        repr->f_norm[idx] = to;
    }


    Vertex::Vertex (                                            ) {}
    Vertex::Vertex ( Representation* representation, size_t idx ) : Repr_Elt(representation,idx) {}
    Vertex::Vertex ( Vertex const & src                         ) : Repr_Elt(src) {}

    Face Vertex::f(size_t i) const {  // size() is num.    array[v->num] the fno's of the neighboring faces         
        return Face(repr,repr->v_f[idx][i]);
    }
    uchar Vertex::num() const {  //  number of neighboring faces                              
        return repr->v_num[idx];
    }
    float Vertex::dist(size_t i) const {  // size() is vtotal.    distance to neighboring vertices based on  xyz   
        return repr->v_dist[idx][i];
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
    float Vertex::nx() const {
        return repr->v_nx[idx];
    }
    float Vertex::ny() const {
        return repr->v_ny[idx];
    }
    float Vertex::nz() const {  //  curr normal
        return repr->v_nz[idx];
    }
    float Vertex::area() const {
        return repr->v_area[idx];
    }
    float Vertex::origarea() const {
        return repr->v_origarea[idx];
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
        CASE(VERTEX_NORMALS,n)
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
    void Vertex::set_area(float to) {
        repr->v_area[idx] = to;
    }
    void Vertex::set_origarea(float to) {
        repr->v_origarea[idx] = to;
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

    int Surface::nvertices() const {  //  # of vertices on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nvertices;
    }
    int Surface::nfaces() const {  //  # of faces on surface, change by calling MRISreallocVerticesAndFaces et al
        return repr->nfaces;
    }
    Vertex Surface::vertices(size_t i) const {
        return Vertex(repr,i);
    }
    Face Surface::faces(size_t i) const {
        return Face(repr,i);
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
        return repr->status;
    }
    MRIS_Status Surface::origxyz_status() const {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        return repr->origxyz_status;
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
    void Surface::set_neg_area(float to) {
        repr->neg_area = to;
    }
    void Surface::set_neg_orig_area(float to) {  //  amount of original surface in folds
        repr->neg_orig_area = to;
    }
    void Surface::set_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane)
        repr->status = to;
    }
    void Surface::set_origxyz_status(MRIS_Status to) {  //  type of surface (e.g. sphere, plane) that this origxyz were obtained from
        repr->origxyz_status = to;
    }


    } // namespace AllM
} // SurfaceFromMRISPV
