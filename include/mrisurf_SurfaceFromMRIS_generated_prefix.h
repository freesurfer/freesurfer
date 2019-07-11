
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
namespace SurfaceFromMRIS {
    typedef MRIS Representation;


    namespace Existence {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace Existence


    namespace Topology {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace Topology


    namespace XYZPosition {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace XYZPosition


    namespace XYZPositionConsequences {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace XYZPositionConsequences


    namespace Distort {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace Distort


    namespace Analysis {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace Analysis


    namespace ExistenceM {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace ExistenceM


    namespace TopologyM {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace TopologyM


    namespace XYZPositionM {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace XYZPositionM


    namespace XYZPositionConsequencesM {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace XYZPositionConsequencesM


    namespace DistortM {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace DistortM


    namespace AnalysisM {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace AnalysisM


    namespace AllM {
        struct Face;
        struct Vertex;
        struct Surface;
    } // namespace AllM
    
    struct Repr_Elt { 
        bool operator==(Repr_Elt const & rhs) const { return repr == rhs.repr && idx == rhs.idx; }
        bool operator!=(Repr_Elt const & rhs) const { return repr != rhs.repr || idx != rhs.idx; }
    protected: 
        Representation* repr; size_t idx; 
        Repr_Elt() : repr(nullptr), idx(0) {}
        Repr_Elt(Representation* repr, size_t idx) : repr(repr), idx(idx) {}
        Repr_Elt(Repr_Elt const & src) : repr(src.repr), idx(src.idx) {}

        friend struct SurfaceFromMRIS::ExistenceM::Face;
        friend struct SurfaceFromMRIS::ExistenceM::Vertex;
        friend struct SurfaceFromMRIS::ExistenceM::Surface;
        friend struct SurfaceFromMRIS::Existence::Face;
        friend struct SurfaceFromMRIS::Existence::Vertex;
        friend struct SurfaceFromMRIS::Existence::Surface;
        friend struct SurfaceFromMRIS::TopologyM::Face;
        friend struct SurfaceFromMRIS::TopologyM::Vertex;
        friend struct SurfaceFromMRIS::TopologyM::Surface;
        friend struct SurfaceFromMRIS::Topology::Face;
        friend struct SurfaceFromMRIS::Topology::Vertex;
        friend struct SurfaceFromMRIS::Topology::Surface;
        friend struct SurfaceFromMRIS::XYZPositionM::Face;
        friend struct SurfaceFromMRIS::XYZPositionM::Vertex;
        friend struct SurfaceFromMRIS::XYZPositionM::Surface;
        friend struct SurfaceFromMRIS::XYZPosition::Face;
        friend struct SurfaceFromMRIS::XYZPosition::Vertex;
        friend struct SurfaceFromMRIS::XYZPosition::Surface;
        friend struct SurfaceFromMRIS::XYZPositionConsequencesM::Face;
        friend struct SurfaceFromMRIS::XYZPositionConsequencesM::Vertex;
        friend struct SurfaceFromMRIS::XYZPositionConsequencesM::Surface;
        friend struct SurfaceFromMRIS::XYZPositionConsequences::Face;
        friend struct SurfaceFromMRIS::XYZPositionConsequences::Vertex;
        friend struct SurfaceFromMRIS::XYZPositionConsequences::Surface;
        friend struct SurfaceFromMRIS::DistortM::Face;
        friend struct SurfaceFromMRIS::DistortM::Vertex;
        friend struct SurfaceFromMRIS::DistortM::Surface;
        friend struct SurfaceFromMRIS::Distort::Face;
        friend struct SurfaceFromMRIS::Distort::Vertex;
        friend struct SurfaceFromMRIS::Distort::Surface;
        friend struct SurfaceFromMRIS::AnalysisM::Face;
        friend struct SurfaceFromMRIS::AnalysisM::Vertex;
        friend struct SurfaceFromMRIS::AnalysisM::Surface;
        friend struct SurfaceFromMRIS::Analysis::Face;
        friend struct SurfaceFromMRIS::Analysis::Vertex;
        friend struct SurfaceFromMRIS::Analysis::Surface;
        friend struct SurfaceFromMRIS::AllM::Face;
        friend struct SurfaceFromMRIS::AllM::Vertex;
        friend struct SurfaceFromMRIS::AllM::Surface;
    };
} // SurfaceFromMRIS
