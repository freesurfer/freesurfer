
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
namespace SurfaceFromMRISPV {
    typedef MRISPV Representation;


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

        friend struct SurfaceFromMRISPV::ExistenceM::Face;
        friend struct SurfaceFromMRISPV::ExistenceM::Vertex;
        friend struct SurfaceFromMRISPV::ExistenceM::Surface;
        friend struct SurfaceFromMRISPV::Existence::Face;
        friend struct SurfaceFromMRISPV::Existence::Vertex;
        friend struct SurfaceFromMRISPV::Existence::Surface;
        friend struct SurfaceFromMRISPV::TopologyM::Face;
        friend struct SurfaceFromMRISPV::TopologyM::Vertex;
        friend struct SurfaceFromMRISPV::TopologyM::Surface;
        friend struct SurfaceFromMRISPV::Topology::Face;
        friend struct SurfaceFromMRISPV::Topology::Vertex;
        friend struct SurfaceFromMRISPV::Topology::Surface;
        friend struct SurfaceFromMRISPV::XYZPositionM::Face;
        friend struct SurfaceFromMRISPV::XYZPositionM::Vertex;
        friend struct SurfaceFromMRISPV::XYZPositionM::Surface;
        friend struct SurfaceFromMRISPV::XYZPosition::Face;
        friend struct SurfaceFromMRISPV::XYZPosition::Vertex;
        friend struct SurfaceFromMRISPV::XYZPosition::Surface;
        friend struct SurfaceFromMRISPV::XYZPositionConsequencesM::Face;
        friend struct SurfaceFromMRISPV::XYZPositionConsequencesM::Vertex;
        friend struct SurfaceFromMRISPV::XYZPositionConsequencesM::Surface;
        friend struct SurfaceFromMRISPV::XYZPositionConsequences::Face;
        friend struct SurfaceFromMRISPV::XYZPositionConsequences::Vertex;
        friend struct SurfaceFromMRISPV::XYZPositionConsequences::Surface;
        friend struct SurfaceFromMRISPV::DistortM::Face;
        friend struct SurfaceFromMRISPV::DistortM::Vertex;
        friend struct SurfaceFromMRISPV::DistortM::Surface;
        friend struct SurfaceFromMRISPV::Distort::Face;
        friend struct SurfaceFromMRISPV::Distort::Vertex;
        friend struct SurfaceFromMRISPV::Distort::Surface;
        friend struct SurfaceFromMRISPV::AnalysisM::Face;
        friend struct SurfaceFromMRISPV::AnalysisM::Vertex;
        friend struct SurfaceFromMRISPV::AnalysisM::Surface;
        friend struct SurfaceFromMRISPV::Analysis::Face;
        friend struct SurfaceFromMRISPV::Analysis::Vertex;
        friend struct SurfaceFromMRISPV::Analysis::Surface;
        friend struct SurfaceFromMRISPV::AllM::Face;
        friend struct SurfaceFromMRISPV::AllM::Vertex;
        friend struct SurfaceFromMRISPV::AllM::Surface;
    };
} // SurfaceFromMRISPV
