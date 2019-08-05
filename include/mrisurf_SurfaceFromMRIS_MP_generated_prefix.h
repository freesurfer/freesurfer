
#pragma once
// GENERATED SOURCE - DO NOT DIRECTLY EDIT
// 
// =======================================
#include "mrisurf_aaa.h"
namespace SurfaceFromMRIS_MP {
    typedef MRIS_MP Representation;


    namespace Existence {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace Existence


    namespace Topology {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace Topology


    namespace XYZPosition {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace XYZPosition


    namespace XYZPositionConsequences {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace XYZPositionConsequences


    namespace Distort {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace Distort


    namespace Analysis {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace Analysis


    namespace ExistenceM {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace ExistenceM


    namespace TopologyM {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace TopologyM


    namespace XYZPositionM {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace XYZPositionM


    namespace XYZPositionConsequencesM {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace XYZPositionConsequencesM


    namespace DistortM {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace DistortM


    namespace AnalysisM {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
        struct Surface;
    } // namespace AnalysisM


    namespace AllM {
        struct Face;
        struct Vertex;
        struct MRIS_MP;
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

        friend struct SurfaceFromMRIS_MP::ExistenceM::Face;
        friend struct SurfaceFromMRIS_MP::ExistenceM::Vertex;
        friend struct SurfaceFromMRIS_MP::ExistenceM::Surface;
        friend struct SurfaceFromMRIS_MP::Existence::Face;
        friend struct SurfaceFromMRIS_MP::Existence::Vertex;
        friend struct SurfaceFromMRIS_MP::Existence::Surface;
        friend struct SurfaceFromMRIS_MP::TopologyM::Face;
        friend struct SurfaceFromMRIS_MP::TopologyM::Vertex;
        friend struct SurfaceFromMRIS_MP::TopologyM::Surface;
        friend struct SurfaceFromMRIS_MP::Topology::Face;
        friend struct SurfaceFromMRIS_MP::Topology::Vertex;
        friend struct SurfaceFromMRIS_MP::Topology::Surface;
        friend struct SurfaceFromMRIS_MP::XYZPositionM::Face;
        friend struct SurfaceFromMRIS_MP::XYZPositionM::Vertex;
        friend struct SurfaceFromMRIS_MP::XYZPositionM::Surface;
        friend struct SurfaceFromMRIS_MP::XYZPosition::Face;
        friend struct SurfaceFromMRIS_MP::XYZPosition::Vertex;
        friend struct SurfaceFromMRIS_MP::XYZPosition::Surface;
        friend struct SurfaceFromMRIS_MP::XYZPositionConsequencesM::Face;
        friend struct SurfaceFromMRIS_MP::XYZPositionConsequencesM::Vertex;
        friend struct SurfaceFromMRIS_MP::XYZPositionConsequencesM::Surface;
        friend struct SurfaceFromMRIS_MP::XYZPositionConsequences::Face;
        friend struct SurfaceFromMRIS_MP::XYZPositionConsequences::Vertex;
        friend struct SurfaceFromMRIS_MP::XYZPositionConsequences::Surface;
        friend struct SurfaceFromMRIS_MP::DistortM::Face;
        friend struct SurfaceFromMRIS_MP::DistortM::Vertex;
        friend struct SurfaceFromMRIS_MP::DistortM::Surface;
        friend struct SurfaceFromMRIS_MP::Distort::Face;
        friend struct SurfaceFromMRIS_MP::Distort::Vertex;
        friend struct SurfaceFromMRIS_MP::Distort::Surface;
        friend struct SurfaceFromMRIS_MP::AnalysisM::Face;
        friend struct SurfaceFromMRIS_MP::AnalysisM::Vertex;
        friend struct SurfaceFromMRIS_MP::AnalysisM::Surface;
        friend struct SurfaceFromMRIS_MP::Analysis::Face;
        friend struct SurfaceFromMRIS_MP::Analysis::Vertex;
        friend struct SurfaceFromMRIS_MP::Analysis::Surface;
        friend struct SurfaceFromMRIS_MP::AllM::Face;
        friend struct SurfaceFromMRIS_MP::AllM::Vertex;
        friend struct SurfaceFromMRIS_MP::AllM::Surface;
    };
} // SurfaceFromMRIS_MP
