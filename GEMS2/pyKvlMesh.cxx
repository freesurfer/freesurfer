#include <pybind11/pybind11.h>
#include "pyKvlMesh.h"

namespace py = pybind11;

KvlMeshCollection::KvlMeshCollection(const std::string &meshFileName) {
    std::cout << "Read mesh: " << meshFileName << std::endl;
    meshCollection = 0;
}
