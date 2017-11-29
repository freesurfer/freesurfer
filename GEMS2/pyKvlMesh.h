#ifndef GEMS_PYKVLMESH_H_H
#define GEMS_PYKVLMESH_H_H

#include "itkObject.h"
#include "kvlAtlasMeshCollection.h"

class KvlMeshCollection {
    kvl::AtlasMeshCollection::Pointer  meshCollection;
public:
    KvlMeshCollection(const std::string &meshFileName);
    bool isEmpty() {
        return meshCollection.IsNull();
    }
};

#endif //GEMS_PYKVLMESH_H_H
