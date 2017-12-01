#include "pyKvlTransform.h"
#include "pyKvlNumpy.h"


py::array_t<double> TransformToNumpy(TransformPointer transform) {
    auto *data = new double[16];
    auto parameters = transform->GetParameters();
    for ( unsigned int row = 0; row < 3; row++ )
    {
        for ( unsigned int col = 0; col < 3; col++ )
        {
            data[ col * 4 + row ] = parameters[ row * 3 + col ];
        }
        data[ 12 + row ] = parameters[ 9 + row ];
    }
    for ( unsigned int col = 0; col < 3; col++ )
    {
        data[ col * 4 + 3 ] = 0.0f;
    }
    data[ 15 ] = 1.0f;
    return createNumpyArrayFStyle({4, 4}, data);
}

TransformPointer NumpyToTransform(py::array_t<double> transform) {
    return NULL;
}

KvlTransform::KvlTransform(const py::array_t<double> &transformMatrix) {
    // Create the ITK transform object and fill in its elements
    m_transform = TransformType::New();
    TransformType::ParametersType  parameters( 12 );
    for ( unsigned int row = 0; row < 3; row++ )
    {
        for ( unsigned int col = 0; col < 3; col++ )
        {
            std::cout << transformMatrix.at(row, col);
            parameters[ row * 3 + col ] = transformMatrix.at(row, col);
        }
        parameters[ 9 + row ] =  transformMatrix.at(row, 4);
    }
    m_transform->SetParameters( parameters );
}