#include "pyKvlTransform.h"
#include "pyKvlNumpy.h"


KvlTransform::KvlTransform(const py::array_t<double> &transformMatrix) {
    m_transform = TransformType::New();
    TransformType::ParametersType  parameters( 12 );
    const double* const data = transformMatrix.data(0);

    for ( unsigned int row = 0; row < 3; row++ )
    {
        for ( unsigned int col = 0; col < 3; col++ )
        {
            parameters[ row * 3 + col ] = data[ col * 4 + row ];
        }
        parameters[ 9 + row ] = data[ 12 + row ];
    }
    m_transform->SetParameters( parameters );
}

py::array_t<double> KvlTransform::AsNumpyArray() const {
    auto *data = new double[16];
    auto parameters = m_transform->GetParameters();
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
