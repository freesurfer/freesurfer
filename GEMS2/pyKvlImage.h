#ifndef GEMS_KVLIMAGE_H
#define GEMS_KVLIMAGE_H

#include "itkImage.h"
#include "pyKvlTransform.h"

typedef itk::Image< float, 3 >  ImageType;
typedef ImageType::Pointer ImagePointer;

class KvlImage {
    ImagePointer imageHandle;
    TransformPointer transform;

public:
    KvlImage(const std::string &imageFileName);
    void greet();
};

void useImage(KvlImage* image);



#endif //GEMS_KVLIMAGE_H
