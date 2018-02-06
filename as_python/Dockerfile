FROM ubuntu:xenial

RUN apt-get update && \
    apt-get install -y build-essential \
                       tcsh \
                       libtool-bin \
                       libtool \
                       automake \
                       gfortran \
                       libglu1-mesa-dev \
                       libfreetype6-dev \
                       uuid-dev \
                       libxmu-dev \
                       libxmu-headers \
                       libxi-dev \
                       libx11-dev \
                       libxml2-utils \
                       libxt-dev \
                       libjpeg62-dev \
                       libxaw7-dev \
                       liblapack-dev \
                       git \
                       gcc-4.8 \
                       g++-4.8 \
                       libgfortran-4.8-dev \
                       curl \
                       python-pip \
                       python3-pip

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 50 && \
    update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.8 50

WORKDIR /home/freesurfer

#Install cmake
RUN curl -O https://cmake.org/files/v3.10/cmake-3.10.2-Linux-x86_64.sh && sh ./cmake-3.10.2-Linux-x86_64.sh --skip-license && cp -r bin /usr/ && cp -r doc /usr/share/ && cp -r man /usr/share/ && cp -r share /usr/

# Install ITK
RUN git clone https://itk.org/ITK.git
RUN cd /home/freesurfer && mkdir ITK-build && cd ITK-build && \
    cmake ../ITK \
          -DITK_BUILD_DEFAULT_MODULES=OFF \
          -DITKGroup_Core=ON \
          -DITKGroup_Filtering=ON \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_FLAGS="-msse2 -mfpmath=sse" \
          -DCMAKE_C_FLAGS="-msse2 -mfpmath=sse" \
    && make -j8
ENV ITK_DIR /home/freesurfer/ITK-build

# Get GEMS python code and requirements
RUN git clone https://github.com/innolitics/freesurfer.git && \
    cd freesurfer && \
    git checkout  nf-gems2-python-port && \
    cd GEMS2 && \
    git clone https://github.com/pybind/pybind11.git
RUN pip3 install --upgrade pip && pip3 install -r /home/freesurfer/freesurfer/as_python/requirements.txt

# Build GEMS
RUN cd /home/freesurfer/freesurfer/GEMS2 && \
    cmake -D CMAKE_CXX_FLAGS="-msse2 -mfpmath=sse -fPIC -fpermissive" \
          -D CMAKE_C_FLAGS="-msse2 -mfpmath=sse -fPIC -fpermissive" \
          -D BUILD_EXECUTABLES=OFF \
          -D BUILD_GUI=OFF \
          -D BUILD_MATLAB=OFF \
          -D BUILD_SHARED_LIBS=OFF \
          -D BUILD_TESTING=OFF \
          -D PYTHON_EXECUTABLE=/usr/bin/python3.5 \
          -D PYTHON_LIBRARY=/usr/lib/x86_64-linux-gnu/libpython3.5m.so \
          . && \
    make -j8

WORKDIR /home/freesurfer/freesurfer
ENV PYTHONPATH .:./GEMS2/bin

