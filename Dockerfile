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
RUN pip3 install --upgrade pip

ARG working_dir=/home/freesurfer
WORKDIR $working_dir

# Get GEMS python code and requirements
COPY GEMS2 $working_dir/GEMS2
COPY as_python $working_dir/as_python

RUN cd $working_dir/GEMS2 && git clone https://github.com/pybind/pybind11.git
RUN pip3 install -r $working_dir/as_python/requirements.txt

#Install cmake
RUN curl -O https://cmake.org/files/v3.10/cmake-3.10.2-Linux-x86_64.sh && sh ./cmake-3.10.2-Linux-x86_64.sh --skip-license && cp -r bin /usr/ && cp -r doc /usr/share/ && cp -r man /usr/share/ && cp -r share /usr/

# Install ITK
RUN git clone https://itk.org/ITK.git
RUN cd $working_dir && mkdir ITK-build && cd ITK-build && \
    cmake ../ITK \
          -DITK_BUILD_DEFAULT_MODULES=OFF \
          -DITKGroup_Core=ON \
          -DITKGroup_Filtering=ON \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_CXX_FLAGS="-msse2 -mfpmath=sse" \
          -DCMAKE_C_FLAGS="-msse2 -mfpmath=sse" \
    && make -j8
ENV ITK_DIR $working_dir/ITK-build

# Build GEMS
RUN cd $working_dir/GEMS2 && \
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

ENV PYTHONPATH .:./GEMS2/bin
ENV SAMSEG_DATA_DIR /data/samseg_data
