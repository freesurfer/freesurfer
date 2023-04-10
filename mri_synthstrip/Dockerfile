FROM ubuntu:18.04

# copy local data
COPY . /external

# shell settings
WORKDIR /freesurfer

# install utils
RUN apt-get update && \
    apt-get install -y build-essential python3 python3-pip python3-dev && \
    rm -rf /var/lib/apt/lists/*

# python packages
RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install numpy torch==1.10.2
RUN python3 -m pip install surfa
RUN python3 -m pip install cache purge

# install synthstrip
RUN cp /external/mri_synthstrip /freesurfer/

# configure model
ENV FREESURFER_HOME /freesurfer
RUN mkdir -p /freesurfer/models
RUN cp /external/synthstrip.*.pt /freesurfer/models/

# clean up
RUN rm -rf /external /root/.cache/pip

ENTRYPOINT ["python3", "/freesurfer/mri_synthstrip"]
