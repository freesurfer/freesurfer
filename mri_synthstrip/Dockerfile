FROM centos:7

# copy local data
COPY . /external

# shell settings
WORKDIR /freesurfer

# install utils
RUN yum -y update
RUN yum -y install libgomp gcc python3 python3-devel
RUN yum clean all

# python packages
RUN python3 -m pip install -U pip
RUN python3 -m pip install scipy torch==1.10.2
RUN python3 -m pip install surfa
RUN python3 -m pip install cache purge

# install synthstrip
RUN cp /external/mri_synthstrip /freesurfer/

# configure model
ENV FREESURFER_HOME /freesurfer
RUN mkdir -p /freesurfer/models
RUN cp /external/synthstrip.1.pt /freesurfer/models/

# clean up
RUN rm -rf /external /root/.cache/pip

ENTRYPOINT ["python3", "/freesurfer/mri_synthstrip"]
