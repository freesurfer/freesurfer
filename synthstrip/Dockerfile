FROM centos:7

# copy local data
COPY . /external

# shell settings
WORKDIR /freesurfer

# install utils
RUN yum -y update
RUN yum -y install libgomp python3
RUN yum clean all

# python packages
RUN pip3 install scipy torch==1.10.2
RUN pip3 install /external/fsmodule

# install synthstrip
RUN cp /external/synthstrip /freesurfer/

# configure model
ENV FREESURFER_HOME /freesurfer
RUN mkdir -p /freesurfer/models
RUN cp /external/synthstrip.1.pt /freesurfer/models/

# setup rest of the env
ENV OS Linux
ENV FSF_OUTPUT_FORMAT nii.gz

# clean up
RUN rm -rf /external /root/.cache/pip

ENTRYPOINT ["python3", "/freesurfer/synthstrip"]
