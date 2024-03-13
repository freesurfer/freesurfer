# base image with Python and SynthStrip on path (eventually)
FROM ubuntu:22.04 AS base
ENV FREESURFER_HOME="/freesurfer"
ENV PYTHONUSERBASE="$FREESURFER_HOME/env"
ENV PATH="$FREESURFER_HOME:$PATH"

RUN apt-get update && \
    apt-get install -y --no-install-recommends python3 && \
    rm -rf /var/lib/apt/lists/*


# intermediate stage with build requirements not needed in final image
FROM base AS build
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    git \
    python3-dev \
    python3-pip

# install packages into user base simple easy COPY
RUN python3 -m pip install -U pip
RUN python3 -m pip install --user \
    torch==2.1.2 \
    git+https://github.com/freesurfer/surfa.git@0d83332351083b33c4da221e9d10a63a93ae7f52

# install synthstrip from local folder
COPY --chmod=0775 mri_synthstrip $FREESURFER_HOME
COPY --chmod=0664 synthstrip.*.pt $FREESURFER_HOME/models/

# export build artifacts
RUN python3 -V >python.txt
RUN python3 -m pip freeze >requirements.txt


FROM scratch AS export
COPY --from=build *.txt /


# only copy files needed and only once to avoid unnecessary layers
FROM base
COPY --from=build $FREESURFER_HOME $FREESURFER_HOME
ENTRYPOINT ["mri_synthstrip"]
