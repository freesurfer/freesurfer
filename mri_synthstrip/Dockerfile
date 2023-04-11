# base image with Python and SynthStrip on path (eventually)
FROM ubuntu:22.04 AS base
ENV FREESURFER_HOME="/freesurfer"
ENV VIRTUAL_ENV="$FREESURFER_HOME/env"
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN apt-get update && \
    apt-get install -y --no-install-recommends python3 && \
    rm -rf /var/lib/apt/lists/*


# intermediate stage with build requirements not needed in final image
FROM base AS build
RUN apt-get update && \
    apt-get install -y build-essential python3-pip python3-dev python3-venv

# install packages into virtual environment as soon as it exists
RUN python3 -m pip install --upgrade pip
RUN python3 -m venv "$VIRTUAL_ENV"
RUN python3 -m pip install torch==2.0.0 surfa

# install synthstrip from local folder
COPY mri_synthstrip $VIRTUAL_ENV/bin
COPY synthstrip.*.pt $FREESURFER_HOME/models/


# only copy files needed for final image
FROM base
COPY --from=build $FREESURFER_HOME $FREESURFER_HOME
ENTRYPOINT ["mri_synthstrip"]
