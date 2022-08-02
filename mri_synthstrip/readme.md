# SynthStrip

Primary source code for the SynthStrip skull-stripping utility. For command documentation, see https://surfer.nmr.mgh.harvard.edu/docs/synthstrip.

## Building New Containers

The `mri_synthstrip` command is automatically built into freesurfer releases, but in order to build and release new standalone SynthStrip containers, run the following:

```bash
# first make sure the models are update and NOT symlinked
git fetch datasrc
git annex get .
git annex unlock synthstrip.*.pt

# then run the docker build, make sure the new version is tagged appropriately
VERSION=X.X
docker build -t freesurfer/synthstrip:${VERSION} .
docker push freesurfer/synthstrip:${VERSION}
```

Once this is built and pushed to the [freesurfer dockerhub](https://hub.docker.com/u/freesurfer), make sure to update the version reference in the `synthstrip-docker` and `synthstrip-singularity` scripts.