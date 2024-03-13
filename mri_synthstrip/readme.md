# SynthStrip

Source code for the SynthStrip skull-stripping utility. For command documentation, see https://w3id.org/synthstrip.


## Building containers

The `mri_synthstrip` command is automatically built into FreeSurfer.
To create a standalone SynthStrip container, first fetch the latest model files and unlock them to replace the `git-annex` symlinks:

```shell
git fetch datasrc
git annex get .
git annex unlock synthstrip.*.pt
```

Then build and push the container to the [Docker Hub](https://hub.docker.com/u/freesurfer), tagging the new version appropriately:

```shell
VERSION=X.X
docker build -t freesurfer/synthstrip:${VERSION} .
docker push freesurfer/synthstrip:${VERSION}
```

Remember to point the default "latest" tag to the new container.

```shell
docker tag freesurfer/synthstrip:${VERSION} freesurfer/synthstrip:latest
docker push freesurfer/synthstrip:latest
```

Finally, update the version reference in the `synthstrip-docker` and `synthstrip-singularity` scripts:

```shell
sed -i "s/^\(version = \).*/\1'$VERSION'/" synthstrip-docker
sed -i "s/\(synthstrip.\)[0-9.]*\(\.\|$\)/\1$VERSION\2/g" synthstrip-singularity
git diff synthstrip-*
```

## Exporting requirements

Export updated requirement files as build artifacts for users who wish to build environments:

```shell
docker build --target export --output env .
```
