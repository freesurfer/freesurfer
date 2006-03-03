The following files, found in the directory:

  /space/freesurfer/build/scripts

are under CVS control:

  build_dev.csh
  build_stablepub.csh
  build_release_type.csh
  create_targz.csh
  exclude_from_targz
  postprocess_targz.csh
  package_maker.tar.gz
  build_scripts_README.txt (this file)

After making any changes to these files, copy them
to your CVS dev/distribution directory and commit them.

The package_maker directory (which contains a README
file of its own explaining the Mac OS package build
process) should be made into a tarball (named
package_maker.tar.gz) and also commited to dev/distribution.
The reason its a tarball is so that the symlinks are
stored, as these are important to the create_dmg script.

Per the files listed in the exclude_from_targz file, it is
important that the grad_unwarp_table directory is excluded from 
public distribution, as MGH is not allowed to distribute it.

#cvs log: '$Id: build_scripts_README.txt,v 1.1 2006/03/03 21:00:19 nicks Exp $'
