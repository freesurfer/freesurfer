#!/usr/bin/env sh

if [ $# -lt 2 ]; then
   echo "Must have at least 2 arguments (<install path> <1 or more binary names>). Exitting..."
   exit 1
fi

# 1st arg = install path
install_path=$1
if [ ! -d "$install_path" ]; then
   echo "Cannot stat install path $install_path. Exitting..."
   exit 1
fi
bin_path=$install_path/bin
mkdir -p $bin_path

# following args are discontinued commands under <install path>/bin
# generate script to run instead that outputs appropriate message
shift 
for binary do
   fpath=$bin_path/$binary
   rm -f $fpath && touch $fpath
   echo "#!/usr/bin/env sh" >> $fpath
   echo "" >> $fpath
   echo "echo \"Please note that the command $binary is no longer built by default.\"" >> $fpath
   echo "echo \"You can inquire to the Freesurfer help list for further assistance, freesurfer@nmr.mgh.harvard.edu\"" >> $fpath
   chmod 755 $fpath
   if [ -e $fpath ]; then 
      echo "Created dummy command $fpath"
   else
      echo "*** Error: Could not create dummy file $fpath"
   fi
   shift
done

