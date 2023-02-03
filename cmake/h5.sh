# Allow large *.h5 model files to live in the annex.
# If *.h5 files exist as soft links to annex files, then they need to be unlocked before or during the install pass.
h5_files=$(find . -maxdepth 1 -name "*.h5")
if [ "$h5_files" == "" ]; then
   # echo "Nothing to do"
   exit 0
fi
for file in $h5_files
do
   if [ -L "$file" ]; then
     # echo "should unlock $file"
     git annex unlock $file
   fi
done
# rm -f .models_unlocked && touch .models_unlocked
