#! /usr/bin/env python

import sys,os,re
import commands

# debug = True
debug = False
whoami = sys.argv[0]

def main():

   if len(sys.argv) > 1:
      target_dir = sys.argv[1]
   else:
      target_dir = '.'
      # pass subdir thru CMakeLists.txt via CMAKE_INSTALL_PREFIX, e.g., "/Volumes/partition_2/freesurfer_src/git/install/Freeview.app/Contents"

   # Find all macho files in the application bundle (binary, libraries, and frameworks)
   # $ (cd <install_base>/Freeview.app/Contents && find . -type f -exec file {} \; | grep -i mach | sed 's;:.*;;' | tr -s '\n' ' ')

   find_macho_files_cmd = "( cd " + target_dir + " && find . -type f -exec file {} \; | grep -i mach | sed 's;:.*;;' | tr -s " + "'\\n'" + " " +  "' '" + ")"
   file_string = ""
   file_string = run_cmd(find_macho_files_cmd)

   # when counting slashes in each files path, assume path starts with leading ./ via find output
   # e.g., ./Frameworks/QtCore.framework/Versions/5/QtCore

   if file_string:
      # if (debug): print("files to process = %s") % file_string
      print("POSTPROCESSING OF MACDEPLOYQT OUTPUT by script %s") % os.path.abspath(whoami)
      print("files to process = %s") % file_string
      for file in file_string.split():
         slash_cnt = len(re.findall("/", file))
         dot_dots = slash_cnt - 1
         rel_path = ''
         for cnt in range(dot_dots,0,-1):
            rel_path = rel_path + '../'
         # Add a source relative RPATH entry to find Qt Frameworks as macdeployqt command may not do so when creating the application bundle.
         # If the RPATH already exists, then install_name_tool will do nothing, but will return status 1 for failure
         # $ (cd <install_base>/Freeview.app/Contents && sudo install_name_tool -add_rpath "@loader_path/../Frameworks" ./MacOS/Freeview)

         rpath_cmd = '( cd ' + target_dir + ' && sudo install_name_tool -add_rpath' + ' ' + '"' + '@loader_path/' + rel_path + 'Frameworks' + '"' + ' ' + file + ')'
         # rpath_cmd_output = run_cmd(rpath_cmd,True,False)
         rpath_cmd_output = run_cmd(rpath_cmd,False,True)
         if rpath_cmd_output:
            print("%s") % rpath_cmd_output


def run_cmd(cmd,fake=False,report_cmd=False):
   if not fake:
      if debug: print "Running cmd: %s" % (cmd)
      cmd_ret = commands.getstatusoutput(cmd)
      if cmd_ret[0] == 0:
         if (report_cmd): print "command succeeded: %s" % (cmd)
         return cmd_ret[1]
      else:
         if (report_cmd): print "*** command failed: %s" % (cmd)
         # sys.exit(1)
         return False
   else:
      print "+++ Would exec: %s" % (cmd)
      return True

if __name__ == "__main__":
   main()

