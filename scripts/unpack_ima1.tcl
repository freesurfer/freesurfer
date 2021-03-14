
##
## unpack_ima1.tcl
## Tthis cript looks at the headers of ima files in a predetermined 
## archive directory. The default archive directory is over-ridden 
## by the environment variable ARCHIVE_DIR. the user selects a 
## session and the path to that session is provided to the script 
## that copies the relevant files via nfs to the local machine and 
## unpacks them locally into b shorts.
##
##
## Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
##
## Terms and conditions for use, reproduction, distribution and contribution
## are found in the 'FreeSurfer Software License Agreement' contained
## in the file 'LICENSE' found in the FreeSurfer distribution, and here:
##
## https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
##
## Reporting: freesurfer@nmr.mgh.harvard.edu
##

#--------------------------------------------------------------------------------------#
#export ALPHA_BIN=/space/annecy/5/users/inverse/freesurfer_alpha/bin
#cp unpack_ima.tcl $ALPHA_BIN/noarch
#export PATH=$PATH:$ALPHA_BIN/noarch
#export PATH=$PATH:$ALPHA_BIN/Linux

#-----------------------------------------------------------------------------------------#
# global variable defined elsewhere

# session(sessionNumber, fieldValue) A 2-D array conataining all that is known about the sessions. 
#   The field Keys are:
#     0 patient_name
#     1 patient_ID
#     2 study_date
#     3 path to first file in session directory
#
# displayedSessions A list containing the "session" indices which are currently displayed 
#    in the list box


#-----------------------------------------------------------------------------------------#

#define needed variables

set targetDir [exec pwd ]
set archiveDir /space/sharbot/1/siemens
set sourceDir $archiveDir

set noArchBin /space/annecy/5/users/inverse/freesurfer_alpha/bin/noarch
set archBin "/space/annecy/5/users/inverse/freesurfer_alpha/bin/[exec uname]"

set mincOnly 0
set copyDone 1

set typeList { \
               {"Siemens file" ".ima" "SIEM" } \
               {"MINC file" ".mnc" "MINC" } \
               {"all" "*" {} } \
              }
set targetDirWritable 0

#-----------------------------------------------------------------------------------------#
# set help strings
set eh "Select a session that you want unpacked. Pick a directory to unpack it. Hit the unpack button. This process can take up to 20 minutes."
set fileBrowserHelp "It seems that this (library) tool only allows for the selection of files. If there's a file in your target directory, it's okay. Otherwise, you have to use the textbox at the bottom. They update each other, so you can browse right up to the directory you want, and finish in the text box. I prolly won't write my own file browser since this app is only a short-term solution."
set commandLineHelp "This script takes 0 or 1 arguments. The argument must be a writeable directory in which the minc files can be unpacked."
set viewProgressHelp "Creates a window to view the output of the unpacker while it is running."
set mincOnlyHelp "Transfers only minc files and does not unpack to b shorts."
set pathHelp "source /space/beowulf/1/users/inverse/freesurfer_alpha/bin/noarch/nmr-std-env to set up your environment."
set bugHelp "none reported 28 March 00."

source $noArchBin/progressbar.tcl


#--------------------------------------------------------------------------------------#
#-------------------------    check for dependencies     ------------------------------#
#--------------------------------------------------------------------------------------#

if { [catch { info exists env(ARCHIVE_DIR)} envArchiveDir ] } \
    { 
	if { [ file isdirectory $env(ARCHIVE_DIR) ] && [ file readable $env(ARCHIVE_DIR)] } \
	    { set archiveDir $env(ARCHIVE_DIR) }
    }


if { ! [file isdirectory $archiveDir] } \
   {
      tk_messageBox -type ok -default ok -title "Error" \
                           -message "could not find archive directory $archiveDir" \
                           -icon error
      exit
   }



catch {exec which mri_info} exitStatus
if { ! [file executable $exitStatus ] } \
 {
     if { [ file isdirectory $exitStatus ] }
    tk_messageBox -type ok -default ok -title "Error" \
                           -message $exitStatus \
                           -icon error
    exit
  }



#--------------------------------------------------------------------------------------#
#--------------------------------    subroutines     ----------------------------------#
#--------------------------------------------------------------------------------------#



proc Dialog_Create {top title args} \
     {
	global dialog
	 if { [winfo exists $top] } \
          {
		switch -- [wm state $top] \
                    {
			normal    {
				    # Raise a buried window
				      raise $top
			          }
			withdrawn -
			iconic    {
			 	    # Open and restore geometry
				     wm deiconify $top
				    catch { wm geometry $top $dialog(geo,$top) }
			          }
		    }

		return 0
	  } \
        else {
		eval {toplevel $top} $args
		wm title $top $title
		return 1
	}
}


proc Dialog_Wait {top varName {focus {}}} \
  {
	upvar $varName var

	# Poke the variable if the user nukes the window
	bind $top <Destroy> [list set $varName cancel]

	# Grab focus for the dialog
	if {[string length $focus] == 0} { set focus $top }
	set old [focus -displayof $top]
	focus $focus
	catch {tkwait visibility $top}
	catch {grab $top}

	# Wait for the dialog to complete
	tkwait variable $varName
	catch {grab release $top}
	focus $old
  }

#--------------------------------------------------------------------------------------#

proc Dialog_Dismiss {top} \
  {
	global dialog
	# Save current size and position
	catch {
		# window may have been deleted
		set dialog(geo,$top) [wm geometry $top]
		wm withdraw $top
	      }
        .commandFrame.viewProgressButton config -state normal
  }

#--------------------------------------------------------------------------------------#

proc CreateOutputMonitor {} \
    {
        global log dialog 
        set f .outputWindow

        if {[Dialog_Create $f "Progress" -borderwidth 10 ] } \
	   {
                             
              set f .outputWindow
              set tf [ frame $f.textBoxFrame ]

              set log [ text $tf.log -width 80 -height 20 -borderwidth 2 -relief raised \
                       -setgrid true ]
              scrollbar $tf.yScrollbar -orient vert -command "$tf.log yview"
              pack $tf.yScrollbar -fill y -side right
              pack $tf.log -side top -fill both -expand true

              set bf [ frame $f.buttonFrame ]
              button $bf.closeButton -text Close -command { DeleteOutputMonitor }
              button $bf.clearButton -text Clear -command { $log delete 1.0 end }
              pack $bf.clearButton -side left
	      pack $bf.closeButton 
       
              pack $tf -fill both -expand true
              pack $bf
              bind $f <Control-c> { DeleteOutputMonitor }
              
          }
        .commandFrame.viewProgressButton config -state disabled
    }

proc DeleteOutputMonitor {} \
    {
        Dialog_Dismiss .outputWindow
        .commandFrame.viewProgressButton config -state normal
    }


#--------------------------------------------------------------------------------------#

proc CreateProgressBarWindow { caption } \
    {
        global log dialog progressValue
        set pbw .myProgressBarWindow

        if {[Dialog_Create $pbw "Progress" -borderwidth 10 ] } \
	   {
                             
              #set pbw .myProgressBarWindow
              
              set pbWiget [::progressbar::progressbar $pbw.pb -variable progressValue]
              label $pbw.progressLabel -text $caption

              pack $pbw.progressLabel
              pack $pbWiget
              
              bind $pbw <Control-c> { DestroyProgressBarWindow }
              
          } \
	 else \
           {
             $pbw.progressLabel configure -text $caption
	   }

        
    }

proc DestroyProgressBarWindow {}   { destroy .myProgressBarWindow }


#--------------------------------------------------------------------------------------#

proc CreateSearchDialog {} \
    {
        global  dialog searchString searchDialogClosed
        set sd .searchDialog

        if {[Dialog_Create $sd "Search" -borderwidth 10 ] } \
	   {
                             
              set sd .searchDialog
              set nameFrame [ frame $sd.nameBoxFrame ]

              label $nameFrame.searchEntryBoxlabel -text "Patient name"
              entry $nameFrame.nameEntryBox -textvariable searchString -relief sunken 

              set spacer1 [ frame $sd.spacer1 -height 20 ]

              set buttonFrame [ frame $sd.buttonFrame ]
              button $buttonFrame.cancelButton -text Cancel -command { set searchString ""
                                                              DestroySearchDialog
                                                            }
              button $buttonFrame.submitButton -text Submit -command { DestroySearchDialog }

    #--- pack ----#
              pack $nameFrame.searchEntryBoxlabel $nameFrame.nameEntryBox   -side left
              pack $buttonFrame.cancelButton      $buttonFrame.submitButton -side left
              
              pack $nameFrame -fill both -expand true
              pack $spacer1
              pack $buttonFrame

   #--- bind ----#
              bind $nameFrame.nameEntryBox <Return> { DestroySearchDialog }
	      bind $sd <Control-c>       { set searchString ""; DestroySearchDialog }
              
          }

    }

proc DestroySearchDialog {} \
  { 
    global searchDialogClosed 
    set searchDialogClosed 1
    destroy .searchDialog 
  }
 

#------------------------------------------------------------------------------------#

proc CreateFileBrowser { dirType callersDir } \
    {
        global dialog newDirSelected archiveDir targetDir fbInfo
        set fbInfo(currentDir) $callersDir
        set fbInfo(currentFile) $fbInfo(currentDir)/.

        #puts "currentFile $fbInfo(currentFile)"
        #puts "currentDir $fbInfo(currentDir)" 

        set fb .fileBrowser

        if {[Dialog_Create $fb "Select $dirType Directory" -borderwidth 10 ] } \
	   {
                             
    #--- wigets ----#
            set fb .fileBrowser
            set dirFrame [ frame $fb.dirFrame ]

            scrollbar $dirFrame.leftScroll -command "$dirFrame.fileList yview"
            set dirPick [listbox $dirFrame.fileList \
                  -yscroll "$dirFrame.leftScroll set" -font fixed \
		  -relief sunken -width 60 -height 15 -setgrid yes ]


            set spacer1 [ frame $fb.spacer1 -height 15 ]

            label $fb.dirEntryBoxlabel -text "$dirType Directory"
            entry $fb.dirEntryBox -textvariable fbInfo(currentDir) -relief sunken 

            set spacer2 [ frame $fb.spacer2 -height 5 ]

            set buttonFrame [ frame $fb.buttonFrame ]
            button $buttonFrame.cancelButton -text Cancel \
                                             -command { set newDirSelected 2 }

            button $buttonFrame.submitButton -text ok \
                                             -command \
                                { if {[ file isdirectory $fbInfo(currentDir) ]} \
                                     {
                                       set newDirSelected 1
                                     } \
				  else {
                                       tk_messageBox -type ok -default ok -title "Error" \
                                         -message "$fbInfo(currentDir) not a directory" \
                                         -icon error
				       }
				}
    #--- pack ----#
            pack $dirFrame.fileList -side left -fill both -expand true 
            pack $dirFrame.leftScroll  -side left -fill y -expand true 

            pack $buttonFrame.cancelButton  -side left
            pack $buttonFrame.submitButton -side left

            pack $dirFrame -side top
            pack $spacer1
            pack $fb.dirEntryBoxlabel -anchor center
            pack $fb.dirEntryBox -fill x -expand true
            pack $spacer2
            pack $buttonFrame
            

   #--- bind ----#

          bind $dirPick <ButtonRelease-1> \
            { 
              set fileIndex [  .fileBrowser.dirFrame.fileList curselection ]
              set subDir [ .fileBrowser.dirFrame.fileList get $fileIndex ]
		if { [string match ".." $subDir ]} \
		     { set fbInfo(currentFile) [ file dirname $fbInfo(currentFile) ] } \
		else { set fbInfo(currentFile) $fbInfo(currentDir)/$subDir }

		if {[ file isdirectory $fbInfo(currentFile) ] && [ file readable $fbInfo(currentFile) ] } \
                 {
		   set fbInfo(currentDir) $fbInfo(currentFile)
                   .fileBrowser.dirFrame.fileList del 0 end

		     if {[ file isdirectory "$fbInfo(currentDir)/.." ] } \
		       { .fileBrowser.dirFrame.fileList insert end ".." }
 
		     foreach listEntry [split [ exec /bin/ls -1 $fbInfo(currentDir) ] \n] \
                    {
                      .fileBrowser.dirFrame.fileList insert end $listEntry
                    }
                 }
            }

	 bind $fb.dirEntryBox <Return> {

             #set fileName [ .fileBrowser.dirEntryBox get ]
             #set fbInfo(currentFile) $fbInfo(currentDir)

	     if {[ file isdirectory $fbInfo(currentDir) ]} \
                 {
                   #set fbInfo(currentDir) $fbInfo(currentFile)

                   .fileBrowser.dirFrame.fileList del 0 end
 
		   foreach listEntry [split [ exec /bin/ls -1 $fbInfo(currentDir) ] \n] \
                    {
                      .fileBrowser.dirFrame.fileList insert end $listEntry
                    }
                 } 

            }

	    bind $fb <Control-c>   { fbInfo(currentDir) $archiveDir; set newDirSelected 2 }
	    bind $fb <Alt-q>       { fbInfo(currentDir) $archiveDir; set newDirSelected 2 }


   #--- main ----#

            if { ! [ file exists "$fbInfo(currentDir)" ] && [string match "Destination" $dirType] } \
		{
		    if { [ file writable [ file dirname $fbInfo(currentDir) ] ] } \
		          { 
                            file mkdir $fbInfo(currentDir)
			  } else { tk_messageBox -type ok \
                                             -title "Error" \
                                             -message "Overwrite contents of $targetDir" \
				             -icon error
                                   set fbInfo(currentDir) $env(HOME)
			         }
		}

	    if {[ file isdirectory "$fbInfo(currentDir)/.." ] } \
		       { .fileBrowser.dirFrame.fileList insert end ".." } \
            else \
                {
                   
                }

	    foreach listEntry [split [ exec /bin/ls -1 $fbInfo(currentDir) ] \n ] \
               {
                 #puts "$listEntry *"
                 .fileBrowser.dirFrame.fileList insert end $listEntry
               }
         
	}

    }

proc DestroyFileBrowser {} \
  { 
    destroy .fileBrowser 
  }

#--------------------------------------------------------------------------------------#

proc Log {pipeName} \
    {
	global  pipeState log 
        if [eof $pipeName] { KillPipe $pipeName } \
        else \
	    {
              gets $pipeName line
		if { $line == 0 } { KillPipe $pipeName }
              $log insert end $line\n
              $log see end
            }
    }


#--------------------------------------------------------------------------------------#

proc KillPipe {pipeName} \
    {
        global pipeState
        catch {close $pipeName}
	set pipeState 1 
    }


#----------------------------------   ReadIndexFile   -----------------------------------#
# reads a text file that has session tab-delimited description info and path

proc ReadIndexFile { indexFile } \
    {
        global session
        set x 0

	if { [file readable $indexFile ] } \
	  {
            set indexFileHandle [ open $indexFile r]
            foreach line [split [read $indexFileHandle ] \n ] \
	      {
		  set fields [ split $line \t ]

	          for {set i 0} {$i < 4} {incr i} \
	             {
	      	        set session($x,$i) [lindex $fields $i]
                        #puts "\$session($x,$i)=$session($x,$i)"
	             }
		  #regsub /local_mount/ $session($x,3) \/ tempPath
                  #set session($x,3) [ file dirname $tempPath ]
	          #puts "$x $session($x,3)"
                  incr x
	      }
           }
	return [expr $x-1]

    }


#--------------------------------  GetHeaderInfo  --------------------------------#


proc GetHeaderInfo { myIMAfile} \
    {
       global headerInfo archBin
       set myIndex 0
	if { [ catch {open "|$archBin/mri_info $myIMAfile" r} IN ] } \
           {puts "$IN"; return 1}
       #set IN [ open "|mri_info $myIMAfile" r ]
	set headerLines [ split [ read $IN ] \n ]
	close $IN

        foreach headerLine $headerLines \
          {
             set rawHeaderIndex [ lindex [ split $headerLine : ] 0 ]
	     regsub -all " " $rawHeaderIndex "_" headerIndex
             set rawHeaderValue [ lindex [ split $headerLine : ] 1 ] 
	     regsub -all {^[" "]} $rawHeaderValue "" headerValue
	     #puts "$headerIndex:$headerValue"
	     set headerInfo($headerIndex) $headerValue
       
          }
        return 1

    }


#--------------------------   GetSessionInfo   --------------------------------#

proc GetSessionInfoFromFile { count } \
    {
	global headerInfo session

        set headerInfo(patient_name) "unknown"
        set headerInfo(patient_id)   "unknown"
	set headerInfo(study_date)   "unknown"

	#puts "Getting $session($count,3)"
        set exitStatus [ GetHeaderInfo $session($count,3) ]

	set session($count,0) $headerInfo(patient_name)
	set session($count,1) $headerInfo(patient_id)
	set session($count,2) $headerInfo(study_date)

        update
	
         return $exitStatus

    }




#--------------------------   GetSessionDescriptors   --------------------------------#
# returns a formatted string of info suitable for displaying in list box

proc GetSessionDescriptors { i } \
    {
       global session numberOfTotalSessions
        
        set IDpadString ""
        set namePadString ""

	if { $i > $numberOfTotalSessions } \
           { set i [ expr $numberOfTotalSessions - 1 ] }
        

        #patient name
	set padSize [expr 30 - [ string length $session($i,0) ] ]
        if { $padSize < 0} { set padSize 0}
        while {$padSize} { append namePadString " "; incr padSize -1 }

        #patient id
        set padSize [expr 10 - [ string length $session($i,1) ] ]
        if { $padSize < 0} { set padSize 0}
        while {$padSize} { append IDpadString " "; incr padSize -1 }

        return "$session($i,0) $namePadString $session($i,1) $IDpadString $session($i,2)"

    }

#--------------------------------    LoadListBox    --------- ---------------------------#
# reads the index file, the incoming directory, and writes info to list box

proc LoadListBox {} \
    {
      global sessionInfoFrame progressValue \
             indexFile archiveDir displayedSessions \
             sessionDescriptor sampleSessionFiles mySessionDescriptor mySampleSessionFiles \
             session numberOfTotalSessions numberOfIndexedSessions numberOfNewDirs


	#-------- LOAD LISTBOX FROM INDEX ---------#

if { [file exists $indexFile ] } \
   { 
      set numberOfIndexedSessions [ReadIndexFile $indexFile]
      set numberOfTotalSessions $numberOfIndexedSessions
      #puts "numberOfIndexedSessions $numberOfIndexedSessions"


      #load index into list box
      for {set index 0} { $index < $numberOfIndexedSessions } {incr index} \
          {
	      $sessionInfoFrame.sessionList insert end [ GetSessionDescriptors $index ]
          }
    }

update

           #------------ GET DIRECTORY LISTING ---------#

set fileList [ exec ls $archiveDir ] 
#puts $fileList

set index 0
foreach fileName $fileList \
   { 
      
      if { [ file isdirectory [ file join $archiveDir $fileName ] ] } \
         { 
           set dirList($index) [ file join $archiveDir $fileName ]
           #puts $dirList($index)
           incr index
         }
   }

set numberOfNewDirs $index; #puts "numberOfDirs $numberOfDirs"
set numberOfTotalSessions [expr $numberOfNewDirs + $numberOfIndexedSessions]


              #------- GET HEADER INFO FOR NEW SESSIONS --------#

set index $numberOfIndexedSessions

if { [ set savedCursor [lindex [$sessionInfoFrame configure -cursor] 4]] != "watch"} \
            { $sessionInfoFrame configure -cursor watch }

CreateProgressBarWindow "Reading incoming directory"


for {set dirNumber 0} { $dirNumber < $numberOfNewDirs } {incr dirNumber} \
  {
      

    set dirContents [ exec ls $dirList($dirNumber) ]
     foreach bareFileName $dirContents   \
      {
	set fileName [ file join $dirList($dirNumber) $bareFileName ]

        #get the first file as a sample
	if { [ file isfile $fileName ] } \
           { 
             #puts $fileName
             set session($index,3) $fileName
             GetSessionInfoFromFile $index
             
             incr index
             break
           }
       } 
    set progressValue [expr $dirNumber*100/$numberOfNewDirs]
    #puts $progressValue
    update
     
  }
  DestroyProgressBarWindow
  $sessionInfoFrame configure -cursor $savedCursor

# enable Option menue items
  .menubar.mOption entryconfigure 1 -state normal
  .menubar.mOption entryconfigure 2 -state normal
  .menubar.mOption entryconfigure 3 -state normal

  update


              #-----------   LOAD LISTBOX   ---------------#

for {set index $numberOfIndexedSessions} { $index < $numberOfTotalSessions } {incr index} \
    {
      #create master listing
#      set sessionDescriptor($index) [ GetSessionDescriptors $index ]

      #push entry onto listbox
      $sessionInfoFrame.sessionList insert end [ GetSessionDescriptors $index ]
    }

 update

           
# the list of displayed sessions contains the entire master array of sessions
 for {set index 0} { $index < $numberOfTotalSessions } {incr index} \
     {
       lappend displayedSessions $index
     }

 update
}


#--------------------------------------------------------------------------------------#

proc CheckDirOK {} \
    {
        global targetDir sourceDir archiveDir targetDirWritable
	set targetDirReady 0
	set sourceDirReady 0

        if { ! $targetDirWritable } \
	    {

        if [ file exists $targetDir ] \
	    {
		if { [ string compare [ exec ls $targetDir ] ""]} \
		    {
			set overwriteAnswer [tk_messageBox -type yesno \
                                             -default yes \
                                             -title "$targetDir not empty" \
                                             -message "Overwrite contents of $targetDir" \
				             -icon question ]
                        if { ![string compare $overwriteAnswer no]} \
                           {
                             set targetDirWritable 0
                             return 0}
                        
		    }

	    }


	if { ! [ file exists $targetDir ]} \
           { 
              if { [ file writable [ file dirname $targetDir ] ] } \
                 { file mkdir $targetDir } \
              else \
	       {
                  tk_messageBox -type ok -default ok -title "Error" \
                       -message "could not create $targetDir" -icon error 
                  set targetDir .
                  return 0
	       }
	   }

	if { ![ file isdirectory $targetDir ] } \
	       {
                  tk_messageBox -type ok -default ok -title "Error" \
                       -message "$targetDir is not a directory" -icon error 
                  set targetDir .
                  return 0
	       }


	if { ! [ file writable $targetDir ] } \
              {
                  tk_messageBox -type ok -default ok -title "Error" \
                       -message "$targetDir is not writeable" -icon error 
                  set targetDir .
                  return 0
	       }
   set targetDirWritable 1

  }

	if { ! [ file isdirectory $sourceDir ] }\
              {
                  tk_messageBox -type ok -default ok -title "Error" \
                       -message "$sourceDir is not a directory" -icon error 
                  set sourceDir $archiveDir
                  return 0
	       }

 	if { ![ file readable $sourceDir ] } \
              {
                  tk_messageBox -type ok -default ok -title "Error" \
                       -message "$sourceDir is not readable" -icon error 
                  set sourceDir $archiveDir
                  return 0
	       }

    return 1

  }

#----------------------------------------------------------------------------#

proc TransferIMAfiles {destinationDir} \
    {
    global sourceDir copyDone log progressValue
    set dirContents [ exec ls $sourceDir ]
    set numberOfFilesTransferred 0
    set dirLength [ llength $dirContents ]
    set progressValue 0

    CreateProgressBarWindow "Transferring IMA files"
    puts "dirLength $dirLength"

    foreach fileName $dirContents \
      {
	#puts "cp ${sourceDir}/${fileName} ${destinationDir}/${fileName}"
        if { $copyDone } \
           { 
             $log insert end "Transfer aborted\n$numberOfFilesTransferred files transferred\n"
             $log see end
             DestroyProgressBarWindow
             update
             return 1 
           }

        if {[catch {exec cp ${sourceDir}/${fileName} ${destinationDir}/${fileName} } \
			 errorMsg ]} {puts stderr $errorMsg}
        incr numberOfFilesTransferred       
	$log insert end "${destinationDir}/${fileName}\n"
        $log see end

        set progressValue [ expr $numberOfFilesTransferred*100/$dirLength ]
        update
      }

    $log insert end "\n$numberOfFilesTransferred files transferred\nTransfer complete\n"
    $log see end
    update
    #set sourceDir $destinationDir
     DestroyProgressBarWindow
    
    }


#----------------------------------------------------------------------------#
proc TransferMINCfiles {destinationDir} \
    {

    global sourceDir copyDone log archBin progressValue
    set dirContents [ exec ls $sourceDir ]
    set numberOfFilesTransferred 0
    set dirLength [ llength $dirContents ]
    set progressValue 0

    CreateProgressBarWindow "Transferring MINC files"
    puts "dirLength $dirLength"

    foreach fileName $dirContents \
      {
        set imaFilePath  ${sourceDir}/${fileName}
        set mincFilePath ${destinationDir}/[ file rootname $fileName ].mnc

	#puts "cp ${sourceDir}/${fileName} ${destinationDir}/${fileName}"
        if { $copyDone } \
           { 
             $log insert end "Transfer aborted\n$numberOfFilesTransferred files transferred\n"
             $log see end
             DestroyProgressBarWindow
             update
             return 1 
           }

        if {[catch {exec $archBin/mri_convert $imaFilePath $mincFilePath }  errorMsg ]} \
           {
             # error!
             puts stderr $errorMsg
             $log insert end "failed: $mincFilePath\n"
             $log insert end $errorMsg
             $log see end
           } \
        else \
	    {
             # successful copy
               incr numberOfFilesTransferred       
	       $log insert end "$mincFilePath\n"
               $log see end
	    }
        set progressValue [ expr $numberOfFilesTransferred*100/$dirLength ]
        update
      }

    $log insert end "\n$numberOfFilesTransferred files transferred\nTransfer complete\n"
    $log see end
    update
    #set sourceDir $destinationDir
     DestroyProgressBarWindow
    
    }


#----------------------------------------------------------------------------#
proc TransferBshortFiles {destinationDir} \
    {

    global sourceDir copyDone log archBin progressValue
    set dirContents [ exec ls $sourceDir ]
    set numberOfFilesTransferred 0
    set dirLength [ llength $dirContents ]
    set progressValue 0

    CreateProgressBarWindow "Transferring FS-FAST files"
    puts "dirLength $dirLength"

    foreach fileName $dirContents \
      {
        set imaFilePath  ${sourceDir}/${fileName}
        set mincFilePath ${destinationDir}/[ file rootname $fileName ].bshort

	#puts "cp ${sourceDir}/${fileName} ${destinationDir}/${fileName}"
        if { $copyDone } \
           { 
             $log insert end "Transfer aborted\n$numberOfFilesTransferred files transferred\n"
             $log see end
             DestroyProgressBarWindow
             update
             return 1 
           }

        if {[catch {exec $archBin/mri_convert $imaFilePath $mincFilePath }  errorMsg ]} \
           {
             # error!
             puts stderr $errorMsg
             $log insert end "failed: $mincFilePath\n"
             $log insert end $errorMsg
             $log see end
           } \
        else \
	    {
             # successful copy
               incr numberOfFilesTransferred       
	       $log insert end "$mincFilePath\n"
               $log see end
	    }
        set progressValue [ expr $numberOfFilesTransferred*100/$dirLength ]
        update
      }

    $log insert end "\n$numberOfFilesTransferred files transferred\nTransfer complete\n"
    $log see end
    update
    #set sourceDir $destinationDir
     DestroyProgressBarWindow
    
    }

#----------------------------------------------------------------------------#

proc CreateMINCfiles {destinationDir} \
    {

        global archBin 
	set imaDir "/tmp/ima[pid]"	
      file mkdir $imaDir
      TransferIMAfiles $imaDir
      set dirContents [ exec ls $imaDir ]
      foreach fileName $dirContents \
	 {
            puts "mri_convert ${imaDir}/${fileName} ${destinationDir}/${fileName}"
	     if {[catch {exec $archBin/mri_convert ${imaDir}/${fileName} ${destinationDir}/[file rootname $fileName].mnc } errorMsg ]} {puts stderr $errorMsg}
	  }
      puts "IMA -> MINC done\n"
      #file delete -force $imaDir
    }


#----------------------------------------------------------------------------#

proc ConvertMINCfiles {mincOnly} \
    {
      global targetDir
      set mincDir "/tmp/minc[pid]"	
      file mkdir $mincDir
      ConvertIMAfiles $mincDir
      set unpackmincdirPipe [open "|unpackmincdir  -src $mincDir -targ $targetDir $mincOnly" ]

      fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
      fconfigure $unpackmincdirPipe -blocking 0
      vwait pipeState
      
   }

#-----------------------   SEARCH DESCRIPTORS FOR MATCH   --------------------------#

proc SearchDescriptors {} \
    {
       global  sessionInfoFrame session searchString displayedSessions

       #puts $searchString
       #puts $displayedSessions
       foreach listIndex $displayedSessions \
	 {
           
           if { [ string match "*$searchString*" $session($listIndex,0) ] } \
	       {
		  lappend hitList $listIndex
                  #puts "$session($listIndex,0)"
	       }
         }

	if { ! [ info exists hitList ] } \
           { 
              tk_messageBox -type ok \
                            -default ok -title "debug" \
                            -message "no match found for \"$searchString\""
              return 1
           }

       # reset displayed list
       unset displayedSessions
       set displayedSessions $hitList
       #puts $displayedSessions
       $sessionInfoFrame.sessionList delete 0 end

       foreach listIndex $displayedSessions \
	 {
           $sessionInfoFrame.sessionList insert end [ GetSessionDescriptors $listIndex ]
           #puts [ GetSessionDescriptors $listIndex ]
	 }

       update

    }


#=============================================================================#
#---------------------------- Define widgets   -------------------------------#
#=============================================================================#


menu .menubar
#attach it to the main window
. config -menu .menubar
foreach menuWidget { File Option Help } \
    {
	set $menuWidget [ menu .menubar.m$menuWidget ]
        .menubar add cascade -label $menuWidget -menu .menubar.m$menuWidget 
    }

$File add command -label "Change archives" -command \
                       { 
                        set newDirSelected 0
                        .menubar.mFile entryconfigure 1 -state disabled
                        .menubar.mFile entryconfigure 2 -state disabled

			CreateFileBrowser "Archive" $archiveDir
			tkwait variable newDirSelected

			   if { $newDirSelected == 1 } \
			      { set archiveDir $fbInfo(currentDir) }

	               if { [CheckDirOK] } { $commandFrame.startButton configure -state normal } \
                        else { $commandFrame.startButton configure -state disabled }


                         DestroyFileBrowser
                        .menubar.mFile entryconfigure 1 -state normal
                        .menubar.mFile entryconfigure 2 -state normal
                       }

$File add command -label "New destination dir" -command \
                       { 
                        set newDirSelected 0
                        .menubar.mFile entryconfigure 1 -state disabled
                        .menubar.mFile entryconfigure 2 -state disabled

			CreateFileBrowser "Destination" $targetDir
			tkwait variable newDirSelected

			   if { $newDirSelected == 1 } \
			      { set targetDir $fbInfo(currentDir) }

	               if { [CheckDirOK] } { $commandFrame.startButton configure -state normal } \
                        else { $commandFrame.startButton configure -state disabled }

                         DestroyFileBrowser
                        .menubar.mFile entryconfigure 1 -state normal
                        .menubar.mFile entryconfigure 2 -state normal
                       }

$File add command -label "View File"  -command \
                       { set fileName [ tk_getOpenFile  -title "Select ima file" \
                                                        -filetypes $typeList \
                                                        -initialdir $archiveDir ]
			
			   if {  [file readable $fileName ] } \
			     { 
			    if [ catch {exec $noArchBin/tkmedit_wrapper $fileName } exitStatus ] \
                                   {  tk_messageBox -type ok \
                                              -default ok -title "Error" \
					      -message $exitStatus
				   }
		             }
		       }

$File add command -label Quit -command {  exit   }


$Option add command -label "Refresh"   -state disabled \
                                       -command { 
                                                  set copyDone 1
                                                  $sessionInfoFrame.sessionList delete 0 end
                                                  LoadListBox 
                                                }

$Option add command -label "Search"  -state disabled \
                                     -command {
                                                set searchDialogClosed 0
                                                CreateSearchDialog
                                                vwait searchDialogClosed
                                                if { [string compare $searchString "" ] } \
				                    {
                                                       #puts $searchString
                                                       SearchDescriptors
                                                    }
				              }

$Option add command -label "debug info" -state disabled \
                                        -command { tk_messageBox -type ok \
                                              -default ok -title "debug" \
					      -message "indexed sessions: $numberOfIndexedSessions\nnew sessions: $numberOfNewDirs\nArchive Dir: $archiveDir\nDest Dir: $targetDir" }

$Help add command -label "Eh?" -command { tk_messageBox -type ok \
                                              -default ok -title "Eh?" \
					      -icon question \
                                              -message $eh }


$Help add command -label "command line" -command { tk_messageBox -type ok \
                                              -default ok -title "command line" \
					      -icon question \
                                              -message $commandLineHelp }

$Help add command -label "view progress" -command { tk_messageBox -type ok \
                                              -default ok -title "view progress" \
					      -icon question \
                                              -message $viewProgressHelp }

$Help add command -label "minc only" -command { tk_messageBox -type ok \
                                              -default ok -title "minc only" \
					      -icon question \
                                              -message $mincOnlyHelp }

$Help add command -label "setting paths" -command { tk_messageBox -type ok \
                                              -default ok -title "required paths" \
					      -icon question \
                                              -message $pathHelp }

$Help add command -label "bugs" -command { tk_messageBox -type ok \
                                              -default ok -title "bugs" \
					      -icon question \
                                              -message $bugHelp }

#-----------------------   Session  --------------------------------#

frame .leftFrame

set sessionInfoFrame [ frame .leftFrame.sessionInfoFrame -borderwidth 1  ]

scrollbar $sessionInfoFrame.leftScroll -command "$sessionInfoFrame.sessionList yview"
set sessionPick [listbox $sessionInfoFrame.sessionList \
                  -yscroll "$sessionInfoFrame.leftScroll set" -font fixed \
		  -relief sunken -width 60 -height 15 -setgrid yes ]


#  pack sessionFrame  #
pack $sessionInfoFrame.sessionList -side left -fill both -expand true 
pack $sessionInfoFrame.leftScroll  -side left -fill y -expand true 


label .leftFrame.targetDirEntryBoxlabel -text "Destination Directory"
entry .leftFrame.targetDirEntryBox -textvariable targetDir -relief sunken 

#  pack leftframe  #
pack $sessionInfoFrame
pack .leftFrame.targetDirEntryBoxlabel -anchor center
pack .leftFrame.targetDirEntryBox -side top -fill x -expand true


#-----------------------------------------------------------------------#

set commandFrame [ frame .commandFrame ]

set commandFrameButtonWidth 8
button $commandFrame.stopButton -text "stop" \
                                -width $commandFrameButtonWidth \
                                -state disabled \
                                -command \
                                {
				  set copyDone 1
				  if { [ info exists unpackmincdirPipe] } \
                                       { KillPipe $unpackmincdirPipe }
                                 
                                 $commandFrame.stopButton config -state disabled
                                 $commandFrame.startButton config -state normal

                                 }

button $commandFrame.viewProgressButton -text "view log" \
                                -width $commandFrameButtonWidth \
                                -command CreateOutputMonitor

frame       $commandFrame.spacer1 -height 15

set tranferTypeFrame [ frame $commandFrame.tranferTypeFrame -borderwidth 2 -relief ridge  ]

label       $tranferTypeFrame.transferTypeLabel -text "Transfer type"
radiobutton $tranferTypeFrame.bshortsRadioButton -text "b shorts" -value bshorts -variable transferType
radiobutton $tranferTypeFrame.mincRadioButton -text "minc" -value minc -variable transferType
radiobutton $tranferTypeFrame.imaRadioButton -text "ima"   -value ima  -variable transferType

$tranferTypeFrame.imaRadioButton select

pack $tranferTypeFrame.transferTypeLabel -pady 5 -padx 5 -anchor c
pack $tranferTypeFrame.bshortsRadioButton  -anchor w
pack $tranferTypeFrame.mincRadioButton     -anchor w
pack $tranferTypeFrame.imaRadioButton      -anchor w


button $commandFrame.startButton -text "start" \
                                -state disabled \
                                -width $commandFrameButtonWidth \
                                -command \
     {   
        set pipeState 0
        set copyDone 0
        $commandFrame.startButton config -state disabled
        $commandFrame.stopButton config -state normal

        switch -exact $transferType \
	 {
	     minc    { TransferMINCfiles $targetDir }

             bshorts { TransferBshortFiles $targetDir }

             ima     { TransferIMAfiles $targetDir }

             default {}
	 }
	   
         
         $commandFrame.stopButton config -state disabled
         $commandFrame.startButton config -state normal 
      }

pack $commandFrame.startButton
pack $commandFrame.stopButton 
pack $commandFrame.viewProgressButton
pack $commandFrame.spacer1
pack $tranferTypeFrame -ipadx 5 -ipady 5


#====================================================#
#-------------------    PACK   ----------------------#
#====================================================#


pack .leftFrame -side left -fill both -expand true
pack .commandFrame -side right -anchor n
update

#=====================================================================#
#-------------------------    BINDINGS  ------------------------------#
#=====================================================================#


bind $sessionPick <ButtonRelease-1> \
  { 
    set sessionIndex [ $sessionInfoFrame.sessionList curselection ]
    incr sessionIndex
      set sourceDir  [ file dirname $session([expr $sessionIndex-1],3) ]
    #puts $sourceDir
    if { [CheckDirOK] } { $commandFrame.startButton configure -state normal }
  }

bind .leftFrame.targetDirEntryBox <Return> \
    {
	if { [CheckDirOK] } { $commandFrame.startButton configure -state normal } \
         else { $commandFrame.startButton configure -state disabled }
    }

bind all <Control-c> {destroy .}
bind all <Alt-q> {destroy .}

#=====================================================================#
#-----------------------------    MAIN   -----------------------------#
#=====================================================================#

#----------------------- COMMAND LINE STUFF ---------------------------------#

if { $argc > 1 } \
  {
    tk_messageBox -type ok -default ok -title "Error" \
                       -message "usage:  unpack_ima.tcl \[distination dir\]" -icon warning
    exit
  }

switch -exact $argc \
    {
      0  { set targetDir . }

      1  {
            if { [file isdirectory $targetDir] && [ file writable  $targetDir] } \
               { set targetDir [ lindex $argv 0 ] }
	 }
 
      default  { 
                  tk_messageBox -type ok -default ok -title "Command Line Help" \
	  			 	                 -icon question \
                                                         -message $commandLineHelp
                set targetDir .
	       }


    }

#----------------------- DEPENDENCIES STUFF ---------------------------------#

if {! [file readable $archiveDir ] } \
    { 
       tk_messageBox -type ok -default ok -title "Command Line Help" \
					      -icon error \
                                              -message "Could not find archive directory $archiveDir"
       exit
    }



#----------------------- READ MASTER INDEX ---------------------------------#

set indexFile  [file join [file dirname $archiveDir] index.txt ]

#puts "indexFile: $indexFile"

# set numberOfIndexedSessions [ ReadIndexFile $indexFile ]
# puts "numberOfIndexedSessions: $numberOfIndexedSessions\n"

# set newDirs [exec ls $archiveDir ]
# set numberOfNewDirs [llength newDirs]
# set numberOfTotalSessions [ expr $numberOfIndexedSessions + $numberOfNewDirs ]

# puts "numberOfNewDirs: $numberOfNewDirs\n"
# puts "numberOfTotalSessions: $numberOfTotalSessions\n\n"

#-----------------------     LOAD LISTBOX    ---------------------------------#

# print the contents of "session" array to listbox
# for { set listBoxLine 0} {$listBoxLine < $numberOfTotalSessions} {incr listBoxLine} \
#     {
#       #puts $listBoxLine
#       set namePadString ""
#       set IDpadString ""
#       set padSize [expr 25 - [ string length $session($listBoxLine,0) ] ]
      
#       while {$padSize} { append namePadString " "; incr padSize -1 }
      
#       set padSize [expr 15 - [ string length $session($listBoxLine,1) ] ]
      
#       while {$padSize} { append IDpadString " "; incr padSize -1 }
      

#       #name ID date path
#       set listBoxString "$session($listBoxLine,0)$namePadString $session($listBoxLine,1)$IDpadString   $session($listBoxLine,2)"
#       #puts "$listBoxString"
#       $sessionInfoFrame.sessionList insert end $listBoxString
#     }

LoadListBox

CreateOutputMonitor
Dialog_Dismiss .outputWindow

#puts "end of program"


#========================================================================#
#------------------------------    TO DO   ------------------------------#
#========================================================================#
set temp " \

"
