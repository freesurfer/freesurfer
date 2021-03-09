
##
## unpack_mnc.tcl
## This script looks at the headers of minc files in a predetermined 
## archive directory. The default archive directory is over-ridden by
## the environment variable ARCHIVE_DIR. The scripts attempts to read
## an index file, "index.txt", which is a sister to the archive 
## directory. This index is produced by "mkmnc_index.tcl" and contains
## header info and paths to sessions in other directries. The user 
## selects a session.The path to that session is provided to the script
## "unpackmincdir," which copies the relevant files via nfs to the local
## machine and unpacks them locally into b shorts.
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
#cp unpack_mnc.tcl $ALPHA_BIN/noarch
#export PATH=$PATH:$ALPHA_BIN/noarch
#export PATH=$PATH:$ALPHA_BIN/Linux
#wrapper program is "browse-sessions"

#define needed variables

if { [ catch "exec pwd" workingDir]} {set targetDir $env(HOME)} else { set targetDir $workingDir}

set archiveDir /space/sharbot/2/minc
set sourceDir /space/sharbot/2/minc
set indexFile  /space/sharbot/2/index.txt
#set archiveDir /space/geneve/1/minctest
#set archiveDir /space/geneve/1/siemenstest
set noArchBin /space/annecy/5/users/inverse/freesurfer_alpha/bin/noarch
set archBin "/space/annecy/5/users/inverse/freesurfer_alpha/bin/[exec uname]"

set pipeState 0
set copyDone 0
set searchString "."
set targetDirWritable 0
set numberOfIndexedSessions 0
set mincOnly 0
set typeList { \
               {"MINC file" ".mnc" "MINC" } \
               {"Siemens file" ".ima" "SIEM" } \
               {"all" "*" {} } \
              }

set eh "Select a session that you want unpacked. Pick a directory to unpack it. Hit the unpack button. This process can take up to 20 minutes."
set fileBrowserHelp " "
set commandLineHelp "This script takes 0 or 1 arguments. The argument must be a writeable directory in which the minc files can be unpacked."
set viewProgressHelp "Creates a window to view the output of the unpacker while it is running."
set mincOnlyHelp "Transfers only minc files and does not unpack to b shorts."
set bugHelp "When copying large files, events (like button presses) are processed only between transfers.\n\nE-mail tony@nmr.mgh.harvard.edu"
set searchHelp "This (temporary) search mechanism tries to match character strings in the display info (like name, time and date), not other stuff in the file headers. Doesn't support logical operators, implied or otherwise."

set monthNames [list ??? Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec ]


# include the progress bar source code

source $noArchBin/progressbar.tcl


#--------------------------------------------------------------------------------------#
#-------------------------    check for dependencies     ------------------------------#
#--------------------------------------------------------------------------------------#

if {[info exists env(ARCHIVE_DIR)]}  \
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


if {! [file readable $archiveDir ] } \
    { 
       tk_messageBox -type ok -default ok -title "Command Line Help" \
					      -icon error \
                                              -message "Could not find archive directory $archiveDir"
       exit
    }



if { ! [file executable $noArchBin/unpackmincdir ] } \
     {

               tk_messageBox -type ok -default ok -title "Error" \
                           -message "Could not find unpackmincdir" \
                           -icon error
	        exit
              
  }

catch {exec which mincinfo} exitStatus
if { ! [ file executable $exitStatus ] } \
 {
     
    tk_messageBox -type ok -default ok -title "Error" \
                           -message $exitStatus \
                           -icon error
    exit
  }



#--------------------------------------------------------------------------------------#
#--------------------------------    subroutines     ----------------------------------#
#--------------------------------------------------------------------------------------#



proc Dialog_Create {top title args} {
	global dialog
	if [winfo exists $top] {
		switch -- [wm state $top] {
			normal {
				# Raise a buried window
				raise $top
			}
			withdrawn -
			iconic {
				# Open and restore geometry
				wm deiconify $top
				catch {wm geometry $top $dialog(geo,$top)
                               }
			}
		}
		return 0
	} else {
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

proc SearchDialog {} \
    {
        global  dialog searchString searchDialogClosed
        set sd .searchDialog

        if {[Dialog_Create $sd "Search" -borderwidth 10 ] } \
	   {
                             
              set sd .searchDialog
              set ef [ frame $sd.entryBoxFrame ]

              label $ef.searchEntryBoxlabel -text "Search"
              entry $ef.entryBox -textvariable searchString -relief sunken 

              set bf [ frame $sd.buttonFrame ]
              button $bf.cancelButton -text Cancel -command { set searchString ""
                                                              set searchDialogClosed 1
                                                              Dialog_Dismiss .sd 
                                                            }
              button $bf.submitButton -text Submit -command { 
                                                              set searchDialogClosed 1
                                                              Dialog_Dismiss .sd 
                                                             }
              pack $ef.searchEntryBoxlabel
              pack $ef.entryBox
              pack $bf.cancelButton -side left
	      pack $bf.submitButton 
       
              pack $ef -fill both -expand true
              pack $bf

              bind $ef.entryBox <Return> { 
                                            set searchDialogClosed 1
                                            Dialog_Dismiss .sd 
                                         }
              bind $sd <Control-c> { Dialog_Dismiss .sd }
              
          }

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

              set log [ text $tf.log -width 80 -height 20 \
                       -borderwidth 2 -relief raised \
                       -setgrid true \
                       -yscrollcommand "$tf.yScrollbar set" ]
              scrollbar $tf.yScrollbar -orient vert -command "$tf.log yview"
              pack $tf.yScrollbar -fill y -side right
              pack $tf.log -side top -fill both -expand true

              set bf [ frame $f.buttonFrame ]
              button $bf.closeButton -text Close -command { Dialog_Dismiss .outputWindow }
              button $bf.clearButton -text Clear -command { $log delete 1.0 end }
              pack $bf.clearButton -side left
	      pack $bf.closeButton 
       
              pack $tf -fill both -expand true
              pack $bf
              bind $f <Control-c> { Dialog_Dismiss .outputWindow }
              
          }
        .commandFrame.viewProgressButton config -state disabled
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
				         set response [tk_messageBox -type yesno \
                                                      -default yes -title "Warning" \
						 -message "make $fbInfo(currentDir)" \
					         -icon warning ]
                                         if {[string match yes $response]} \
					     { file mkdir $fbInfo(currentDir)}
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

#--------------------------   GetSessionDescriptors   --------------------------------#

proc GetSessionDescriptors { mincFile } \
    {
	if {[ catch { exec mincinfo $mincFile -attvalue patient:full_name } patientName]} \
	    { set patientName "unknown" }

	if {[ catch { exec mincinfo $mincFile -attvalue study:start_time } startTime]} \
	    { set startTime 000000.000 }

	if {[ catch { exec mincinfo $mincFile -attvalue study:start_date } startDate]} \
	    { set startDate 00000000 }

        set padSize [expr 30 - [ string length $patientName ] ]
        while {$padSize} { append padString " "; incr padSize -1 }
        
        #puts "$startDate $startTime"
        set time [ GetNiceDate $startDate $startTime ]
        set sessionDescriptor "$patientName $padString $time"
       return $sessionDescriptor
    }


#--------------------------   Convert date   --------------------------------#


proc GetNiceDate {date time} \
    {
        global monthNames

	set year   [ string range $date 0 3 ]
        set month  [ string range $date [ expr 5 - [ string index $date 4 ] ] 5 ]
        set day    [ string range $date 6 7 ]
        set hour   [ string range $time 0 1 ]
        set minute [ string range $time 2 3 ]

        return "$hour:$minute    $day [lindex $monthNames $month ]  $year"
    }

#--------------------------   ReadIndexFile   --------------------------------#

proc ReadIndexFile { indexFile } \
    {
        global sessionDescriptor sampleSessionFiles
        set x 0

	if { [file readable $indexFile ] } \
	  {
            set indexFileHandle [ open $indexFile r]
            foreach line [split [read $indexFileHandle ] \n ] \
	      {
		  set fields [ split $line \t ]
	          set sessionDescriptor($x)  [lindex $fields 0]
                  set sampleSessionFiles($x) "[lindex $fields 1]/null"
	          #puts "$x $sessionDescriptor($x) $sampleSessionFiles($x)"
                  incr x
	      }
           }
	return [expr $x-1]

    }


#--------------------   Check src and dest dirs for goodness   --------------------#

proc CheckDirOK {} \
    {
        global targetDir sourceDir targetDirWritable
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
                             return 0
                           }
                        
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


#----------------------------------------- ---------------------------------#


proc LoadListBox {} \
    {
      global sessionInfoFrame progressValue \
             indexFile archiveDir \
             sessionDescriptor sampleSessionFiles mySessionDescriptor mySampleSessionFiles \
             numberOfTotalSessions numberOfIndexedSessions numberOfNewDirs


	#-------- LOAD LISTBOX FROM INDEX ---------#

if { [file exists $indexFile ] } \
   { 
      set numberOfIndexedSessions [ReadIndexFile $indexFile]

      #load index into list box
      for {set index 0} { $index < $numberOfIndexedSessions } {incr index} \
          {
            $sessionInfoFrame.sessionList insert end $sessionDescriptor($index)
          }

    }
update

           #------------ GET DIRECTORY LISTING ---------#

set fileList [ exec ls $archiveDir ] 
#puts $fileList

set index 0
foreach fileName $fileList \
   { 
      
      if { [ file isdirectory [ file join $archiveDir $fileName ] ] && \
           [ file readable    [ file join $archiveDir $fileName ] ] } \
         { 

           set dirList($index) [ file join $archiveDir $fileName ]
           #puts $dirList($index)
           incr index
         }
   }

set numberOfNewDirs $index
set numberOfTotalSessions [expr $numberOfNewDirs + $numberOfIndexedSessions]
# puts "numberOfNewDirs $numberOfNewDirs"
# puts "numberOfIndexedSessions $numberOfIndexedSessions"
# puts "numberOfTotalSessions $numberOfTotalSessions"

              #------- GET HEADER INFO FOR NEW SESSIONS --------#

if { [ set savedCursor [lindex [$sessionInfoFrame configure -cursor] 4]] != "watch"} \
            { $sessionInfoFrame configure -cursor watch }

CreateProgressBarWindow "Reading incoming directory"

set index $numberOfIndexedSessions
for {set dirNumber 0} { $dirNumber < $numberOfNewDirs } {incr dirNumber} \
  {
    set dirContents [ exec ls $dirList($dirNumber) ]
     foreach bareFileName $dirContents   \
      {
	set fileName [ file join $dirList($dirNumber) $bareFileName ]
	  if { [ file isfile $fileName ] && [ file readable $fileName ] } \
           { set sampleSessionFiles($index) $fileName; incr index; break }
      }
     set progressValue [expr $dirNumber*100/$numberOfNewDirs]
     update
  }

set numberOfTotalSessions $index
#puts "numberOfTotalSessions $numberOfTotalSessions"
#.myProgressBarWindow.progressLabel configure -text "loading listbox"

              #-----------   LOAD LISTBOX   ---------------#

for {set index $numberOfIndexedSessions} { $index < $numberOfTotalSessions } {incr index} \
    {
      #create master listing
	
      set sessionDescriptor($index) [ GetSessionDescriptors $sampleSessionFiles($index) ]

      #push entry onto listbox
      $sessionInfoFrame.sessionList insert end $sessionDescriptor($index)
      #set progressValue [expr $index*100/$numberOfTotalSessions]
    }

 update
DestroyProgressBarWindow

           #-----------   SYNCHRONIZE MASTER & DISPLAY ARRAYS   ------------#

for {set index 0} { $index < $numberOfTotalSessions } {incr index} \
    {
      set mySessionDescriptor($index) $sessionDescriptor($index)
      set mySampleSessionFiles($index) $sampleSessionFiles($index)
      
    }

           #-----------   turn on disabled features   ------------#

  $sessionInfoFrame configure -cursor $savedCursor

  .menubar.mOption entryconfigure 1 -state normal
  .menubar.mOption entryconfigure 2 -state normal
  .menubar.mOption entryconfigure 3 -state normal
}


#-----------------------   SEARCH DESCRIPTORS FOR MATCH   --------------------------#

proc SearchDescriptors {} \
    {
       global sessionInfoFrame searchString numberOfTotalSessions \
              sessionDescriptor sampleSessionFiles mySessionDescriptor mySampleSessionFiles 

       set x 0
       $sessionInfoFrame.sessionList delete 0 end

       for {set index 0} { $index < $numberOfTotalSessions } {incr index} \
          {
   	   if { [regexp -nocase $searchString $sessionDescriptor($index)] } \
	     {
               set mySessionDescriptor($x) $sessionDescriptor($index)
               set mySampleSessionFiles($x) $sampleSessionFiles($index)
               $sessionInfoFrame.sessionList insert end $mySessionDescriptor($x)
               incr x
	      }
	  }
    }


#----------------------------------------------------------------------------#

proc TransferMINCfiles {destinationDir} \
  {
    global sourceDir copyDone log
    set dirContents [ exec ls $sourceDir ]
    set numberOfFilesTransferred 0

    foreach fileName $dirContents \
      {
	#puts "cp ${sourceDir}/${fileName} ${destinationDir}/${fileName}"
        if { $copyDone } \
           { 
             $log insert end "Transfer aborted\n$numberOfFilesTransferred files transferred\n"
             $log see end
             update
             return 1 
           }

        if {[catch {exec cp ${sourceDir}/${fileName} ${destinationDir}/${fileName} } \
			 errorMsg ]} {puts stderr $errorMsg}
        incr numberOfFilesTransferred       
	$log insert end "${destinationDir}/${fileName}\n"
        $log see end
        update
      }

    $log insert end "\nnumberOfFilesTransferred files transferred\nTransfer complete\n"
    $log see end
    update
    set sourceDir $destinationDir
  }



#==============================================================================#
#--------------------------    Define widgets    ------------------------------#
#==============================================================================#


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
			      { 
                                set archiveDir $fbInfo(currentDir) 

	                         if { [CheckDirOK] } \
                                    { 
                                      $commandFrame.unpackButton configure -state normal 
                                    } \
                                 else { $commandFrame.unpackButton configure -state disabled }
			      }

                         DestroyFileBrowser
                         $sessionInfoFrame.sessionList delete 0 end
                         LoadListBox 
 
                        .menubar.mFile entryconfigure 1 -state normal
                        .menubar.mFile entryconfigure 2 -state normal
                       }

$File add command -label "New Destination" -command \
                       { 
			set savedTargetDir targetDir
                        set newDirSelected 0

                        .menubar.mFile entryconfigure 1 -state disabled
                        .menubar.mFile entryconfigure 2 -state disabled

			CreateFileBrowser "Destination" $targetDir
			tkwait variable newDirSelected

			   if { $newDirSelected == 1 } \
			      { set targetDir $fbInfo(currentDir)

	                         if { [CheckDirOK] } \
                                    { $commandFrame.unpackButton configure -state normal } \
                                  else { $commandFrame.unpackButton configure -state disabled 
                                          set targetDir $savedTargetDir
                                       }
			      }

                         DestroyFileBrowser
 
                        .menubar.mFile entryconfigure 1 -state normal
                        .menubar.mFile entryconfigure 2 -state normal

                       }

$File add command -label Quit -command {  exit   }


$Option add command -label "Refresh" -state disabled -command { 
                                                $sessionInfoFrame.sessionList delete 0 end
                                                LoadListBox 
                                              }

$Option add command -label "Search" -state disabled -command {
                                      set searchDialogClosed 0
                                      SearchDialog
                                      vwait searchDialogClosed
                                     if { [string compare $searchString "" ] } \
				          {
                                             SearchDescriptors

                                           }
                                        Dialog_Dismiss .searchDialog
                                        }

$Option add command -label "debug info" -command { tk_messageBox -type ok \
                                              -default ok -title "debug" \
					      -message "indexed sessions: $numberOfIndexedSessions\nnew sessions: $numberOfNewDirs\nArchive Dir: $archiveDir\nDest Dir: $targetDir" }

$Help add command -label "Eh?" -command { tk_messageBox -type ok \
                                              -default ok -title "Eh?" \
					      -icon question \
                                              -message $eh }

$Help add command -label "file browser" -command { tk_messageBox -type ok \
                                              -default ok -title "file browser?" \
					      -icon question \
                                              -message $fileBrowserHelp }

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

$Help add command -label "search" -command { tk_messageBox -type ok \
                                                -default ok -title "search" \
					        -icon question \
					        -message $searchHelp
                                             }

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



label .leftFrame.targetDirEntryBoxlabel -text "Destination Directory"
entry .leftFrame.targetDirEntryBox -textvariable targetDir -relief sunken 


#-----------------------------------------------------------------------#

set commandFrame [ frame .commandFrame ]

set commandFrameButtonWidth 14
button $commandFrame.stopButton -text "stop unpacking" \
                                -width $commandFrameButtonWidth \
                                -state disabled \
                                -command \
                                {
				  set copyDone 1
				  if { [ info exists unpackmincdirPipe] } \
                                       { KillPipe $unpackmincdirPipe }
                                  
                                  $commandFrame.stopButton config -state disabled
                                  $commandFrame.unpackButton config -state normal
                                 }

button $commandFrame.viewProgressButton -text "view progress" \
                                -width $commandFrameButtonWidth \
                                -command CreateOutputMonitor

checkbutton $commandFrame.mincOnlyCheckButton -text "minc only" -variable mincOnly

button $commandFrame.unpackButton -text "start unpacking" \
                                -state disabled \
                                -width $commandFrameButtonWidth \
                                -command \
     {  
       if { [CheckDirOK] } \
       { 
        set pipeState 0
        set copyDone 0
        $commandFrame.unpackButton config -state disabled
        $commandFrame.stopButton config -state normal

         switch -exact $transferType \
	 {
	     minc {
                   #puts "unpackmincdir -src $sourceDir -targ $targetDir -minconly"
                   set unpackmincdirPipe [open "|$noArchBin/unpackmincdir -src $sourceDir -targ $targetDir -minconly" ] 
	           fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
	           fconfigure $unpackmincdirPipe -blocking 0
	           vwait pipeState  
                   $commandFrame.stopButton config -state disabled
                   $commandFrame.unpackButton config -state normal 
	          }

	     bshorts {
                   #puts "unpackmincdir -src $sourceDir -targ $targetDir"
                   set unpackmincdirPipe [open "|$noArchBin/unpackmincdir -src $sourceDir -targ $targetDir" ] 
	        fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
	        fconfigure $unpackmincdirPipe -blocking 0
	        vwait pipeState  
                $commandFrame.stopButton config -state disabled
                $commandFrame.unpackButton config -state normal 
	           }

	     rawminc { TransferMINCfiles $targetDir }

        }
       }
     }
        #----------    Transfer type   --------------#

frame       $commandFrame.spacer1 -height 15

set tranferTypeFrame [ frame $commandFrame.tranferTypeFrame -borderwidth 2 -relief ridge  ]

label       $tranferTypeFrame.transferTypeLabel -text "Transfer type"
radiobutton $tranferTypeFrame.bshortsRadioButton -text "b shorts" \
                                                 -value bshorts -variable transferType
radiobutton $tranferTypeFrame.mincRadioButton -text "minc(sorted)" \
                                              -value minc -variable transferType
radiobutton $tranferTypeFrame.rawMincRadioButton -text "minc"   \
                                             -value rawminc  -variable transferType

$tranferTypeFrame.bshortsRadioButton select



#===================================================================================#
#------------------------------        PACK        ---------------------------------#
#===================================================================================#

#  pack transfer radio buttons  #
pack $tranferTypeFrame.transferTypeLabel   -pady 5
pack $tranferTypeFrame.bshortsRadioButton  -anchor w
pack $tranferTypeFrame.mincRadioButton     -anchor w
pack $tranferTypeFrame.rawMincRadioButton  -anchor w


#  pack sessionFrame  #
pack $sessionInfoFrame.sessionList -side left -fill both -expand true 
pack $sessionInfoFrame.leftScroll  -side left -fill y -expand true 

#  pack leftframe  #
pack $sessionInfoFrame
pack .leftFrame.targetDirEntryBoxlabel -anchor center
pack .leftFrame.targetDirEntryBox -side top -fill x -expand true

# pack command column #
pack $commandFrame.unpackButton
pack $commandFrame.stopButton 
pack $commandFrame.viewProgressButton
#pack $commandFrame.mincOnlyCheckButton -anchor w
pack $commandFrame.spacer1
pack $tranferTypeFrame -ipadx 5 -ipady 5


#pack the two topmost frames
pack .leftFrame -side left -fill both -expand true
pack .commandFrame -side right -anchor n


#==================================================================#
#-----------------------    BINDINGS  -----------------------------#
#==================================================================#

bind $sessionPick <ButtonRelease-1> \
  { 
    set sessionIndex [ $sessionInfoFrame.sessionList curselection ]
    set sourceDir [ file dirname $mySampleSessionFiles($sessionIndex) ]
    #puts $sourceDir
      #if { [CheckDirOK] } { $commandFrame.unpackButton configure -state normal }
    $commandFrame.unpackButton configure -state normal
  }

bind .leftFrame.targetDirEntryBox <Return> \
    {
        set targetDirWritable 0
	if { [CheckDirOK] } { $commandFrame.unpackButton configure -state normal } \
         else { $commandFrame.unpackButton configure -state disabled }
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
                       -message "usage:  unpack_mnc.tcl \[distination dir\]" -icon warning
    exit
  }

if { $argc == 1 } \
  {
    set arg [ lindex $argv 0 ]
      if { [ string match {*help*} $arg ] } \
         {
          tk_messageBox -type ok -default ok -title "Command Line Help" \
					      -icon question \
                                              -message $commandLineHelp
	  } \
     else { if {[file writable $arg]} { set targetDir $arg }  }
  }



LoadListBox

CreateOutputMonitor

#wm withdraw .outputWindow 

#$log insert end "New feature: Display fewer sessions by searching (under \"option\" pull-down menu)\n\n"
$log see end


#===========================================================================#
#-------------------------------    TO DO   --------------------------------#
#===========================================================================#
set temp " \
(1) take greatcoat to cleaners \
(2) Boston Edison              \
(3) get minc on demand         \
(4) search sql database         \
"
