#!/usr/bin/wish
# this cript looks at the headers of ima files in a predetermined archive directory. The default archive directory is over-ridden by the environment variable ARCHIVE_DIR. the user selcets a session and the path to that session is proveded to the script that copies the relevant files via nfs to the local machine and unpacks them locally into b shorts.
#--------------------------------------------------------------------------------------#
#Structure of this script:

#global variables
#subroutines
#   dialogs
#    functions
#widgets
#bindings
#main
#--------------------------------------------------------------------------------------#
#export ALPHA_BIN=/space/annecy/5/users/inverse/freesurfer_alpha/bin
#cp unpack_ima.tcl $ALPHA_BIN/noarch
#export PATH=$PATH:$ALPHA_BIN/noarch
#export PATH=$PATH:$ALPHA_BIN/Linux

#----------------------------------------------------------------------------------------#
# global variable defined elsewhere

# sessionNumber: integer 0-N, numerical order of sessions read
#
# session(sessionNumber, fieldValue) A 2-D array containing all that is known about the sessions. 
#   The field Keys are:
#     0 patient_name
#     1 patient_ID
#     2 study_date
#     3 registration_date
#     4 registration_time
#     5 experimenter
#     6 experiment_name
#     7 path to first file in session directory
#
# displayedSessions A list containing the "session" indices which are currently displayed 
#    in the list box


#-----------------------------------------------------------------------------------------#

#define needed variables

set targetDir $env(HOME)
set archiveDir /space/sharbot/1/siemens
#set archiveDir /space/bourget/1/siemens
set sourceDir $archiveDir
set indexFile  [file join [file dirname $archiveDir] index2.txt ]


set env(MRI_DIR) /space/annecy/5/users/inverse/freesurfer_alpha
set noArchBin "$env(MRI_DIR)/bin/noarch"
set archBin "$env(MRI_DIR)/bin/[exec uname]"


set mincOnly 0
set copyDone 1

set typeList { \
               {"Siemens file" ".ima" "SIEM" } \
               {"MINC file" ".mnc" "MINC" } \
               {"bshort file" ".bshort" "BSHT" } \
               {"all" "*" {} } \
              }
array set knownSequences {
                            scout scout
                            ep2d  bold
                            epi   bold
                            ep3d  bold
                            mpr   3danat
                            se    t1conv
	                    tse7  t2conv
                            tse5  protoden
                            flash flash
                            default default
	                  }
set sessionIndex 0;         #default selection in listbox
set targetDirWritable 0;    #default permission
set sequenceDialogClosed 0; #default state: dialog is open?
set lastSessionScanned -1;  # which session is the current sequnece info for?
set percentOfSequencesScanned 0
set patientName ""
        set afterDate   01
        set afterMonth  01
        set afterYear   1990
        set currentTime [clock seconds]
        set beforeDate  [clock format $currentTime -format %d ]
        set beforeMonth [clock format $currentTime -format %m ]
        set beforeYear  [clock format $currentTime -format %Y ]

#--------------------------------  dicom send strings   ----------------------------------#

set dicomSendCommand "${archBin}/siemens_dicom_send"
set mosaicString "-mosaic mni/mosaic64_ 64 64 -mosaic mos64_ 64 64 -mosaic ep2d_fid_ts_20b2604 64 64 -mosaic ep2d_fid_ts_15b3125 64 64 -mosaic ep2d_se_ts_39b2604 64 64 -mosaic ep2d_se_ts_15b3125 64 64 -mosaic ep2d_irm_ts_29b3125 64 64 -mosaic ep2d_irm_ts_39b2604 64 64 -mosaic ep2d_asl_s 64 64"
set portNumber 50082
set archiveHost bourget

#------------------------------    help strings    ---------------------------------------#

set eh "\
Select a session that you want unpacked.
Pick a directory to unpack it. Hit the 
unpack button."

set commandLineHelp "\
This script takes 0 or 1 arguments.
The argument must be a writeable directory
in which the minc files can be unpacked."

set changeDestinationDirHelp "\
Use this dialog to select the directory where
you want session files to be transferred."

set changeArchiveDirHelp "\
Use this dialog to select the directory where
you want to transfer session files from. This
can be any directory containing ima session
directories, including a local one.
Use \"View Recent Pushes\" in order to view that
directory.

If the new directory is on an AMD-mediated
NFS-mounted  remote filesystem, you not be able
navigate by mouse. You will need to
type at least the first three directory levels, 
     e.g. /space/bourget/1

This is because AMD cannot see filesystems that
have not yet been mounted."

set viewRecentPushesHelp "\
View the contents of the archive directory that receives
pushes from the scanners. These sessions might or might
not be in the archive database. No session is deleted
from the incoming until copies have been validated on 
at least two other filesystems. This policy is unrelated
to the scheduling of database updates."

set viewSliceHelp "\
Activates tkmedit. You select a file in a series
that you wish to view."

set viewProgressHelp "\
Creates a window to view the output of the
transfer/translation programs while they are running.
The window can be dismissed, but it continues to 
log information."

set viewHeaderHelp "\
Left-click on a subject to select.

Right-click to see header info.

The output is a reading of pre-selected header
fields/values"

set viewSequencesHelp "\
View the sequences used in the selected session.
Transferring in b-short format causes images to be
separated in the destination directory according to
sequence class.

The count shows the number images generated by a given
sequence.

The destination directory for each sequence is under the
destination directory for the session. This cannot be 
changed. The destination directory for the sequence can 
be changed by pushing on the destination button for the 
sequence."

set loadTextDBhelp "\
There is a text file on the archive that contains 
information about all the sessions on all archive
file systems. This file is used is used to generate 
the database. If the database server is unavailable, 
this file can be read directly."

set stateInfoHelp "\
Primarily for debugging purposes."

set transferIMAhelp "\
This option will copy ima files directly, with no triage."

set transferBshortHelp "\
Temporarily disabled

Translate ima files to bshort format. The integrity of the
output has not yet been independently verified.
Bshort files are sorted according to sequence type."

set transferMINChelp "\
Temporarily disabled

Translate ima files to minc format. There is a known bug
in this translation."

set transferdicomMinchelp "\
Temporarily disabled

Pushes the ima files to a dicom server on the archive host.
The server writes minc files to /space/incoming/minc
To access these files, use the utility \"browse-sessions\""

set bugHelp "\
This is a pre-release version.
E-mail tony@nmr.mgh.harvard.edu
to report errors."

#-----------------------------------------------------------------------------------------#
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


if { ! [file executable $archBin/mri_info ] } \
    {
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

        if {[Dialog_Create $f "Event Log" -borderwidth 10 ] } \
	   {
                             
              set f .outputWindow
              set tf [ frame $f.textBoxFrame ]

              set log [ text $tf.log -width 80 -height 20 \
                                     -borderwidth 2 -relief raised \
                                     -setgrid true \
                                     -yscrollcommand "$tf.yScrollbar set" ]
              scrollbar $tf.yScrollbar -orient vert -command "$log yview"
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
        global  dialog searchString searchDialogClosed \
                patientName \
                afterDate afterMonth afterYear beforeDate beforeMonth beforeYear

        set sd .searchDialog

        #puts "$beforeDate $beforeMonth $beforeYear"

        if {[Dialog_Create $sd "Search" -borderwidth 10 ] } \
	   {
                             
              set sd .searchDialog

              set nameFrame [ frame $sd.nameBoxFrame ]
              label $nameFrame.searchEntryBoxlabel -text "Patient name"
              entry $nameFrame.nameEntryBox -textvariable patientName -relief sunken 

              set dateFrame  [ frame $sd.dateBoxFrame ]
              set titleLabel [ label $dateFrame.titleLabel -text ""]
              set dateLabel  [ label $dateFrame.dateLabel  -text "DD"]
	      set monthLabel [ label $dateFrame.monthLabel -text "MM"]
	      set yearLabel  [ label $dateFrame.yearLabel  -text "YY"]
              grid  $titleLabel $dateLabel $monthLabel $yearLabel

              label $dateFrame.afterLabel -text "After"
              entry $dateFrame.afterDate -textvariable afterDate -relief sunken -width 2
              entry $dateFrame.afterMonth -textvariable afterMonth -relief sunken -width 2
              entry $dateFrame.afterYear -textvariable afterYear -relief sunken -width 4
              grid $dateFrame.afterLabel $dateFrame.afterDate \
                   $dateFrame.afterMonth $dateFrame.afterYear

              label $dateFrame.beforeLabel -text "Before"
              entry $dateFrame.beforeDate -textvariable beforeDate -relief sunken -width 2
              entry $dateFrame.beforeMonth -textvariable beforeMonth -relief sunken -width 2
              entry $dateFrame.beforeYear -textvariable beforeYear -relief sunken -width 4
              grid $dateFrame.beforeLabel $dateFrame.beforeDate \
                   $dateFrame.beforeMonth $dateFrame.beforeYear

              set spacer1 [ frame $sd.spacer1 -height 20 ]

              set buttonFrame [ frame $sd.buttonFrame ]
              button $buttonFrame.cancelButton -text Cancel -command { set searchString ""
                                                              DestroySearchDialog
                                                            }
              button $buttonFrame.submitButton -text Search -command { QueryDB }

    #--- pack ----#
              pack $nameFrame.searchEntryBoxlabel $nameFrame.nameEntryBox   -side left
              pack $buttonFrame.cancelButton      $buttonFrame.submitButton -side left
              
              pack $nameFrame -fill both -expand true
              pack $dateFrame
              pack $spacer1
              pack $buttonFrame

   #--- bind ----#
              bind $nameFrame.nameEntryBox <Return> { QueryDB }
	      bind $sd <Control-c>       { set searchString ""; DestroySearchDialog }
              
          }

    }

proc DestroySearchDialog {} \
  { 
    global searchDialogClosed 
    set searchDialogClosed 1
    destroy .searchDialog 
  }
 

#--------------------------------------------------------------------------------------#

proc CreateSessionInfoView {myIMAfile} \
    {
        global  dialog session sessionInfoViewClosed archBin
        set siv .sessionInfoView

        if {[Dialog_Create $siv "Session Info" -borderwidth 10 ] } \
	   {
                             
              set siv .sessionInfoView
              set infoFrame [ frame $siv.infoFrame ]
              set t [ text $infoFrame.t -setgrid true -wrap word -width 60 -height 30 \
			  -yscrollcommand "$infoFrame.sy set " ]
              scrollbar $infoFrame.sy -orient vert -command "$infoFrame.t yview"

              button $siv.closeButton -text Close -command { DestroySessionInfoView }

# -- pack  ---#
              pack $infoFrame.sy -side right -fill y
              pack $infoFrame.t -side left -fill both -expand true

              pack $infoFrame
              pack $siv.closeButton -anchor center

# -- bind  ---#
              bind $siv <Control-c> DestroySessionInfoView
#              bind $siv <Alt-w>     DestroySessionInfoView

# -- tags  ---#
              $t tag configure bold -font {times 12 bold}

# -- read file  ---#
	      if { [ catch {open "|$archBin/mri_info $myIMAfile" r} IN ] } \
                 {$t insert end "$IN"; return 1}
       
	      set headerLines [ split [ read $IN ] \n ]
	      close $IN

#--  print file header in text box  ---#
# 		    set start [ $t search -count cnt -regexp "\:.*" 1.0 end ]
# 	            $t tag add bold [expr $start + 1] $cnt

              foreach line $headerLines \
		  {  
		      set linePart [ split $line : ]
		      $t insert end "[ lindex $linePart 0 ] : " 
		      $t insert end "[ lindex $linePart 1 ]\n" bold
		}
           #end of dialog create
	  }
     #end of function
    }

proc DestroySessionInfoView {} \
  { 
    global sessionInfoViewClosed 
    set sessionInfoViewClosed 1
    destroy .sessionInfoView
  }


#--------------------------------------------------------------------------------------#

proc CreateAlertDialog {title alertMessage} \
    {
        global  dialog session alertDialogClosed archBin
        set ad .alertDialog
        set longestLine 0
        set numberOfLines 30

# -- determine size  ---#
	set textLines [ split $alertMessage \n ]
        foreach line $textLines \
	    {
		set lineLength [string length $line ]
		if { $lineLength > 60} { set lineLength 60 }
                if { $lineLength > $longestLine } { set longestLine $lineLength }
	    }
	if { $numberOfLines > [llength $textLines ] } { set numberOfLines [llength $textLines ] }

        if {[Dialog_Create $ad $title -borderwidth 10 ] } \
	   {
                             
              set ad .alertDialog
              set msgFrame [ frame $ad.msgFrame ]
              set t [ text $msgFrame.t -setgrid true -wrap word \
			               -width [expr $longestLine + 2] \
                                       -height [expr $numberOfLines + 1] \
			         -yscrollcommand "$msgFrame.sy set" ]
              scrollbar $msgFrame.sy -orient vert -command "$msgFrame.t yview"

              button $ad.closeButton -text Close -command { DestroyAlertDialog }

# -- pack  ---#
              pack $msgFrame.sy -side right -fill y
              pack $msgFrame.t -side left -fill both -expand true

              pack $msgFrame
              pack $ad.closeButton

# -- bind  ---#
              bind $ad <Control-c> DestroyAlertDialog
              bind $ad <Alt-w>     DestroyAlertDialog
              bind $ad <Return>     DestroyAlertDialog
# -- tags  ---#
              $t tag configure bold -font {times 12 bold}


#--  print msg in text box  ---#


              foreach line $textLines \
		  {  

                    $t insert end "$line\n"
		}
           #end of dialog create
	  }
     #end of function
    }

proc DestroyAlertDialog {} \
  { 
    global alertDialogClosed
    .alertDialog.msgFrame.t delete 1.0 end
    set alertDialogClosed 1
    destroy .alertDialog
  }


#-----------------------------   CreateSequenceDialog   ---------------------------------#


proc CreateSequenceDialog { sessionNumber} \
    {
        global  dialog archBin targetDir newDirSelected fbInfo \
                session sequenceFiles sequenceTypes knownSequences \
                percentOfSequencesScanned sequenceProgressLabel \
                sequenceDialogClosed lastSessionScanned copyDone
                
       #sequenceFiles:   array of lists. each key is a seqtype. Each list contains filenames
       #sequenceTypes:  list of sesionsequence types like ep2d, mpr, se, etc
       #knownSequences: associative array of sequence types to dirNames
	#fbInfo:         associative array; keys are currentDir and currentFile
       #puts "lastSessionScanned $lastSessionScanned"
       #puts "sessionNumber $sessionNumber"
       set sequenceProgressLabel "Percent of\nsession scanned"
       set rowHeight 40
       set columnWidth 100
       set seqd .sequenceDialog


        if {[Dialog_Create $seqd Sequences -borderwidth 10 ] } \
	   {
              set seqd .sequenceDialog
              set sequenceDialogClosed 0
	      set progressFrame [frame $seqd.progressFrame ]
              label $progressFrame.progressLabel -textvariable sequenceProgressLabel
              label $progressFrame.progressValue -textvariable percentOfSequencesScanned
              button $progressFrame.cancelButton -text Cancel \
                         -command { 
                                    set copyDone 1
                                    DestroySequenceDialog 
                                   }
              pack $progressFrame.progressLabel -side top
              pack $progressFrame.progressValue
              pack $progressFrame.cancelButton
              pack $progressFrame -side right
              update

	      if {$lastSessionScanned != $sessionNumber} \
                 { 
		     if { [GetSequences $sessionNumber] } \
			 { puts "GetSequences failed"; DestroySequenceDialog; return 1}
	          }


              pack forget $progressFrame
              #set sequenceProgressLabel ""
              #set percentOfSequencesScanned ""
              

              set tableFrame [ frame $seqd.tableFrame -relief ridge -borderwidth 4 ]
               
# -- table header  ---#
	      set sequenceLabelFrame [ frame $tableFrame.sequenceLabelFrame \
                                       -relief ridge -borderwidth 2 -bg grey50 \
					   -height $rowHeight -width $columnWidth ]
              grid $sequenceLabelFrame -sticky news -row 0 -column 0
              label $sequenceLabelFrame.label -text Sequence
              pack  $sequenceLabelFrame.label -fill x -fill y

	      set countLabelFrame [ frame $tableFrame.countLabelFrame \
                                    -relief ridge -borderwidth 2 -bg grey50 \
					-height $rowHeight -width $columnWidth ]
              grid $countLabelFrame -sticky news -row 0 -column 1
              label $countLabelFrame.label -text Count
              pack  $countLabelFrame.label -fill x -fill y

	      set dirLabelFrame [ frame $tableFrame.dirLabelFrame \
                                  -relief ridge -borderwidth 2 -bg grey50 \
				      -height $rowHeight -width $columnWidth ]
              grid $dirLabelFrame -sticky news -row 0 -column 2
              label $dirLabelFrame.label -text Dir
              pack  $dirLabelFrame.label -fill both -expand true

# 	      set transferLabelFrame [ frame $tableFrame.transferLabelFrame \
#                                   -relief ridge -borderwidth 2 -bg grey50 \
# 				      -height $rowHeight -width $columnWidth ]
#               grid $transferLabelFrame -sticky news -row 0 -column 3
#               label $transferLabelFrame.label -text Transfer
#               pack  $transferLabelFrame.label -fill both -expand true
              
              #grid  $sequenceLabelFrame $countLabelFrame $dirLabelFrame 

set row 1
# -- table content  ---#
              foreach sequenceType $sequenceTypes \
		{
		  set l [frame $tableFrame.${sequenceType}LabelFrame \
                                   -relief ridge -borderwidth 2 \
                                    -height $rowHeight -width $columnWidth ]
                  grid $l -sticky news -row $row -column 0
                  label  $l.label -text $sequenceType
                  grid $l.label -sticky news

		  set c [frame $tableFrame.${sequenceType}CountFrame \
                                       -relief ridge -borderwidth 2 \
                                       -height $rowHeight -width $columnWidth ]
                  grid $c -sticky news -row $row -column 1
		  label  $c.label \
                      -text [llength $sequenceFiles($sequenceType) ]
                  grid $c.label -sticky news

		  set d [frame $tableFrame.${sequenceType}dirFrame \
                                             -relief ridge -borderwidth 2  \
                                             -height $rowHeight -width $columnWidth ]
                  grid $d -sticky news -row $row -column 2
                  button $d.button -text "$knownSequences($sequenceType)" \
                                   -width 8 \
		                   -command "GetSeqDir $sequenceType"
                  grid $d.button -sticky news

#                   set ch [frame $tableFrame.${sequenceType}checkbuttonFrame \
#                                              -relief ridge -borderwidth 0  \
#                                              -height $rowHeight -width $columnWidth ]
#                   grid $ch -sticky news -row $row -column 3
#                   set ch.checkbuttonState 1
#                   checkbutton $ch.checkbutton -variable ch.checkbuttonState \
# 		      -command "if \${ch.checkbuttonState} \
#                        \"lappend sequenceTypes $sequenceType\" \
#                        else \"DeleteSequenceListItem $sequenceType\"
#                        "
#                   set ch.checkbuttonState 1
#                   grid $ch.checkbutton -sticky news


                  grid  $l $c $d
                  incr row
		}

# --- controls  ---#
              set buttonFrame [ frame $seqd.buttonFrame ]
              button $buttonFrame.okButton -text Okay -command { DestroySequenceDialog }
              button $buttonFrame.cancelButton -text Cancel \
                         -command { 
                                    set copyDone 1
                                    DestroySequenceDialog 
                                   }
# --- misc ---#
              set spacer1 [ frame $seqd.spacer1Frame -height 20]

# -- pack  ---#

              pack $buttonFrame.okButton -side left
              #pack $buttonFrame.cancelButton

              pack $tableFrame
              pack $spacer1
              pack $buttonFrame

# -- bind  ---#
              bind $seqd <Control-c> DestroySequenceDialog
              bind $seqd <Alt-w>     DestroySequenceDialog


           
	  }; #end of dialog create

     return 0
    }; #end of function

proc DestroySequenceDialog {} \
  { 
    global sequenceDialogClosed
    set sequenceDialogClosed 1
    destroy .sequenceDialog
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


            if {  ! [ file exists "$fbInfo(currentDir)" ] && [string match "Destination" $dirType] } \
		{
		    if {[catch "file mkdir $fbInfo(currentDir)" errorMsg ] } \
			{
                          CreateAlertDialog $errorMsg	   
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
# 		if { ! [string compare $line "Done sending files." ] } \
#                    { 
#                        $log insert end "log: $line\n"
#                        $log see end
#                       if {[catch {close $pipeName} errorMsg ]} \
#                          {
# 			     CreateAlertDialog "Error: Log" $errorMsg
# 	                 }
#                       set  pipeState 1
#                       return 0
#                    }

#                if { $line == 0 } { KillPipe $pipeName }

              $log insert end $line\n
              $log see end

              
            }
    }


#--------------------------------------------------------------------------------------#

proc KillPipe {pipeName} \
    {
        global pipeState log


      if {[catch {close $pipeName} errorMsg ]} \
	    {
              $log insert end $errorMsg\n
              $log see end

#              foreach pid [pid $pipeName] \
# 	      {

#                 if {[catch {exec kill $pid} errorMsg ]} \
# 	           { 
#                     CreateAlertDialog "Error" "Could not close process $pid\n$errorMsg"
#                     return 1
# 		   }
#               }
	    }
	set pipeState 1
        return 0 
    }


#----------------------------------   ReadIndexFile   -----------------------------------#
# reads a text file that has session tab-delimited description info and path

proc ReadIndexFile { indexFile } \
    {
        global session numberOfSessions
        set x 0

	if { [file readable $indexFile ] } \
	  {
            set indexFileHandle [ open $indexFile r]
            foreach line [split [read $indexFileHandle ] \n ] \
	      {
                  # if line has less than 8 fields and 7 delimiters, skip
                  if {[string length $line] < 15 } { continue }

                  # break up fields into a list
		  set fields [ split $line \t ]
                  if { [ llength $fields ] != 8 } \
                     { puts "$indexFile:line $x: [llength $fields] fields read, 8 expected"
                       continue
                     }

                  #load 2-D array with each field
	          for {set i 0} {$i < 8} {incr i} \
	             {
	      	        set session($x,$i) [lindex $fields $i]
                        #puts "\$session($x,$i)=$session($x,$i)"
	             }

	          #puts "$x $session($x,7)"
                  incr x
	      }
           }

        set numberOfSessions  [expr $x-1]
        LoadListBox
	return 0

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


#--------------------------   GetSessionInfoFromFile   --------------------------------
# fills the 2-D array "session" with values that this script uses elsewhere.
# It assumes that the 8th field has been filled with the name of a file whose
# header can be used.

proc GetSessionInfoFromFile { sessionIndex } \
    {
	global headerInfo session

# set defaults
        set headerInfo(patient_name)     "unknown"
        set headerInfo(patient_id)       "unknown"
	set headerInfo(study_date)       "unknown"
	set headerInfo(registration_date) "00000000"
        set headerInfo(registration_time) "000000"
        set headerInfo(experimenter)      "unknown"
        set headerInfo(experiment_name)   "unknown"

 # get header key:values
        if { [ catch { GetHeaderInfo $session($sessionIndex,7) } exitStatus ] } \
            { CreateAlertDialog Error $exitStatus; return 1}

	set session($sessionIndex,0) $headerInfo(patient_name)
	set session($sessionIndex,1) $headerInfo(patient_id)
	set session($sessionIndex,2) $headerInfo(study_date)
	set session($sessionIndex,3) $headerInfo(registration_date)
	set session($sessionIndex,4) $headerInfo(registration_time)
	set session($sessionIndex,5) $headerInfo(experimenter)
	set session($sessionIndex,6) $headerInfo(experiment_name)

        update
	
         return 0

    }




#--------------------------   GetSessionDescriptors   --------------------------------#

# returns a formatted string of info suitable for displaying in list box

proc GetSessionDescriptors { i } \
    {
       global session numberOfSessions
        
        set IDpadString ""
        set namePadString ""

	if { $i >= $numberOfSessions } \
           { set i [ expr $numberOfSessions - 1 ] }
        

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



#------------------------------    ReadArchiveIndex    ---------------------------------#

proc ReadIncomingDir {} \
    {

      global sessionInfoFrame progressValue \
             indexFile archiveDir displayedSessions \
             sessionDescriptor sampleSessionFiles mySessionDescriptor mySampleSessionFiles \
             session numberOfSessions 

set fileList [ exec ls $archiveDir ] 
#puts $fileList

# set the full path of each session dir into the list dirList
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

set numberOfNewDirs $index
# puts "numberOfSessions $numberOfSessions"

              #------- GET HEADER INFO FOR NEW SESSIONS --------#

set index 0

if { [ set savedCursor [lindex [$sessionInfoFrame configure -cursor] 4]] != "watch"} \
            { $sessionInfoFrame configure -cursor watch }

CreateProgressBarWindow "Reading incoming directory"


for {set dirNumber 0} { $dirNumber < $numberOfNewDirs } {incr dirNumber} \
  {
      

    set dirContents [ exec ls -1 $dirList($dirNumber) ]
     foreach bareFileName $dirContents   \
      {
	set fileName [ file join $dirList($dirNumber) $bareFileName ]

        #get the first scan as a sample
	if { [ file isfile $fileName ] } \
           { 
             #puts $fileName
             if {[catch {exec mri_info $fileName | grep sequence } exitStatus]} \
                   { 
                     #could not find "sequnce" in the header
                     continue
                   } \
               else \
		   {
                      #found a sequencename, now search for "scout"
                      #puts "$exitStatus"
                      if {[string match "*scout*" $exitStatus]} \
		           { 
                             # ignore this scout
                             #puts "[file tail $fileName]: scout" 
                           } \
                      else \
                         { 
                           #found a session scan
                           #puts "$index: [file tail $fileName] is scan"
                           set session($index,7) $fileName
                           GetSessionInfoFromFile $index
                           incr index
                           break 
                         }
		  }; #end of found sequence name
            }; # end of is file
      
      }; # end of foreach 

    set numberOfSessions $index
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
  LoadListBox
  return 0
}


#-----------   LOAD LISTBOX   ---------------#

proc LoadListBox {} \
    {
       global session sessionInfoFrame numberOfSessions

        $sessionInfoFrame.sessionList delete 0 end
        for {set index 0} { $index < $numberOfSessions } {incr index} \
           {
              $sessionInfoFrame.sessionList insert end [ GetSessionDescriptors $index ]
            }

        update
        return 0
 
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
     set copyDone 1
    }


#----------------------------------------------------------------------------#

proc TransferMINCfiles {destinationDir} \
    {

       global sourceDir copyDone log archBin unpackmincdirPipe

       #set dirContents [ exec ls $sourceDir ]
       #set dirLength [ llength $dirContents ]
       #puts "dirLength $dirLength"

	#puts "cp ${sourceDir}/${fileName} ${destinationDir}/${fileName}"
        if { $copyDone } \
           { 
             $log insert end "Transfer aborted\n"
             $log see end
             update
             return 1 
           }

       if { ! [ file exists $archBin/ima2mnc ] } \
	   {
	        CreateAlertDialog "Error" \
                   "Couldn't find $archBin/ima2mnc"
             return 1
	   }
#       puts "$archBin/ima2mnc -range 1 3 -host bourget -aetitle bay2fmri -port 50082 $sourceDir $destinationDir"
       #puts "$archBin/ima2mnc  $sourceDir $destinationDir"

	set unpackmincdirPipe [open "|$archBin/ima2mnc $sourceDir $destinationDir" ]

        fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
        fconfigure $unpackmincdirPipe -blocking 0

        vwait pipeState

       $log insert end "\n\n"
       $log see end
       update
       
       set copyDone 1
    }


#-----------------------------   DicomPushMinc   ----------------------------------#

proc DicomPushMinc {destinationDir} \
    {
       global session sessionIndex sourceDir copyDone log pipeState unpackmincdirPipe \
              noArchBin archiveHost env \
              dicomSendCommand mosaicString portNumber

        set bareFileName [ file tail $session($sessionIndex,7) ]

	if {[regexp {^[0-9]+} $bareFileName sessionNumber]} \
           { } \
        else \
           {
             CreateAlertDialog "Error" \
                   "Couldn't get session number from $bareFileName"
             return 1
           }
        

        if {[string match "martyrium*" $env(HOSTNAME) ]} \
           {puts "$dicomSendCommand $mosaicString $archiveHost $portNumber bay2fmri $sessionNumber -dir $sourceDir -max_outstanding 0"}

	set unpackmincdirPipe [open "|$dicomSendCommand $mosaicString $archiveHost $portNumber bay2fmri $sessionNumber -dir $sourceDir -max_outstanding 0" ]

      fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
      fconfigure $unpackmincdirPipe -blocking 0
      vwait pipeState

      set copyDone 1
#/home/dicom/bin/siemens_dicom_send bourget 50082 bay2fmri 480 -dir /local_mount/space/bourget/3/siemens/sonata-21006-20001010-123121-480  -max_outstanding 0
    }

#-----------------------------   TransferBshortFiles   ----------------------------------#

proc TransferBshortFiles {destinationDir} \
    {
       #sequenceFiles:   array of lists. each key is a seqtype. Each list contains filenames
       #sequenceTypes:  list of sesionsequence types like ep2d, mpr, se, etc
       #knownSequences: associative array of sequence types to dirNames

    global sequenceFiles sequenceTypes knownSequences \
           commandFrame sourceDir copyDone log archBin progressValue \
           lastSessionScanned sessionIndex
    
    if { $lastSessionScanned != $sessionIndex } \
	{ 
	  if { [CreateSequenceDialog $sessionIndex]} \
	      { 
                 if { $copyDone } {return 1} \
                 else { CreateAlertDialog "Error" "Session Scan failed"; return 1 }
	      }
	}

    set numberOfFilesToTransfer 0
    foreach sequenceName $sequenceTypes \
        { incr numberOfFilesToTransfer [llength $sequenceFiles($sequenceName)] }

    set dirContents [ exec ls $sourceDir ]
    set numberOfFilesTransferred 0
    set dirLength [ llength $dirContents ]
    set progressValue 0

    CreateProgressBarWindow "Transferring b-short files"
    

    foreach sequenceName $sequenceTypes \
      {
          # make dest dir if non-existant
	  if { ! [ file exists ${destinationDir}/$knownSequences($sequenceName)]} \
	      {
		  if {[ catch "file mkdir [file join $destinationDir $knownSequences($sequenceName)]" errormsg] } \
                      { 
                         CreateAlertDialog Error $errormsg
                         continue
		      }

	      }

        foreach fileName $sequenceFiles($sequenceName) \
          {
              set imaFilePath  [file join $sourceDir $fileName ]
              set bshortFilePath [file join $destinationDir [file join $knownSequences($sequenceName) [ file rootname $fileName ].bshort ] ]

	     puts "cp $imaFilePath $bshortFilePath"
             if { $copyDone } \
                { 
                $log insert end "Transfer aborted\n$numberOfFilesTransferred files transferred\n"
                $log see end
                DestroyProgressBarWindow
                update
                return 1 
                }

              if {[catch {exec $archBin/mri_convert $imaFilePath $bshortFilePath }  errorMsg ]} \
                {
                  # error!
                   #puts stderr $errorMsg
                   $log insert end "failed: $bshortFilePath\n"
                   $log insert end $errorMsg
                   $log see end
                 } \
              else \
	        {
                  # successful copy
                   incr numberOfFilesTransferred       
	          $log insert end "$bshortFilePath\n"
                  $log see end
	        }
              set progressValue [ expr $numberOfFilesTransferred*100/$numberOfFilesToTransfer ]
             update
	  }; #end fileName $sequenceFiles($sequenceName)
      }; #foreach sequenceName $sequenceTypes

    $log insert end "\n$numberOfFilesTransferred files transferred\nTransfer complete\n"
    $log see end
    update
    #set sourceDir $destinationDir
     DestroyProgressBarWindow
     set copyDone 1
    
      }; #end of proc 

#-----------------------------   TransferBshortFiles   ----------------------------------#

proc Ima2sessions {destinationDir} \
    {
       global sourceDir copyDone log noArchBin archBin ima2sessionsPipe

	#puts "cp ${sourceDir}/${fileName} ${destinationDir}/${fileName}"
        if { $copyDone } \
           { 
             $log insert end "Transfer aborted\n"
             $log see end
             update
             return 1 
           }

       if { ! [ file exists $noArchBin/unpackimadir ] } \
	   {
	        CreateAlertDialog "Error" \
                   "Couldn't find $noArchBin/unpackimadir"
             return 1
	   }

       #puts "$noArchBin/unpackimadir  -src $sourceDir -targ $destinationDir"

	set ima2sessionsPipe [open "|$noArchBin/unpackimadir -src $sourceDir -targ $destinationDir" ]

        fileevent $ima2sessionsPipe readable { Log $ima2sessionsPipe }
        fconfigure $ima2sessionsPipe -blocking 0

        vwait pipeState

       $log insert end "\n\n"
       $log see end
       update
       
       set copyDone 1
    
      }; #end of proc 

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

#-----------------------   GetSequences   --------------------------#
# fills sequenceFiles (an array of lists containing filenames)

proc GetSequences {sessionNumber} \
    {
       #knownSequences: associative array pairing sequences to dirnames
       #sequenceFiles:   array of lists. each key is a seq type. Each list contains filenames

       global session archBin archiveDir copyDone \
              lastSessionScanned percentOfSequencesScanned \
              knownSequences sequenceFiles sequenceTypes env

       #foreach {key value} [array get knownSequences] {puts "$key $value"}
       if {[info exists sequenceFiles]} {unset sequenceFiles }
       set copyDone 0
       set sequenceTypes [list default]

       set sessionDir [ file dirname $session($sessionNumber,7) ]
       set fileNames  [ split [exec ls -1 $sessionDir] "\n" ]
       #puts "sessionDir $sessionDir"
       #puts "number of images [llength $fileNames]"
        set numberOfFiles [llength $fileNames]
        set currentFileNumber 1

       foreach fileName $fileNames \
	 {
             if { $copyDone } {return 2}

             #update sequence dialog, if open
             set percentOfSequencesScanned [ format "%2.0f" [ expr {$currentFileNumber*1.00 / $numberOfFiles*1.00 } * 100 ] ]

             #puts "$percentOfSequencesScanned"

             set fullFileName [file join $sessionDir $fileName]
             #puts "$fullFileName"
	     if {[catch {exec $archBin/mri_info $fullFileName}   exitStatus]} \
               { 
	         set response [tk_messageBox -type yesno \
                                  -default no -title "Error" \
	                          -icon error \
				  -message "$exitStatus\n\nContinue?" ]
                  if {[string match yes $response]} {continue} else { return 1}

               }

	    set sequenceLine "sequence name: unknown"
            set lines [ split $exitStatus "\n" ]

            foreach line $lines \
              {
                  if {[string match {sequence name*} $line]} \
                     {
                       set sequenceLine $line
                       break
                     } 
	      }

            #puts "$fileName: line is $line"
            set sequenceTypeFound 0
            set temp [split  $sequenceLine : ]
	    set sequencePath [ lindex $temp 1 ]
	    set sequenceProgramName [file tail $sequencePath]
            #puts "$fileName: sequence is $sequenceString"

	    foreach sequenceType [array names knownSequences ] \
	       {
		   if {[ string match "*${sequenceType}*" $sequenceProgramName ]} \
		       {
                           #add this type to current list of types
                           if {[lsearch $sequenceTypes $sequenceType] == -1} \
			      { lappend sequenceTypes $sequenceType }

                         #collect the filename for this type
                         #puts "$sequenceType found"
                         lappend sequenceFiles($sequenceType) $fileName
                         set sequenceTypeFound 1
                         break
                       }
	       }

	     if {! $sequenceTypeFound} {lappend sequenceFiles(default) $fileName}


         incr currentFileNumber
         update
	 }; #end foreach fileNames
       
# remove default from lists if unused
       if { ! [info exists sequenceFiles(default)]} \
	   {
             puts "default dir is empty"
	     set sequenceTypes [ lreplace $sequenceTypes [lsearch $sequenceTypes default] \
				      [lsearch $sequenceTypes default] ]
	   }

# print results - be sure to comment out
   if {[string match "martyrium.*" $env(HOSTNAME)]} \
	   {
	       puts "My dirs"
	       foreach dirName [array names sequenceFiles] \
		   {
		       puts "\t$dirName: [llength $sequenceFiles($dirName)]"
		       if {[catch "open /home/tony/Dicom/tempdir/$dirName w" OUT ]} \
			   {puts $OUT; continue}
		       foreach fileName $sequenceFiles($dirName) {puts $OUT $fileName}
		       close $OUT
		   }

	       puts "My sequences: $sequenceTypes"
	       foreach sequenceType $sequenceTypes \
		   {
		       puts "\t$sequenceType: $knownSequences($sequenceType)"
		   }
	   }

       set lastSessionScanned $sessionNumber
       return 0
    }



#--------------------------    startTransfer    -------------------------------#

proc StartTransfer {} \
     {   
        global commandFrame transferType pipeState copyDone targetDir

        if { ! [CheckDirOK] } \
           { 
             $commandFrame.startButton configure -state disabled 
             return 1
           }
        set pipeState 0
        set copyDone 0
        $commandFrame.startButton config -state disabled
        $commandFrame.stopButton config -state normal

        switch -exact $transferType \
	 {
	     minc    { TransferMINCfiles $targetDir }

             #bshorts { TransferBshortFiles $targetDir }
             bshorts { Ima2sessions $targetDir }

             dicom_minc     { DicomPushMinc $targetDir }

             ima     { TransferIMAfiles $targetDir }

             default { TransferIMAfiles $targetDir }
	 }
	   
         tkwait variable copyDone

         $commandFrame.stopButton config -state disabled
         $commandFrame.startButton config -state normal 
         #puts "StartTransfer exited"
      }


#--------------------------    GetSeqDir    -------------------------------#
#updates the destination diectory associated with a given sequence

proc GetSeqDir { sequenceType } \
    {
        global knownSequences fbInfo targetDir newDirSelected

	set newDirSelected 0
	#disable change archive/destdir options
	.menubar.mFile entryconfigure 1 -state disabled
	.menubar.mFile entryconfigure 2 -state disabled

	CreateFileBrowser "Destination" "${targetDir}/$knownSequences($sequenceType)"
	tkwait variable newDirSelected

	DestroyFileBrowser

        if { $newDirSelected == 2 } { return $knownSequences($sequenceType) }

        set newDir [file tail $fbInfo(currentDir)]
	set knownSequences($sequenceType) $newDir
	#puts "$newDir selected"

        .sequenceDialog.tableFrame.${sequenceType}dirFrame.button configure -text  $newDir
        set $knownSequences($sequenceType) $newDir

	.menubar.mFile entryconfigure 1 -state normal
	.menubar.mFile entryconfigure 2 -state normal
	#"set newDir \[ CreateFileBrowser destination $targetDir \]"
    }


#--------------------------    DeleteSequenceListItem    ------------------------#

proc DeleteSequenceListItem { listItem } \
    {
        global sequenceTypes 
	return [ lreplace $sequenceTypes [ lsearch sequenceTypes $listItem ] ]
    }


#--------------------------    QueryDB    ------------------------#
proc QueryDB {  } \
    {
        global session sessionInfoFrame numberOfSessions \
               auto_path env \
               patientName \
               afterDate afterMonth afterYear beforeDate beforeMonth beforeYear

        .menubar.mFile entryconfigure 1 -state disabled; #change archive
        .menubar.mFile entryconfigure 2 -state disabled; #change destination

        for {set index 0} {$index <= [ .menubar.mView index end] } {incr index} \
	    {
               .menubar.mView entryconfigure $index -state disabled
	    }


        set afterTimeStamp  "${afterYear}${afterMonth}${afterDate}"
        set beforeTimeStamp "${beforeYear}${beforeMonth}${beforeDate}"

        #load the sql package 
        #the file pkgIndex.tcl must be in same dir as sql.so
        lappend auto_path $env(MRI_DIR)/local/bin/[exec uname]
        if {[catch {package require Sql} errorMsg]} \
	    {
              CreateAlertDialog Error $errorMsg
              return 1
            }

       # Connect to the database
	if {[catch {set conn [sql connect bourget query]} errorMsg]} \
	    {
              CreateAlertDialog Error $errorMsg
              return 1
            }


       # select the database to be used.
       sql selectdb $conn archive

       #execute search command
       set searchCommand "SELECT * FROM header_info WHERE patient_name LIKE \'\%$patientName\%\' AND registration_date BETWEEN $afterTimeStamp AND $beforeTimeStamp"

       #puts "$searchCommand"
       sql query $conn "$searchCommand"

       #update internal database
       set index 0
       while {[set row [sql fetchrow $conn]] != ""} \
         {
            set x 0
	     foreach item $row { set session($index,$x) $item; incr x }
            #puts "$session($index,0)"
            incr index

         }

       set numberOfSessions $index
       sql endquery $conn
       sql disconnect $conn

       update
       LoadListBox

        .menubar.mFile entryconfigure 1 -state normal; #change archive
        .menubar.mFile entryconfigure 2 -state normal; #change destination

        for {set index 0} {$index <= [ .menubar.mView index end] } {incr index} \
	    {
               .menubar.mView entryconfigure $index -state normal
	    }


       return 0
    }



#=============================================================================#
#---------------------------- Define widgets   -------------------------------#
#=============================================================================#


menu .menubar
#attach it to the main window
. config -menu .menubar
foreach menuWidget { File View Option Help } \
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


$File add command -label Quit -command {  exit   }


         #--------- View ------------#

$View add command -label "Recent Pushes"   -state normal \
                                           -command { ReadIncomingDir  }

$View add command -label "All sessions"   -state normal \
                                       -command { 
                                    set afterDate   01
                                    set afterMonth  01
                                    set afterYear   1990
                                    set currentTime [clock seconds]
                                    set beforeDate  [clock format $currentTime -format %d ]
                                    set beforeMonth [clock format $currentTime -format %m ]
                                    set beforeYear  [clock format $currentTime -format %Y ]
                                    set patientName ""
                                    QueryDB
				       }


$View add command -label "Slice"  -command \
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

$View add command -label "Sequences" -command {
                                                   #puts "sessionIndex $sessionIndex"
                                                    CreateSequenceDialog $sessionIndex
                                                 }

         #--------- Option ------------#

$Option add command -label "Load text DB"   -state disabled \
                                       -command { 
                                                  set copyDone 1
                                                  ReadIndexFile
                                                }

$Option add command -label "Search"   -state normal \
                                       -command { CreateSearchDialog
				       }

$Option add command -label "state info" -state disabled \
                                        -command { 
                                              CreateAlertDialog "State" \
					       "sessions: $numberOfSessions\nArchive Dir: $archiveDir\nDest Dir: $targetDir" }


         #--------- Help ------------#

$Help add command -label "Eh?" -command { tk_messageBox -type ok \
                                              -default ok -title "Eh?" \
					      -icon question \
                                              -message $eh }


$Help add command -label "command line" \
                  -command { CreateAlertDialog "Command line arguments" $commandLineHelp }


$Help add cascade -label file -menu $Help.file
set FileHelp [menu $Help.file -tearoff 0 ]

$FileHelp add command -label "Change destination dir" \
                  -command { CreateAlertDialog "Change destination dir" \
                             $changeDestinationDirHelp }

$FileHelp add command -label "Change archive dir" \
                  -command { CreateAlertDialog "Change archive dir" \
                             $changeArchiveDirHelp }



$Help add cascade -label View -menu $Help.view

set ViewHelp [menu $Help.view -tearoff 0 ]

$ViewHelp add command -label "view recent pushes" \
                  -command { CreateAlertDialog "view recent pushes" $viewRecentPushesHelp }

$ViewHelp add command -label "view slice" \
                  -command { CreateAlertDialog "view slice" $viewSliceHelp }

$ViewHelp add command -label "view log" \
                  -command { CreateAlertDialog "view log" $viewProgressHelp }

$ViewHelp add command -label "Viewing headers" \
                  -command { CreateAlertDialog "Viewing Headers" $viewHeaderHelp }

$ViewHelp add command -label "Viewing Sequences" \
                  -command { CreateAlertDialog "Viewing Sequences" $viewSequencesHelp }



$Help add cascade -label Option -menu $Help.option

set OptionHelp [menu $Help.option -tearoff 0 ]

$OptionHelp add command -label "Load text DB" \
                  -command { CreateAlertDialog "Load text DB" $loadTextDBhelp }

$OptionHelp add command -label "state info" \
                  -command { CreateAlertDialog "state info" $stateInfoHelp }



$Help add cascade -label Transfer -menu $Help.transfer

set TransferHelp [menu $Help.transfer -tearoff 0 ]

$TransferHelp add command -label "IMA" \
                  -command { CreateAlertDialog "IMA" $transferIMAhelp }

$TransferHelp add command -label "Sessions" \
                  -command { CreateAlertDialog "b short" $transferBshortHelp }

$TransferHelp add command -label "minc" \
                  -command { CreateAlertDialog "minc" $transferMINChelp }

$TransferHelp add command -label "dicom minc" \
                  -command { CreateAlertDialog "dicom minc" $transferdicomMinchelp }



$Help add command -label "bugs" \
                  -command { CreateAlertDialog "bug" $bugHelp }



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

button $commandFrame.startButton -text "start" \
                                -state disabled \
                                -width $commandFrameButtonWidth \
                                -command StartTransfer


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
radiobutton $tranferTypeFrame.bshortsRadioButton -text "Sessions" -value bshorts -variable transferType -state normal
radiobutton $tranferTypeFrame.mincRadioButton -text "MINC" -value minc -variable transferType -state normal

radiobutton $tranferTypeFrame.imaRadioButton -text "IMA"   -value ima  -variable transferType

$tranferTypeFrame.bshortsRadioButton select



#=====================================================================#
#----------------------------    PACK   ------------------------------#
#=====================================================================#

#  pack transfer radio buttons  #
pack $tranferTypeFrame.transferTypeLabel -pady 5 -padx 5 ; #-anchor c
pack $tranferTypeFrame.imaRadioButton      -anchor w
pack $tranferTypeFrame.bshortsRadioButton  -anchor w
pack $tranferTypeFrame.mincRadioButton     -anchor w


# pack command column #
pack $commandFrame.startButton
pack $commandFrame.stopButton 
pack $commandFrame.viewProgressButton
pack $commandFrame.spacer1
pack $tranferTypeFrame -ipadx 5 -ipady 5



pack .leftFrame -side left -fill both -expand true
pack .commandFrame -side right -anchor n
update

#=====================================================================#
#-------------------------    BINDINGS  ------------------------------#
#=====================================================================#


bind $sessionPick <ButtonRelease-1> \
  { 
    set sessionIndex [ $sessionInfoFrame.sessionList curselection ]
    set sourceDir  [ file dirname $session($sessionIndex,7) ]
#   puts "sessionIndex: $sessionIndex"
#   puts $sourceDir
    $commandFrame.startButton configure -state normal
  }

bind $sessionPick <ButtonRelease-3> \
  { 
    set sessionIndex [ $sessionInfoFrame.sessionList curselection ]
      if { $sessionIndex > -1 } \
	{ CreateSessionInfoView $session($sessionIndex,7)}
  }

bind .leftFrame.targetDirEntryBox <Return> \
    {
	if { [CheckDirOK] } { $commandFrame.startButton configure -state normal } \
         else { $commandFrame.startButton configure -state disabled }
    }

bind all <Control-c> { exit }
bind all <Alt-q> { exit }
bind all <Control-q> { exit }

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
      0  { set targetDir $env(HOME) }

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

# foreach sequenceType [array names knownSequences] \
#    {puts "$sequenceType: $knownSequences($sequenceType)"}

wm title . "NMR Scanner Archive: Retrieve"

if {[QueryDB]} \
    {
        CreateAlertDialog Error "Could not connect to archive database\nAttempting to access secondary index"
	if {[ReadIndexFile $indexFile]} \
	    {
              CreateAlertDialog Error "Could not read secondary index\nAttempting to read incoming directory"
		if {[ReadIncomingDir]} \
	         {
                    CreateAlertDialog Error "Could not read incoming directory"
                    exit
	         }
            }
    }

CreateOutputMonitor
.menubar.mOption entryconfigure 1 -state normal
.menubar.mOption entryconfigure 3 -state normal

#Dialog_Dismiss .outputWindow

#puts "end of program"


#========================================================================#
#------------------------------    TO DO   ------------------------------#
#========================================================================#
set temp " \
Rewrite index maker to pick scan file instead of scout.
Shell script to create an hourly index of incoming directory, including dir size.
Shell script to rcp directories that haven't changed size since last hour.
Regenerate remote archive index.
"
