
##
## unpack_ima.tcl
##
## Original Author: Tony Harris
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

#-------------------- NOTICE ------------------------------------#
# This program is under revision control. Do not edit it without
# going through the proper checkin/checkout steps!
#----------------------------------------------------------------#

# this script looks at the headers of ima files in a predetermined 
# archive directory. The default archive directory is over-ridden 
# by the environment variable ARCHIVE_DIR. the user selects a session 
# and the path to that session is proveded to the script that copies 
# the relevant files via nfs to the local machine and unpacks them 
# locally into b shorts.

# Copyright (C) 2001 Tony Harris tony@nmr.mgh.harvard.edu

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  

#----------------------------------------------------------------#
#Structure of this script:

#global variables
#subroutines
#   dialogs
#   functions
#       session info
#       search
#       transfer
#       utils
#widgets
#bindings
#main
#-----------------------------------------------------------------#

# Alphabetical listing of functions (updated 07 Aug 2001) :

# CheckDirOK
# CheckScratchDirOK
# CheckSourceDirOK
# CheckTargetDirOK
# ConvertMINCfiles
# CreateAlertDialog
# CreateFileBrowser
# CreateMINCfiles
# CreateOutputMonitor
# CreateProgressBarWindow
# CreateSearchDialog
# CreateSequenceDialog
# CreateSessionInfoView
# DeleteOutputMonitor
# DeleteSequenceListItem
# DestroyAlertDialog
# DestroyFileBrowser
# DestroyProgressBarWindow
# DestroySearchDialog
# DestroySequenceDialog
# DestroySessionInfoView
# Dialog_Create
# Dialog_Dismiss
# Dialog_Wait
# DicomPushMinc
# GetHeaderInfo
# GetMincSisterSession
# GetSeqDir
# GetSequences
# GetSessionDescriptors
# GetSessionInfoFromFile
# Ima2sessions
# KillPipe
# LoadListBox
# LocalSearch
# Log
# MountCDROM
# QueryDB
# ReadIncomingDir
# ReadIndexFile
# RestoreSessions
# Run_unpackmincdir
# SaveSessions
# SelectCDROM
# StartTransfer
# TransferIMAfiles
# TransferIMAtoMINC
# TransferSortMincFiles
# ViewCDROM


#----------------------------------------------------------------------------------------#
# global variable defined elsewhere

# sessionNumber: integer 0-N, numerical order of sessions read
#
# session(sessionNumber, fieldValue) 
#   A 2-D array containing a subset of header info from each session.
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
# functions which write to "session" array: QueryDB ReadIncomingDir ReadIndexFile
set numberOfSessionKeys 8

#-----------------------------------------------------------------------------------------#

#define global variables

set targetDir $env(HOME)
#set archiveDir /space/sharbot/1/siemens
set archiveDir /space/bourget/1/siemens
set mincDir /space/bourget/1/minc
set sourceDir $archiveDir
#set indexFile  [file join [file dirname $archiveDir] index2.txt ]
set indexFile  /space/newdata/1/ima_index.txt

set env(FREESURFER_HOME) /space/lyon/1/home/inverse/freesurfer_alpha
set noArchBin "$env(FREESURFER_HOME)/bin/noarch"
set archBin "$env(FREESURFER_HOME)/bin/[exec uname]"
set cdromDir "/mnt/cdrom"

set mincOnly 0
set mincOnlyString ""
set copyDone 1; #copy routines check this after each iteration

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
set sequenceDialogClosed 1; #default state: closed
set searchDialogClosed 1
set lastSessionScanned -1;  # which session is the current sequence info for?
set percentOfSequencesScanned 0
set patientName ""
set afterDate   01;         # search for sessions made after this date
set afterMonth  01
set afterMonthLabel Jan
set afterYear   1990
set currentTime [clock seconds]
set beforeDate  [clock format $currentTime -format %d ]
set beforeMonth [clock format $currentTime -format %m ]
set beforeMonthLabel [clock format $currentTime -format %b ]
set beforeYear  [clock format $currentTime -format %Y ]

# list of stuff to kill when "stop" button is pressed
set mySubProcesses [list ima2mnc ima2mnc.nofork unpackmincdir unpackimadir unpackimadir2 mri_convert dicomserver dicomserver.nofork siemens_dicom_send ]

#--------------------------------  dicom send strings   ----------------------------------#

set dicomSendCommand "${archBin}/siemens_dicom_send"
set mosaicString "-mosaic mni/mosaic64_ 64 64 -mosaic mos64_ 64 64 -mosaic ep2d_fid_ts_20b2604 64 64 -mosaic ep2d_fid_ts_15b3125 64 64 -mosaic ep2d_se_ts_39b2604 64 64 -mosaic ep2d_se_ts_15b3125 64 64 -mosaic ep2d_irm_ts_29b3125 64 64 -mosaic ep2d_irm_ts_39b2604 64 64 -mosaic ep2d_asl_s 64 64"
set portNumber 50082
set archiveHost bourget

#------------------------------    help strings    ---------------------------------------#

set basicHelp "\
Select a session that you want transferred/translated.
Pick a directory to save it to. Hit the start button.

Documentation is at:
http://surfer.nmr.mgh.harvard.edu/archive"

set commandLineHelp "\
This script takes 0 or 1 arguments.
The argument must be a writeable directory."

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
NFS-mounted remote filesystem, you not be able
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
Convert the Siemens IMA data for the selected session 
into the Sessions Format used by FS-FAST. Functional 
data is converted into bshort format, and structural 
data is converted into COR format."

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
To report errors send email to analysis-bugs@nmr.mgh.harvard.edu. 
Include the following: 
(1) Date and approximate time of scan
(2) Scanner used
(3) Subject name
(4) Number of runs and type of each run
(5) destination directory
(6) log file
(7) the nature of the problem."

#-----------------------------------------------------------------------------------------#
# include the progress bar source code

if { [ file exists $noArchBin/progressbar.tcl ] } { source $noArchBin/progressbar.tcl } \
else \
   {
      tk_messageBox -type ok -default ok -title "Error" \
                    -message "could not find  $noArchBin/progressbar.tcl" \
                    -icon error
      exit
   }

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


#---------------------------------     Dialog_Dismiss     ------------------------------#

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



#----------------------------      CreateOutputMonitor       ----------------------------#
# Logging output goes here
# Just write to the pipe "log" and the output appears in this window

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


#-------------------------------     CreateProgressBarWindow     -------------------------------#

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


#---------------------------------    CreateSearchDialog    ----------------------------------#

proc CreateSearchDialog { } \
    {
        global  dialog searchString searchDialogClosed dbType \
                afterMonthLabel beforeMonthLabel \
                patientName afterDate afterMonth afterYear beforeDate beforeMonth beforeYear

        if { ! $searchDialogClosed } { return 1 }
        set sd .searchDialog


        if {[Dialog_Create $sd "Search" -borderwidth 10 ] } \
           {
                             
              set sd .searchDialog

              set nameFrame [ frame $sd.nameBoxFrame ]
              label $nameFrame.searchEntryBoxlabel -text "Patient name"
              entry $nameFrame.nameEntryBox -textvariable patientName -relief sunken 

              set dateFrame  [ frame $sd.dateBoxFrame ]
              set titleLabel [ label $dateFrame.titleLabel -text "" ]
              set dateLabel  [ label $dateFrame.dateLabel  -text "Day" ]
              set monthLabel [ label $dateFrame.monthLabel -text "Month" ]
              set yearLabel  [ label $dateFrame.yearLabel  -text "Year" ]
              grid  $titleLabel $dateLabel $monthLabel $yearLabel

              label $dateFrame.afterLabel -text "After"
              entry $dateFrame.afterDate -textvariable afterDate -relief sunken -width 2 
              #entry $dateFrame.afterMonth -textvariable afterMonth -relief sunken -width 2
              menubutton $dateFrame.afterMonthMenu -textvariable afterMonthLabel \
                                        -menu $dateFrame.afterMonthMenu.menuList
        entry $dateFrame.afterYear -textvariable afterYear -relief sunken -width 4
              grid $dateFrame.afterLabel $dateFrame.afterDate \
                   $dateFrame.afterMonthMenu $dateFrame.afterYear

              label $dateFrame.beforeLabel -text "Before"
              entry $dateFrame.beforeDate -textvariable beforeDate -relief sunken -width 2
              #entry $dateFrame.beforeMonth -textvariable beforeMonth -relief sunken -width 2
              menubutton $dateFrame.beforeMonthMenu -textvariable beforeMonthLabel \
                                                   -menu $dateFrame.beforeMonthMenu.menuList
              entry $dateFrame.beforeYear -textvariable beforeYear -relief sunken -width 4
              grid $dateFrame.beforeLabel $dateFrame.beforeDate \
                   $dateFrame.beforeMonthMenu $dateFrame.beforeYear


              set todayFrame [ frame $sd.todayFrame  -width 50 ]
        set spacer2 [ frame  $todayFrame.spacer2 -height 40 -width 40 ]
              button $todayFrame.todayButton -text today -command {
      set currentTime [clock seconds]
      set beforeDate  [clock format $currentTime -format %d ]
      set beforeMonth [clock format $currentTime -format %m ]
      set beforeMonthLabel [clock format $currentTime -format %b ]
      set beforeYear  [clock format $currentTime -format %Y ]

      set afterDate $beforeDate
      set afterMonth $beforeMonth
      set afterMonthLabel $beforeMonthLabel
      set afterYear $beforeYear
                }


              set spacer1 [ frame $sd.spacer1 -height 20 ]

              

              set buttonFrame [ frame $sd.buttonFrame ]
              button $buttonFrame.cancelButton -text Cancel \
                                               -command { 
                                                           set searchString ""
                                                           DestroySearchDialog
                                                         }
              button $buttonFrame.submitButton -default active -text Search -command \
                  { 
         #if {[ string equal $dbType "sql" ]} { QueryDB } else { LocalSearch }
                  if { [ string compare $dbType "sql" ] == 0 } { QueryDB } else { LocalSearch }            }

  # month menu items #
        set afterMonthMenuList [menu $dateFrame.afterMonthMenu.menuList ]
              $afterMonthMenuList add command -label Jan \
                                              -command { set afterMonthLabel "Jan"; set afterMonth 01}
              $afterMonthMenuList add command -label Feb \
                                              -command { set afterMonthLabel "Feb"; set afterMonth 02}
              $afterMonthMenuList add command -label Mar \
                                              -command { set afterMonthLabel "Mar"; set afterMonth 03}
              $afterMonthMenuList add command -label Apr \
                                              -command { set afterMonthLabel "Apr"; set afterMonth 04}
              $afterMonthMenuList add command -label May \
                                              -command { set afterMonthLabel "May"; set afterMonth 05}
              $afterMonthMenuList add command -label Jun \
                                              -command { set afterMonthLabel "Jun"; set afterMonth 06}
              $afterMonthMenuList add command -label Jul \
                                              -command { set afterMonthLabel "Jul"; set afterMonth 07}
              $afterMonthMenuList add command -label Aug \
                                              -command { set afterMonthLabel "Aug"; set afterMonth 08}
              $afterMonthMenuList add command -label Sep \
                                              -command { set afterMonthLabel "Sep"; set afterMonth 09}
              $afterMonthMenuList add command -label Oct \
                                              -command { set afterMonthLabel "Oct"; set afterMonth 10}
              $afterMonthMenuList add command -label Nov \
                                              -command { set afterMonthLabel "Nov"; set afterMonth 11}
              $afterMonthMenuList add command -label Dec \
                                              -command { set afterMonthLabel "Dec"; set afterMonth 12}

        set beforeMonthMenuList [menu $dateFrame.beforeMonthMenu.menuList ]
              $beforeMonthMenuList add command -label Jan \
                                              -command { set beforeMonthLabel "Jan"; set beforeMonth 01}
              $beforeMonthMenuList add command -label Feb \
                                              -command { set beforeMonthLabel "Feb"; set beforeMonth 02}
              $beforeMonthMenuList add command -label Mar \
                                              -command { set beforeMonthLabel "Mar"; set beforeMonth 03}
              $beforeMonthMenuList add command -label Apr \
                                              -command { set beforeMonthLabel "Apr"; set beforeMonth 04}
              $beforeMonthMenuList add command -label May \
                                              -command { set beforeMonthLabel "May"; set beforeMonth 05}
              $beforeMonthMenuList add command -label Jun \
                                              -command { set beforeMonthLabel "Jun"; set beforeMonth 06}
              $beforeMonthMenuList add command -label Jul \
                                              -command { set beforeMonthLabel "Jul"; set beforeMonth 07}
              $beforeMonthMenuList add command -label Aug \
                                              -command { set beforeMonthLabel "Aug"; set beforeMonth 08}
              $beforeMonthMenuList add command -label Sep \
                                              -command { set beforeMonthLabel "Sep"; set beforeMonth 09}
              $beforeMonthMenuList add command -label Oct \
                                              -command { set beforeMonthLabel "Oct"; set beforeMonth 10}
              $beforeMonthMenuList add command -label Nov \
                                              -command { set beforeMonthLabel "Nov"; set beforeMonth 11}
              $beforeMonthMenuList add command -label Dec \
                                              -command { set beforeMonthLabel "Dec"; set beforeMonth 12}

    #--- pack ----#
              pack $nameFrame.searchEntryBoxlabel $nameFrame.nameEntryBox   -side left
              pack $buttonFrame.cancelButton      $buttonFrame.submitButton -side left
              pack $todayFrame.spacer2            $todayFrame.todayButton   -side left

              pack $nameFrame -fill both -expand true
              pack $dateFrame
              pack $todayFrame
              pack $spacer1
        pack $buttonFrame

   #--- bind ----#
              bind $nameFrame.nameEntryBox <Return> {  
#                 if {[ string equal $dbType "sql" ]} { QueryDB } else { LocalSearch }
                  if { [ string compare $dbType "sql" ] == 0 } { QueryDB } else { LocalSearch }
                 }
        bind $sd <Control-c>       { set searchString ""; DestroySearchDialog }
                      

   #--- main ----#
        set searchDialogClosed 0
        SaveSessions
    }

    }

proc DestroySearchDialog {} \
  { 
    global searchDialogClosed
    RestoreSessions
    LoadListBox
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

              pack $infoFrame -fill both -expand true
              pack $siv.closeButton -anchor center

# -- bind  ---#
              bind $siv <Control-c> DestroySessionInfoView
#              bind $siv <Alt-w>     DestroySessionInfoView

# -- tags  ---#
              $t tag configure bold -font {times 12 bold}

# -- read file  ---#
        if { ! [ file readable $myIMAfile ] } \
                 {
         tk_messageBox -type ok -default ok -title "Error" \
                                   -message "Could not read $myIMAfile" \
                                   -icon error 
                     DestroySessionInfoView
         return 1
                  }

        set infoProg mri_info
        if [string match {*.dcm} $myIMAfile] { set infoProg dcm_info }
              if { [ catch {open "|$archBin/$infoProg $myIMAfile" r} IN ] } \
                   {$t insert end "$IN"; return 1}

       
        set headerLines [ split [ read $IN ] \n ]
        close $IN

#--  print file header in text box  ---#
#         set start [ $t search -count cnt -regexp "\:.*" 1.0 end ]
#               $t tag add bold [expr $start + 1] $cnt

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


#----------------------------       CreateAlertDialog      -----------------------------------#
# general purpose dialog for displaying informationl strings
# some limited ability to expand according to string size, to wrap & stuff like that
proc CreateAlertDialog {title alertMessage} {
  global  dialog session alertDialogClosed archBin
  set ad .alertDialog
  set longestLine 0
  set numberOfLines 30

# -- determine size  ---#
  set textLines [ split $alertMessage \n ]
  foreach line $textLines {
    set lineLength [string length $line ]
    if { $lineLength > 60} { set lineLength 60 }
    if { $lineLength > $longestLine } { set longestLine $lineLength }
  }
  if { $numberOfLines > [llength $textLines ] } { set numberOfLines [llength $textLines ] }

  if {[Dialog_Create $ad $title -borderwidth 10 ] } {
     set ad .alertDialog
     set msgFrame [ frame $ad.msgFrame ]
     set h [expr $numberOfLines + 1]
     if {$h < 6} {set h 6}
     set t [ text $msgFrame.t -setgrid true -wrap word \
             -width [expr $longestLine + 2] -height $h \
             -yscrollcommand "$msgFrame.sy set" ]
     scrollbar $msgFrame.sy -orient vert -command "$msgFrame.t yview"

     button $ad.closeButton -text Close -command { DestroyAlertDialog }

# -- pack  ---#
     pack $msgFrame.sy -side right -fill y
     pack $msgFrame.t -side left -fill both -expand true

     pack $msgFrame -side left -fill both -expand true
     pack $ad.closeButton

# -- bind  ---#
     bind $ad <Control-c> DestroyAlertDialog
     bind $ad <Alt-w>     DestroyAlertDialog
     bind $ad <Return>     DestroyAlertDialog
# -- tags  ---#
     $t tag configure bold -font {times 12 bold}


#--  print msg in text box  ---#


     foreach line $textLines {
        $t insert end "$line\n"
     }
     ::tk::PlaceWindow $ad
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

#         set transferLabelFrame [ frame $tableFrame.transferLabelFrame \
#                                   -relief ridge -borderwidth 2 -bg grey50 \
#               -height $rowHeight -width $columnWidth ]
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
#           -command "if \${ch.checkbuttonState} \
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
                    { 
      if { [ catch "file mkdir $fbInfo(currentDir)" errorMsg ] } \
          {
        tk_messageBox -type ok -message $errorMsg
          }
                    }
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



#-------------------------------------------------------------------------------------#
#                              GET INFO ROUTINES                                      #
#-------------------------------------------------------------------------------------#


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
       set infoProg mri_info
       if [string match {*.dcm} $myIMAfile] { set infoProg dcm_info }

  if { [ catch {open "|$archBin/$infoProg $myIMAfile" r} IN ] } \
           {puts "$IN"; return 1}
       
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



#------------------------------    ReadIncomingDir    ---------------------------------#
# Extract header info from the first file in each session dir on specified directory
# Load "session" array
# called by main, ViewCDROM,  view local directory, view recent pushes

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
             set infoProg mri_info
             if [string match {*.dcm} $fileName] { set infoProg dcm_info }

             if {[catch {exec $infoProg $fileName | grep sequence } exitStatus]} \
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


#-------------------------------   LocalSearch   -------------------------------#
 # Search the "session" array for $searchString
 # rewrites 'session" and updates listbox
 # is called by CreateSearchDialog
proc LocalSearch {} \
    {
       global  session numberOfSessions numberOfSessionKeys searchString patientName

       set searchString $patientName
       set x 0

  for {set listIndex 0} {$listIndex < $numberOfSessions} {incr listIndex} \
         {
           
           #if { [ string match "*$searchString*" $session($listIndex,0) ] } 
           if { [regexp -nocase $searchString $session($listIndex,0) ] } \
         {
                 set matchedSessions($x)  $listIndex
                 incr x
                 lappend hitList $listIndex
                 #puts "$session($listIndex,0)"
              }
         }

       if { ! $x } \
           { 
              tk_messageBox -type ok \
                            -default ok -title "debug" \
                            -message "no match found for \"$searchString\""
              return 1
           }

       update

       # reset session array to just the hits

       set numberOfSessions $x

  for {set id 0} { $id < $numberOfSessions} {incr id } \
     {
         for {set field 0} {$field < $numberOfSessionKeys} {incr field } \
     {
        set session($id,$field) $session($matchedSessions($id),$field)
           }
     }
       update
       LoadListBox


    }

#-----------------------   GetSequences   --------------------------#
# fills sequenceFiles (an array of lists containing filenames)
# called by CreateSequenceDialog

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

           set infoProg mri_info
           if [string match {*.dcm} $fullFileName] { set infoProg dcm_info }

           if {[catch {exec $archBin/$infoProg $fullFileName}   exitStatus]} \
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
#            if {[catch "open /home/tony/Dicom/tempdir/$dirName w" OUT ]} \
#          {puts $OUT; continue}
#            foreach fileName $sequenceFiles($dirName) {puts $OUT $fileName}
#            close $OUT
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
        lappend auto_path $env(FREESURFER_HOME)/local/bin/[exec uname]
        if {[catch {package require Sql} errorMsg]} \
           {
              CreateAlertDialog Error $errorMsg
              return 1
            }

       # Connect to the database
  if {[catch {set conn [sql connect bourget.nmr.mgh.harvard.edu query]} errorMsg]} \
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

proc MountCDROM {} \
{
  global  cdromDir
  set foundFlag 0

# get list of mounted FSs
   set mountOutput [ split [ exec mount ] \n ] 

# is $cdromDir already mounted?
   foreach mountString $mountOutput \
     { 
        if { [ string match  *${cdromDir}* $mountString ] } {  set foundFlag 1 }
     }

# if mounted then return true
    if { $foundFlag } {return 1 }


# try to mount cdrom
     if { [ catch "exec mount $cdromDir" errorMsg ] } \
        {
     tk_messageBox -type ok -default ok -title "Error" \
                   -message "$errorMsg" -icon error 
     return 0
   }

    return 1 ;#mount state is true
}


#--------------------------------------------------------------------------------------#

proc SelectCDROM {} \
{
   global  cdromDir  newDirSelected fbInfo

   set newDirSelected 0
   set initialDir "/"

    if { [ file isdirectory $cdromDir ] } { set initialDir $cdromDir}

   .menubar.mFile entryconfigure 1 -state disabled
   .menubar.mFile entryconfigure 2 -state disabled
   .menubar.mFile entryconfigure 3 -state disabled

   CreateFileBrowser "CDROM" $initialDir
   tkwait variable newDirSelected

     DestroyFileBrowser
    .menubar.mFile entryconfigure 1 -state normal
    .menubar.mFile entryconfigure 2 -state normal
    .menubar.mFile entryconfigure 3 -state normal    
#ok
   if { $newDirSelected == 1 } \
      { set cdromDir $fbInfo(currentDir) }
# cancel
    if { $newDirSelected == 2 } { return 0 ;# not selected }


# try to mount
    if { [MountCDROM] } {return 1 } else {return 0}
# return true if selected and mounted
}

#--------------------------------------------------------------------------------------#

proc ViewCDROM {} \
{
  global  cdromDir archiveDir 
 
 

#is $cdromDir already mounted?
   if { ! [MountCDROM] } \
      {
        # user selects new cdrom Dir
         if { ! [SelectCDROM]} { return 1}
      }

    set archiveDir $cdromDir
    ReadIncomingDir
    return 0

}

#--------------------------------------------------------------------------------------#
proc CheckDirOK {} \
    {
  if { ! [ CheckTargetDirOK ] } { return 0 }
  if { ! [ CheckSourceDirOK ] } { return 0 }
  if { ! [ CheckScratchDirOK ] } { return 0 }
        return 1; # return true if okay
    }

#--------------------------------------------------------------------------------------#
proc CheckTargetDirOK {} \
    {
       global targetDir 

       if { ! [ file exists $targetDir ]} \
           { 
              if { [ file writable [ file dirname $targetDir ] ] } \
                 { file mkdir $targetDir } \
              else \
                {
                  tk_messageBox -type ok -default ok -title "Error" \
                       -message "could not create $targetDir" -icon error 
        set targetDir [pwd]
                  return 0
                }
           }
# targetDir exists, is it a directory?     
  
       if { ![ file isdirectory $targetDir ] } \
          {
             tk_messageBox -type ok -default ok -title "Error" \
                       -message "$targetDir is not a directory" -icon error 
             set targetDir [pwd]
             return 0
          }

# targetDir is a directory, is it empty?

             if { [ string compare [ exec ls $targetDir ] ""] } \
                 {
                     set overwriteAnswer [ tk_messageBox -type yesno \
                                             -default yes \
                                             -title "$targetDir not empty" \
                                             -message "Overwrite contents of $targetDir" \
                                             -icon question ]
                     if { ! [string compare $overwriteAnswer no] }  { return 0  }
                        
                  }


      if { ! [ file writable $targetDir ] } \
              {
                  tk_messageBox -type ok -default ok -title "Error" \
                       -message "$targetDir is not writeable" -icon error 
                  set targetDir [pwd]
                  return 0
         }

    

  return 1
  }

#----------------------------------------------------------------------------#

proc CheckSourceDirOK {} \
{
  global  sourceDir archiveDir

  if { ! [ file isdirectory $sourceDir ] }\
              {
                 set errmsg \
"$sourceDir is not a directory.  It is possible if this data is over one \
year old that the archive is no longer online.  You will need to either use \
the CD made at the time the data was taken or contact \
help@nmr.mgh.harvard.edu and ask that $sourceDir be recovered from \
the archive tapes"
      CreateAlertDialog {ERROR: missing source} $errmsg
                  #tk_messageBox -type ok -default ok -title "Error" \
                  #     -message $errmsg -icon error 
                  set sourceDir $archiveDir
                  return 0
         }

   if { ! [ file readable $sourceDir ] } \
              {
                  tk_messageBox -type ok -default ok -title "Error" \
                       -message "$sourceDir is not readable" -icon error 
                  set sourceDir $archiveDir
                  return 0
         }

   if { ! ( [ string match {*/Sonata*} $sourceDir ] 
    || [ string match {*/Allegra*} $sourceDir ] 
    || [ string match {*/TRIO*} $sourceDir ] ) } {
  return 1
   } else {
  set errmsg {\
This subject cannot be unpacked through browse-sessions. However it can \
be unpacked using a program called unpacksdcmdir. Run 'unpacksdcmdir \
-help' to get extensive documentation on how to use it. The source \
directory can either be the location of the data on a CD or the location \
of the data in the archive. For the subject you just selected, the location is }
  append errmsg $sourceDir
  append errmsg "\n\n\n\n\n\n\n\n\n\n"
  CreateAlertDialog {ERROR: UNSUPPORTED DATA} $errmsg
  set sourceDir $archiveDir
  return 0
   }
  
   # NEVER REACH HERE

  }

#----------------------------------------------------------------------------#

proc CheckScratchDirOK {} \
{
  global  scratchDir newDirSelected fbInfo

# if scratch space was never set
    if { [string compare $scratchDir "no scratch dir set" ] == 0 } \
      {
           set newDirSelected 0
           .menubar.mFile entryconfigure 1 -state disabled
           .menubar.mFile entryconfigure 2 -state disabled
           .menubar.mFile entryconfigure 3 -state disabled

    CreateFileBrowser "Scratch Space" [exec pwd]
           tkwait variable newDirSelected

             DestroyFileBrowser
            .menubar.mFile entryconfigure 1 -state normal
            .menubar.mFile entryconfigure 2 -state normal
            .menubar.mFile entryconfigure 3 -state normal    
#ok
           if { $newDirSelected == 1 } \
              { set scratchDir $fbInfo(currentDir) }
# cancel
           if { $newDirSelected == 2 } { return 1 }
      }

  if { ! [ file exists $scratchDir ] } \
      {
    if { ! [ file writable [file dirname $scratchDir] ] } \
        {
      tk_messageBox -type ok -default ok -title "Error" \
          -message "Cannot create $scratchDir in [file dirname $scratchDir]" 
                      -icon error 
                  return 0
        }
          file mkdir $scratchDir
          
      }

  if { ! [ file isdirectory $scratchDir ] } \
          {
            
             tk_messageBox -type ok -default ok -title "Error" -icon error \
                                -message "$scratchDir is not a directory" 
             return 0
            
          }

  if { ! [ file writable $scratchDir ] } \
          {
             tk_messageBox -type ok -default ok -title "Error" -icon error \
                                -message "$scratchDir is not writable" 
             return 0
            
          }

#puts "$scratchDir good to go [ file writable $scratchDir ]"
    return 1

  }

#-------------------------------------------------------------------------------------#
#                              TRANSFER ROUTINES                                      #
#-------------------------------------------------------------------------------------#



proc TransferIMAfiles {destinationDir} \
    {
    global sourceDir copyDone log progressValue

    set dirContents [ exec ls $sourceDir ]
    set numberOfFilesTransferred 0
    set dirLength [ llength $dirContents ]
    set progressValue 0

    CreateProgressBarWindow "Transferring IMA files"
    #puts "dirLength $dirLength"

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



#-----------------------      TransferIMAtoMINC   --------------------------------------#
# creates minc from ima using ima2mnc
# called by StartTransfer

proc TransferIMAtoMINC {destinationDir} \
    {

       global sourceDir copyDone log archBin unpackmincdirPipe commandFrame progressValue

       set imaDirContents [ exec ls $sourceDir ]
       set imaDirLength [ llength $imaDirContents ]
       $log insert end "Source dir: $sourceDir\n"
       $log insert end "Destination dir: $destinationDir\n"
       $log insert end "Converting $imaDirLength ima files to minc\n"
       $log insert end "\n\n"
       $log see end

       if { $copyDone } \
           { 
             $log insert end "ima->minc conversion aborted\n"
             $log see end
             update
             set copyDone 1
       return 1 
           }

  # if minc already exist, skip ima->minc conversion
  set mincSession [ GetMincSisterSession $sourceDir ]
       #puts "$mincSession"
       if { [ string length $mincSession ] > 0 } \
     {
         set mincDirPath "/space/bourget/1/minc/$mincSession"
         set dirContents [ exec ls $mincDirPath ]
         set numberOfFilesTransferred 0
         set dirLength [ llength $dirContents ]
         set progressValue 0
         CreateProgressBarWindow "Transferring directly from minc archive"
         foreach fileName $dirContents \
       {
           if { [ catch "exec cp $mincDirPath/$fileName $destinationDir" errorMsg ] } \
         {
             set response [ tk_messageBox -type yesno -default no -title "Error" \
              -message "$errorMsg\nContinue?" \
              -icon error ]

             if { [ string match $response "no" ] } \
           {
               set copyDone 1
               DestroyProgressBarWindow
               return 0
           }
         }
           incr numberOfFilesTransferred
           set progressValue [ expr $numberOfFilesTransferred*100/$dirLength ]
           update
       }
         DestroyProgressBarWindow
         set copyDone 1
         return 0
     }

       if { ! [ file exists $archBin/ima2mnc ] } \
           {
              CreateAlertDialog "Error" \
                   "Couldn't find $archBin/ima2mnc"
        set copyDone 1
              return 1
           }
#       puts "$archBin/ima2mnc -range 1 3 -host bourget -aetitle bay2fmri -port 50082 $sourceDir $destinationDir"
       #puts "$archBin/ima2mnc  -servparent -nofork $sourceDir $destinationDir"

  set unpackmincdirPipe [open "|$archBin/ima2mnc  -wait $sourceDir $destinationDir" ]

        fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
        fconfigure $unpackmincdirPipe -blocking 0

        vwait pipeState

       $log insert end "\n\n"
       $log see end
       update
       
       set copyDone 1
       return 0
   }

#(\x1a) -----------------------      Run_unpackmincdir   --------------------------------------#
# creates bshorts using unpackmincdir
# source dir is a minc dir
# called by StartTransfer

proc Run_unpackmincdir {sourceDir destinationDir} \
    {

       global  copyDone log noArchBin unpackmincdirPipe commandFrame mincOnlyString

       set dirContents [ exec ls $sourceDir ]
       set dirLength [ llength $dirContents ]
       $log insert end "$dirLength minc scans\n"
       $log insert end "\n\n"
       $log see end

  #puts "cp ${sourceDir}/${fileName} ${destinationDir}/${fileName}"
        if { $copyDone } \
           { 
             $log insert end "unpacking minc aborted\n"
             $log see end
             update
             return 1 
           }

       if { ! [ file exists $noArchBin/unpackmincdir ] } \
     {
          CreateAlertDialog "Error" \
                   "Couldn't find $archBin/unpackmincdir"
             return 1
     }
  #puts "$noArchBin/unpackmincdir -src $sourceDir -targ $destinationDir $mincOnlyString"
  set unpackmincdirPipe [open "|$noArchBin/unpackmincdir -src $sourceDir -targ $destinationDir $mincOnlyString" ]

        fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
        fconfigure $unpackmincdirPipe -blocking 0 -eofchar \x00

        vwait pipeState

       $log insert end "\n\n"
       $log see end
       update
       
       set copyDone 1
       return 0
    }




#-----------------------------   DicomPushMinc   ----------------------------------#
# calls dicom server to translate ima to minc on the dicom server host
# called by StartTransfer
# this function has been deprecated

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

#-----------------------------   TransferSortMincFiles   ----------------------------------#
# copy over ima files, sorted into sequence directories. translate each ima sequence dir to a minc sequence dir using ima2mnc
#called by StartTransfer
proc TransferSortMincFiles {destinationDir} \
    {
       #sequenceFiles:   array of lists. each key is a seqtype. Each list contains filenames
       #sequenceTypes:  list of sesionsequence types like ep2d, mpr, se, etc
       #knownSequences: associative array of sequence types to dirNames

    global sequenceFiles sequenceTypes knownSequences \
           commandFrame sourceDir archBin progressValue \
           lastSessionScanned sessionIndex sequenceDialogClosed \
           log unpackmincdirPipe  pipeState copyDone scratchDir
    
    if { ! [ CheckScratchDirOK ] } {set copyDone 1; return 1 }

    if { $lastSessionScanned != $sessionIndex } \
      { 
         if { [ CreateSequenceDialog $sessionIndex]} \
            { 
                if { $copyDone } {return 1} \
                else { CreateAlertDialog "Error" "Session Scan failed"; return 1 }
             }
         vwait sequenceDialogClosed
       }
# set up temp dir for ima files
    set imaDir "$scratchDir/ima[pid]"
    if {[ catch "file mkdir $imaDir" errormsg ] } \
       { 
          CreateAlertDialog Error $errormsg
          DestroyProgressBarWindow
    set copyDone 1
          return 1 
       }
    
# Create sequence directories
    foreach sequenceName $sequenceTypes \
      {
#       puts "\nsequenceName: $sequenceName"
#       puts "sequenceDir: $knownSequences($sequenceName)"

       set imaSequenceDir  [file join $imaDir $knownSequences($sequenceName) ]
       set destSequenceDir [file join $destinationDir $knownSequences($sequenceName) ]

      # mkdir will make all parent dirs needed
      if { [ catch "file mkdir $destSequenceDir $imaSequenceDir" errormsg ] } \
          { 
             CreateAlertDialog Error $errormsg
             DestroyProgressBarWindow
             return 1 
           }

# get number of files to transfer
    set numberOfFilesTransferred 0
    set progressValue 0
    set dirContents [ exec ls $sourceDir ]
    set numberOfFilesToTransfer [ llength $dirContents ]

    CreateProgressBarWindow "Transferring files"

# copy all the files from the current sequence to the local ima sequence dir   
        foreach fileName $sequenceFiles($sequenceName) \
          {
             if { $copyDone } \
                { 
                   $log insert end "Transfer aborted\n$numberOfFilesTransferred files transferred\n"
                    $log see end
                    DestroyProgressBarWindow
                    update
                    return 1 
                }

             set sourceFilePath  [file join $sourceDir $fileName ]
             #puts "cp $imaFilePath $tempFilePath"
             if {[catch {exec cp $sourceFilePath $imaSequenceDir }  errorMsg ]} \
                {
                  # error!
                   $log insert end $errorMsg
                   $log insert end "\n"
       $log see end
                 } \
              else \
                {
                  # successful copy
                   incr numberOfFilesTransferred       
                   $log insert end "Transferring $sequenceName: $fileName\n"
                   $log see end
               }
              set progressValue [ expr $numberOfFilesTransferred*100/$numberOfFilesToTransfer ]
             update
        }; #end fileName $sequenceFiles($sequenceName)
      }; #foreach sequenceName $sequenceTypes

#check all files transferred
    if { $numberOfFilesToTransfer != $numberOfFilesTransferred } \
  {
     CreateAlertDialog "Error" \
                "Number of files to translate: $numberOfFilesToTransfer
                 Number of files transferred: $numberOfFilesTransferred

                 ima directory has been preserved at:
                 $imaDir"
            DestroyProgressBarWindow
            set copyDone 1
            return 1
  }
# normal termination of transfer
    $log insert end "\n$numberOfFilesTransferred files transferred\nTransfer complete\n"
    $log see end
    update
    DestroyProgressBarWindow
    set progressValue 0

# write minc files to destination session dirs
    CreateProgressBarWindow "Translating to minc"
    foreach sequenceName $sequenceTypes \
      {
       set imaSequenceDir [file join $imaDir $knownSequences($sequenceName) ]
       set mincSequenceDir [file join $destinationDir $knownSequences($sequenceName) ]

#  puts "$archBin/ima2mnc -servparent -nofork $tempFilePath $destFilePath"
  set unpackmincdirPipe [open "|$archBin/ima2mnc $imaSequenceDir $mincSequenceDir" ]

        fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
        fconfigure $unpackmincdirPipe -blocking 0

        vwait pipeState
  incr numberOfFilesTransferred [ llength $sequenceFiles($sequenceName) ]
        set progressValue [ expr $numberOfFilesTransferred*100/$numberOfFilesToTransfer ]
        $log insert end "\n\n"
        $log see end
        update

      }
    if { $numberOfFilesToTransfer != $numberOfFilesTransferred } \
  {
     CreateAlertDialog "Error" \
                "Number of files to translate: $numberOfFilesToTransfer
                 Number of files transferred: $numberOfFilesTransferred

                 ima directory has been preserved at:
                 $imaDir"

  } \
   else  {
          if {[ catch "exec rm -r $imaDir" errormsg] } { CreateAlertDialog Error $errormsg }
  }

     DestroyProgressBarWindow
     set copyDone 1
     return 0
  }; #end of proc 

#-----------------------------   Ima2sessions   ----------------------------------#
#called by StartTransfer
proc Ima2sessions {destinationDir} \
    {
       global copyDone log noArchBin archBin ima2sessionsPipe \
              sourceDir unpackimadir_cmd seqcfgoption
  #puts "cp ${sourceDir}/${fileName} ${destinationDir}/${fileName}"
        if { $copyDone } \
           { 
             $log insert end "Transfer aborted\n"
             $log see end
             update
             return 1 
           }

  # if minc already exist, skip ima->minc conversion
  set mincSession [ GetMincSisterSession $sourceDir ]
  if { [string length $mincSession ] > 0 } \
     {
         set mincDirPath "/space/bourget/1/minc/$mincSession"
         Run_unpackmincdir $mincDirPath $destinationDir
         set copyDone 1
         return 0
     } 


       if { ! [ file exists $noArchBin/$unpackimadir_cmd ] } \
     {
          CreateAlertDialog "Error" \
                   "Couldn't find $noArchBin/$unpackimadir_cmd"
             return 1
     }

       #puts "$noArchBin/unpackimadir  -src $sourceDir -targ $destinationDir"

  set ima2sessionsPipe [open "|$noArchBin/$unpackimadir_cmd -src $sourceDir -targ $destinationDir $seqcfgoption" ]

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
       # Converts individual ima files to minc files directly
       # This function has been deprecated because its dependency, "mri_convert", has a known bug

        global archBin scratchDir
        set imaDir "$scratchDir/ima[pid]"  
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


#-----------------------         ConvertMINCfiles    --------------------------------------#
# performs triage on minc files using unpackmincdir

proc ConvertMINCfiles {mincOnly} \
    {
      global targetDir scratchDir
  set mincDir "$scratchDir/minc[pid]"  
      file mkdir $mincDir
      ConvertIMAfiles $mincDir
      set unpackmincdirPipe [open "|unpackmincdir  -src $mincDir -targ $targetDir $mincOnly" ]

      fileevent $unpackmincdirPipe readable {Log $unpackmincdirPipe}
      fconfigure $unpackmincdirPipe -blocking 0
      vwait pipeState
      
   }




#--------------------------    startTransfer    -------------------------------#
# coordinates the stransfer/translate routines
# is called by the start button
proc StartTransfer {} \
     {   
        global commandFrame transferType pipeState copyDone  mincOnlyString \
               targetDir sourceDir scratchDir

        if { ! [CheckTargetDirOK] } \
           { 
             $commandFrame.startButton configure -state disabled 
             set copyDone 1
             return 1
           }

        if { ! [CheckSourceDirOK] } \
           { 
             $commandFrame.startButton configure -state disabled 
             set copyDone 1
             return 1
           }
        set pipeState 0
        set copyDone 0
        $commandFrame.startButton config -state disabled
        $commandFrame.stopButton config -state normal

        switch -exact $transferType \
          {
              #for bshorts
        experimental {
                              if { ! [CheckScratchDirOK] } \
                                  { 
                                    $commandFrame.startButton configure -state disabled 
                                    set copyDone 1
                                    return 1
                                  }
                  set mincDir "$scratchDir/minc[pid]"
                  file mkdir $mincDir
                  #puts "minc dir: $mincDir"
                  if { [TransferIMAtoMINC $mincDir] == 0 } \
                     { 
             set copyDone 0; #reset global flag
             set mincDirPath [file join $mincDir [ exec ls $mincDir ] ]
             #puts "minc dir: $mincDirPath"; #dicom_send creates a subdir
                                     Run_unpackmincdir $mincDirPath $targetDir 
                                 } \
          else {puts "TransferIMAtoMINC failed"}
      exec rm -rf $mincDir "$scratchDir/ima[pid]"
        }



              minc    { TransferIMAtoMINC $targetDir }

        customMincSort { TransferSortMincFiles $targetDir  }

        mincSorted {
      set mincOnlyString "-minconly"
      set mincSession [ GetMincSisterSession $sourceDir ]
      if { [string length $mincSession ] > 0 } \
          {
        # if minc already exist, skip ima->minc conversion
        set mincDirPath "/space/bourget/1/minc/$mincSession"
        Run_unpackmincdir $mincDirPath $targetDir
          } \
      else {
                              if { ! [CheckScratchDirOK] } \
                                  { 
                                    $commandFrame.startButton configure -state disabled 
                                    set copyDone 1
                                    return 1
                                  }
               
                  set mincDir "$scratchDir/minc[pid]"
                  file mkdir $mincDir
                  #puts "minc dir: $mincDir"
                  if { [TransferIMAtoMINC $mincDir] == 0 } \
                     { 
             set copyDone 0; #reset global flag
             #dicom_send creates a subdir
             set mincDirContents [list [ exec ls $mincDir ] ]
                                     #the dicomserver will create exactly one session dir
             if { [llength $mincDirContents ] == 1 } \
           {
                    set mincDirPath [ file join $mincDir $mincDirContents ]
                    #puts "minc dir: $mincDirPath"; 
                                             Run_unpackmincdir $mincDirPath $targetDir 
           }
                                 } \
          else { puts "TransferIMAtoMINC failed" }
                 exec rm -fr $mincDir "${scratchDir}/ima[pid]"

                           
      }
      set mincOnlyString ""
                          }

             #bshorts { TransferBshortFiles $targetDir }
             bshorts { Ima2sessions $targetDir }

             dicom_minc     { DicomPushMinc $targetDir }

             ima     { TransferIMAfiles $targetDir }

             default { TransferIMAfiles $targetDir }
          }
     
         #tkwait variable copyDone

         $commandFrame.stopButton configure -state disabled
         $commandFrame.startButton configure -state normal 
         #puts "StartTransfer exited"
     }



#-------------------------------------------------------------------------------------#
#                              UTILITY ROUTINES                                      #
#-------------------------------------------------------------------------------------#


proc Log {pipeName} \
    {
  global  pipeState log 
        if [eof $pipeName] { KillPipe $pipeName } \
        else \
      {
               gets $pipeName line
#     if { ! [string compare $line "Done sending files." ] } \
#                    { 
#                        $log insert end "log: $line\n"
#                        $log see end
#                       if {[catch {close $pipeName} errorMsg ]} \
#                          {
#            CreateAlertDialog "Error: Log" $errorMsg
#                    }
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
#         {

#                 if {[catch {exec kill $pid} errorMsg ]} \
#              { 
#                     CreateAlertDialog "Error" "Could not close process $pid\n$errorMsg"
#                     return 1
#        }
#               }
      }
  set pipeState 1
        return 0 
    }


#---------------------------------     SaveSessions    -------------------------------------#
# save the "session" array before starting a search
# called by CreateSearchDialog
proc SaveSessions { } \
    {
        global session numberOfSessions numberOfSessionKeys savedSessions numberOfSavedSessions

        set numberOfSavedSessions $numberOfSessions

  for {set id 0} { $id < $numberOfSessions} {incr id } \
     {
         for {set field 0} {$field < $numberOfSessionKeys} {incr field } \
     {
        set savedSessions($id,$field) $session($id,$field) 
           }
     }


    }

proc RestoreSessions { } \
    {
        global session numberOfSessions numberOfSessionKeys savedSessions numberOfSavedSessions

        set numberOfSessions $numberOfSavedSessions

  for {set id 0} { $id < $numberOfSessions} {incr id } \
     {
         for {set field 0} {$field < $numberOfSessionKeys} {incr field } \
     {
        set session($id,$field) $savedSessions($id,$field) 
           }
     }


    }

#---------------------------------     GetMincSisterSession    --------------------------------#
# find out if there is a minc session that corresponds to a given ima session

proc GetMincSisterSession { imaDir } \
    {
        global mincDir

  set baseIMAdirName [ file tail $imaDir ]
  set dirNameSegments [ split $baseIMAdirName - ]
  set scanNumberString [ lindex $dirNameSegments end ]
  set scanNumberStringLength [ expr [string length $scanNumberString ] + 2 ]
  set baseIMAdir [ string range $baseIMAdirName 0  [expr [string length $baseIMAdirName ] - $scanNumberStringLength ]  ] 
        #puts "$baseIMAdir"
        
  #set baseMincDirname [ exec ls -d $mincDir/\*$baseIMAdir ]; # doesn't work
  set mincSessions   [ split [ exec ls -1 $mincDir ] "\n" ]
  foreach mincSession $mincSessions \
  {
      if { [ string match *$baseIMAdir $mincSession ] } \
      {
    #puts "found $mincSession"
    return $mincSession
      }
   }
  return ""
    }




#=============================================================================#
#---------------------------- Define widgets   -------------------------------#
#=============================================================================#


menu .menubar
#attach it to the main window
. config -menu .menubar
foreach menuWidget { File View Search Option Help } \
    {
        set $menuWidget [ menu .menubar.m$menuWidget ]
        .menubar add cascade -label $menuWidget -menu .menubar.m$menuWidget 
    }

$File add command -label "Select Archive" -command \
        { 
          set newDirSelected 0
    .menubar.mFile entryconfigure 1 -state disabled ; #disable other file browser items
          .menubar.mFile entryconfigure 2 -state disabled
          .menubar.mFile entryconfigure 3 -state disabled

           CreateFileBrowser "Archive" $archiveDir
           tkwait variable newDirSelected

         if { $newDirSelected == 1 } \
            { set archiveDir $fbInfo(currentDir) }

          DestroyFileBrowser
         .menubar.mFile entryconfigure 1 -state normal
         .menubar.mFile entryconfigure 2 -state normal
         .menubar.mFile entryconfigure 3 -state normal

       if { ! [CheckSourceDirOK] } { $commandFrame.startButton configure -state disabled }

       }
                       

$File add command -label "Select Destination" -command \
    { 
      set newDirSelected 0
      .menubar.mFile entryconfigure 1 -state disabled
      .menubar.mFile entryconfigure 2 -state disabled
      .menubar.mFile entryconfigure 3 -state disabled

      CreateFileBrowser "Destination" $targetDir
      tkwait variable newDirSelected

         if { $newDirSelected == 1 } \
            { set targetDir $fbInfo(currentDir) }

               if { ! [CheckTargetDirOK] } { $commandFrame.startButton configure -state disabled }

             }
         DestroyFileBrowser
         .menubar.mFile entryconfigure 1 -state normal
         .menubar.mFile entryconfigure 2 -state normal
         .menubar.mFile entryconfigure 3 -state normal


$File add command -label "Select Scratch Dir" -command \
       { 
           set newDirSelected 0
           set oldScratchDir $scratchDir

           .menubar.mFile entryconfigure 1 -state disabled
           .menubar.mFile entryconfigure 2 -state disabled
           .menubar.mFile entryconfigure 3 -state disabled

           if { [string compare "no scratch dir set" $scratchDir ] == 0 } \
               { CreateFileBrowser "Scratch Space" [pwd] } \
           else \
               { CreateFileBrowser "Scratch Space" $scratchDir }

           tkwait variable newDirSelected

           if { $newDirSelected == 1 } \
              { set scratchDir $fbInfo(currentDir) }


             DestroyFileBrowser
            .menubar.mFile entryconfigure 1 -state normal
            .menubar.mFile entryconfigure 2 -state normal
            .menubar.mFile entryconfigure 3 -state normal

     CheckScratchDirOK
      
       }

$File add command -label "Select CDROM Dir" -command { SelectCDROM }

$File add command -label Quit -command {  exit   }


         #--------- View ------------#


$View add command -label "cdrom"   -state normal \
                                       -command { ViewCDROM  }

$View add command -label "Local Directory"   -state normal \
                                           -command { 
#      if { [ string equal $archiveDir "/space/bourget/1/siemens" ] } 
      if { [ string compare $archiveDir "/space/bourget/1/siemens" ] == 0 } \
    { 
              set newDirSelected 0
              .menubar.mFile entryconfigure 1 -state disabled
              .menubar.mFile entryconfigure 2 -state disabled

              CreateFileBrowser "Archive" $archiveDir
              tkwait variable newDirSelected

              if { $newDirSelected == 1 } \
                 { set archiveDir $fbInfo(currentDir) }

              if { [CheckSourceDirOK] } \
      { ReadIncomingDir; $commandFrame.startButton configure -state normal } \
              else { $commandFrame.startButton configure -state disabled }


               DestroyFileBrowser
               .menubar.mFile entryconfigure 1 -state normal
               .menubar.mFile entryconfigure 2 -state normal
    } \
     else { ReadIncomingDir }
 
   }



$View add command -label "Recent Pushes"   -state normal \
                                           -command { set archiveDir /space/bourget/1/siemens
                                                      ReadIncomingDir  }

$View add command -label "Entire database"   -state normal \
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
                                                  vwait sequenceDialogClosed
                                              }

         #--------- Search ------------#

$Search add command -label "Archive database"   -state normal \
                    -command { 
                                .menubar.mSearch entryconfigure 1 -state disabled
                                .menubar.mSearch entryconfigure 2 -state disabled
                               set dbType sql
                               CreateSearchDialog
                                .menubar.mSearch entryconfigure 1 -state normal
                                .menubar.mSearch entryconfigure 2 -state normal
                             }

$Search add command -label "Window Contents"   -state normal \
                    -command { 
                                .menubar.mSearch entryconfigure 1 -state disabled
                                .menubar.mSearch entryconfigure 2 -state disabled
              set dbType local
                                CreateSearchDialog 
                                .menubar.mSearch entryconfigure 1 -state normal
                                .menubar.mSearch entryconfigure 2 -state normal

                             }

         #--------- Option ------------#

$Option add command -label "Load text DB"   -state disabled \
                                       -command { 
                                                  set copyDone 1
                                                  ReadIndexFile $indexFile
                                                }

$Option add command -label "state info" -state normal \
                                        -command { 
                                              CreateAlertDialog "State" \
                 "sessions: $numberOfSessions\nArchive Dir: $archiveDir\nDest Dir: $targetDir\nScratch Dir: $scratchDir\ncdrom Dir: $cdromDir"
 }

$Option add command -label "minc sorted"   -state normal \
                                       -command { 
                                                  set transferType customMincSort
                    StartTransfer
                                                }

$Option add cascade  -label "do not use" -menu $Option.experimental
    set expMenu [ menu $Option.experimental -tearoff 1 ]

$expMenu add command -label "find minc"   -state normal \
                                       -command { 
                                                  GetMincSisterSession $sourceDir
                                                }

$expMenu add command -label "bshorts"   -state normal \
                                       -command { 
                                                  set transferType experimental
                    StartTransfer
                                                }

$expMenu add command -label "kill pipe"   -state normal \
                                       -command { 
                                                  KillPipe $unpackmincdirPipe
                                                }


         #--------- Help ------------#

$Help add command -label "Basic Help" -command { CreateAlertDialog "Basic Help"  $basicHelp }


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



$Help add cascade -label "Output Format" -menu $Help.transfer

set TransferHelp [menu $Help.transfer -tearoff 0 ]

$TransferHelp add command -label "Sessions" \
                  -command { CreateAlertDialog "b short" $transferBshortHelp }

$TransferHelp add command -label "IMA" \
                  -command { CreateAlertDialog "IMA" $transferIMAhelp }

$TransferHelp add command -label "minc" \
                  -command { CreateAlertDialog "minc" $transferMINChelp }

# $TransferHelp add command -label "dicom minc" \
#                   -command { CreateAlertDialog "dicom minc" $transferdicomMinchelp }



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

          foreach subProcess $mySubProcesses \
      {
              catch "exec killall $subProcess" errorMsg
    #puts $errorMsg
            }
                                 
                                 $commandFrame.stopButton config -state disabled
                                 $commandFrame.startButton config -state normal

                                 }

button $commandFrame.viewProgressButton -text "view log" \
                                -width $commandFrameButtonWidth \
                                -command CreateOutputMonitor

frame       $commandFrame.spacer1 -height 15

set tranferTypeFrame [ frame $commandFrame.tranferTypeFrame -borderwidth 2 -relief ridge  ]

label       $tranferTypeFrame.transferTypeLabel  -text "Output Format"
radiobutton $tranferTypeFrame.bshortsRadioButton -text "Sessions\n(bshort)" -value bshorts -variable transferType -state normal
radiobutton $tranferTypeFrame.mincRadioButton    -text "MINC" -value minc -variable transferType -state normal
radiobutton $tranferTypeFrame.mincSortedRadioButton    -text "MINC (sorted)" -value mincSorted -variable transferType -state normal
radiobutton $tranferTypeFrame.imaRadioButton     -text "IMA"   -value ima  -variable transferType

$tranferTypeFrame.bshortsRadioButton select



#=====================================================================#
#----------------------------    PACK   ------------------------------#
#=====================================================================#

#  pack transfer radio buttons  #
pack $tranferTypeFrame.transferTypeLabel -pady 5 -padx 5 ; #-anchor c
pack $tranferTypeFrame.bshortsRadioButton  -anchor w
pack $tranferTypeFrame.imaRadioButton      -anchor w
pack $tranferTypeFrame.mincRadioButton     -anchor w
pack $tranferTypeFrame.mincSortedRadioButton     -anchor w


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

bind all <Alt-f> {
                    .menubar.mSearch entryconfigure 1 -state disabled
                    .menubar.mSearch entryconfigure 2 -state disabled
                     set dbType sql
                     CreateSearchDialog
                     .menubar.mSearch entryconfigure 1 -state normal
                     .menubar.mSearch entryconfigure 2 -state normal
                 }
bind all <Alt-t> {set patientName "xfer"; QueryDB }

bind all <Control-c> { exit }
bind all <Alt-q> { exit }
bind all <Control-q> { exit }

#=====================================================================#
#-----------------------------    MAIN   -----------------------------#
#=====================================================================#

#----------------------- COMMAND LINE STUFF ---------------------------------#

set targetDir .
set unpackimadir_cmd unpackimadir
set seqcfgoption ""

# Parse Command line options #
set ntharg 0
while  { $ntharg < $argc } {

  set option [lindex $argv $ntharg];
  #puts "option $ntharg $option"
  incr ntharg 1

  # number or arguments remaining
  set nargrem [expr $argc - $ntharg]; 

  switch -exact -- $option {
    -unpackimadir2 {
       set seqcfgfile "$env(FREESURFER_HOME)/scanseq.unpackcfg"
       set unpackimadir_cmd unpackimadir2
       set seqcfgoption "-seqcfg $seqcfgfile"
       #puts stdout "SeqCfgFile $seqcfgfile"
    }
    -unpackimadir {
       set unpackimadir_cmd unpackimadir
       set seqcfgoption ""
    }
    default { 
       set targetDir $option
    }
  }
}

#----------------------- DEPENDENCIES STUFF ---------------------------------#

if {! [file readable $archiveDir ] } \
    { 
       tk_messageBox -type ok -default ok -title "Error" \
                     -icon error \
                         -message "Could not find archive directory $archiveDir"
       exit
    }

if { [file writable "/scratch" ] }  { set scratchDir "/scratch" } \
    else { set scratchDir "no scratch dir set" }

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
