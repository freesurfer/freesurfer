
##
## browse-minc-header.tcl
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

set binDir /space/beowulf/1/users/inverse/freesurfer_alpha/bin/noarch
#set binDir /home/tony/MINC
set mincFile "Browse Minc Headers"
set currentDir .
set typelist { \
               {"MINC file" ".mnc" "MINC" } \
               {"all" "*" {} } \
              }



#====================================================#
#--------------  getMincVariables  -----------------#
#====================================================#
proc getMincVariables {} \
 {
   global mincFile .listBoxFrame.leftList .listBoxFrame.rightList mincFile
   .listBoxFrame.leftList delete  0 end
   .listBoxFrame.rightList delete  0 end
   foreach mincvariables [exec mincinfo $mincFile -varnames]\
       {
	   .listBoxFrame.leftList insert end $mincvariables
       }
    
 }


#----------------------------------------------------#

#====================================================#
#--------------  getMincAttributes  -----------------#
#====================================================#
proc getMincAttributes { } \
 {
  # A variable has been selected. Get list of attributes
  # and return the name of the variable.
 
     global mincFile outputText .listBoxFrame.rightList .listBoxFrame.leftList
     .listBoxFrame.rightList delete  0 end
     set variableIndex [ .listBoxFrame.leftList curselection ]
     

     set variable [ .listBoxFrame.leftList get $variableIndex ]
          
     foreach mincVariable [ exec mincinfo $mincFile -varatts $variable ] \
       {
	   .listBoxFrame.rightList insert end $mincVariable
       }
    return $variable
 }


#====================================================#
#-----------------  getMincValue  -------------------#
#====================================================#
proc getMincValue { variable  } \
  {
    # An attribute has been selected. Print the value and
    # return the name of the attribute

    global mincFile outputText .listBoxFrame.rightList .listBoxFrame.leftList
    

    set attributeIndex [ .listBoxFrame.rightList curselection ]
    set attribute [ .listBoxFrame.rightList get $attributeIndex ]
    set outputText [ exec  mincinfo $mincFile -attvalue $variable:$attribute ]
    

    return $attribute
  }

#====================================================#
#---------------  getMincFile  -----------------#
#====================================================#
proc getMincFile {   } \
  {
    # copy the browsed file to a different location
    
    global mincFile currentDir 
    
    set baseMincFileName [ file tail $mincFile ]
    set destinationFile [ tk_getSaveFile -title "Get MINC file" \
                                    -defaultextension .mnc     \
                                    -initialdir $currentDir     \
		                    -initialfile $baseMincFileName ]
    
    if {! [ string compare $destinationFile "" ] }  { return ""}
    set currentDir [ file dirname $mincFile ]
#    if {[ file exists $destinationFile ]} \
\#	{
#	   set overwriteAnswer [tk_messageBox -type okcancel -default cancel \
\#                             -title "Overwrite?" \
\#                             -message "overwrite existing $destinationFile?"\
\#				    -icon warning ]
#	    if { [string compare $overwriteAnswer "ok"] } {return ""}  
#	}
     file copy -force $mincFile $destinationFile

    }

#====================================================#
#---------------  getMincDir  -----------------#
#====================================================#
proc getMincDir {   } \
  {
    # copy the browsed file to a different location
    
    global mincFile currentDir outputText
    
    set fileNames [ exec ls $currentDir ]
    set baseMincDir [ file tail $currentDir ]
    set destinationDirName [ tk_getSaveFile -title "Copy MINC Dir" \
                                    -initialdir $currentDir  ]
    
    if {! [ string compare $destinationDirName "" ] }  { return ""}
    
    set destinationDirParent [ file dirname $destinationDirName ]
    if { ! [ file writable $destinationDirParent ] } \
       { 
         tk_messageBox -type ok -default ok -title "Error" \
              -message "No write permission in $destinationDirParent " \
              -icon error 
         return ""
       }

    file mkdir $destinationDirName

     foreach fileName $fileNames \
      {  
        #set outputText $destinationDirName/$fileName 
        file copy  $currentDir/$fileName $destinationDirName/$fileName
      }

      set currentDir [ file dirname $destinationDirName ]

    }


#====================================================#
#--------------  Create the command bar -------------#
#====================================================#



frame .commandFrame 

frame .commandFrame.buttonFrame 
button .commandFrame.buttonFrame.fileButton \
    -text "View File" -width 8 \
    -command { set newFile [ tk_getOpenFile -title "Select MINC file" \
                                             -filetypes $typelist \
                                             -initialdir $currentDir ] 

               if { [ string compare $newFile ""] } { set mincFile $newFile }
               set currentDir [ file dirname $mincFile ]
               wm title . $mincFile
               getMincVariables 
              }

button .commandFrame.buttonFrame.getFileButton -text "Copy File" -width 8 \
	-command { getMincFile}

button .commandFrame.buttonFrame.getDirButton -text "Copy Dir" -width 8 \
	-command { getMincDir }

#---------- Create the output frame ---------------#

frame .commandFrame.outputFrame
label .commandFrame.outputFrame.variableBox -textvariable mincVariable
label .commandFrame.outputFrame.attributeBox -textvariable mincAttribute
label .commandFrame.outputFrame.outputBox \
                   -width 30 -wraplength 200 -font roman \
                   -textvariable outputText


#------------------  Create Exit button -----------------#


  if { [file readable $binDir/Susan.gif] }\
     {
        set susan [image create photo \
          -file $binDir/Susan.gif -format gif \
          -height 120 -width 130 -palette 256 ]

         set exitIcon [ image create photo ]
         $exitIcon copy $susan -subsample 2 
         set exitIconType "image"
     } \
  else { set exitIconType "text"
         set exitIcon "Quit"
        }

button .commandFrame.susanButton -$exitIconType $exitIcon -command exit


#====================================================#
# ----------------  Create list boxes ---------------#
#====================================================#

frame .labelFrame
label .labelFrame.variableLabel -width 20  -height 2 \
                                 -justify left -text "\nVariables"
label .labelFrame.attributeLabel -width 20 -height 2 \
                                 -justify right -text "\nAttributes"


frame .listBoxFrame

scrollbar .listBoxFrame.leftScroll -command ".listBoxFrame.leftList yview"
set variablePick [listbox .listBoxFrame.leftList \
                  -yscroll ".listBoxFrame.leftScroll set" \
		      -relief sunken -width 20 -height 15 -setgrid yes ]

scrollbar .listBoxFrame.rightScroll -command ".listBoxFrame.rightList yview"
set attributePick [ listbox .listBoxFrame.rightList  \
                 -yscroll ".listBoxFrame.rightScroll set" \
                 -relief sunken -width 20 -height 15 -setgrid yes ]


#====================================================#
# -------------  Pack main canvas widgets -----------#
#====================================================#
pack .listBoxFrame -side bottom

pack .listBoxFrame.rightScroll .listBoxFrame.rightList  -side right  \
        -fill both -expand yes -anchor w
pack .listBoxFrame.leftList .listBoxFrame.leftScroll  -side left   \
        -fill both -expand yes -anchor w

pack .labelFrame -side bottom
pack .labelFrame.variableLabel -side left
pack .labelFrame.attributeLabel -side left

pack .commandFrame -side top


pack .commandFrame.buttonFrame -side left
pack .commandFrame.buttonFrame.fileButton 
pack .commandFrame.buttonFrame.getFileButton 
pack .commandFrame.buttonFrame.getDirButton 

pack .commandFrame.outputFrame -side left
pack .commandFrame.outputFrame.variableBox
pack .commandFrame.outputFrame.attributeBox
pack .commandFrame.outputFrame.outputBox -side left

pack .commandFrame.susanButton -side right


#====================================================#
#----------------    BINDINGS  ----------------------#
#====================================================#

bind $variablePick <ButtonRelease-1> \
  { 
    set mincVariable [ getMincAttributes  ]
    set mincAttribute ""
    set outputText ""
  }
bind $attributePick <ButtonRelease-1>  \
	{ set mincAttribute [ getMincValue $mincVariable  ] }



# Set up bindings for the file manager.  Control-c closes the window.
bind all <Control-c> {destroy .}


#====================================================#
#-------------------    MAIN   ----------------------#
#====================================================#



#--------------- check for mincinfo   ------------------#


catch {exec which mincinfo} exitStatus
if { ! [file executable $exitStatus ] } \
 {
    tk_messageBox -type ok -default ok -title "Error" \
                       -message $exitStatus \
                       -icon error
    exit
  }


#--------------- check minc file   ------------------#

if { $argc > 1 } \
  {
    tk_messageBox -type ok -default ok -title "Error" \
                       -message "usage:  bmh.tcl \[mincfile\]" -icon warning
    exit
  }


if { $argc == 1 } \
  {
    set mincFile [ lindex $argv 0 ]
    if ![ file readable $mincFile ] \
      {
         tk_messageBox -type ok -default ok -title "Error" \
                       -message "$mincFile not readable" -icon error 
         
         set mincFile [ tk_getOpenFile  -title "Select MINC file"\
                        -filetypes $typelist \
                        -initialdir $currentDir ]
      }
  } \
else \
    { 
       set mincFile [ tk_getOpenFile -title "Select MINC file" \
                        -filetypes $typelist \
                        -initialdir $currentDir ]
     }

if { ! [ string compare $mincFile "" ] } { exit }


#------------------- body   ----------------------#

set currentDir [ file dirname $mincFile ]
wm title . $mincFile
getMincVariables

