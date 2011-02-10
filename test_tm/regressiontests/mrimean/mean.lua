-- -*- lua -*-


-- A common file for the MRImean tests

module( ..., package.seeall )


require( "tools.freesurfer" )

-- Make access to the FreeSurfer module quicker
local FS = tools.freesurfer


-- Filter widths to use
widthVals = { 1, 3, 17 }


-- Function to return the list of files to use
function GetFiles( pattern )
   return FS.GetFiles( FS.MRIdir(), pattern )
end



-- Function to turn an input filename and mean width into a test name
local function testName( input, width )
   local inputItems = FS.split( input, "-" )

   return inputItems[1].."-"..inputItems[2].."-"..width
end



-- Function to generate a list of tests
function testGen( inputs, widths, tol )
   local testTable = {}

   -- inputs is a table of input file names
   -- widthss is a table of mean filter widths
   -- tol is the tolerance for mri_diff
   -- Generates nInputs*nWidths tests

   for i,input in ipairs(inputs) do
      for w,width in ipairs(widths) do
	 local tName = testName( input, width )
	 table.insert( testTable, { id=tName, tol=tol, input=FS.MRIdir()..input, width=width } )
      end
   end

   return testTable
end



-- Script to run
script = [[
      $(submit JOBNAME="$(id)", TIME="01:00:00" )

      export MAINCMDS=" --input=$(input) --width=$(width) --repeats=1"
      export CPUOUT="$(outputDir)/$(id).cpu.mgz"
      export GPUOUT="$(outputDir)/$(id).gpu.mgz"

      ${TM_BIN_DIR}/meanfilter_test $MAINCMDS \--output=$CPUOUT
      ${TM_BIN_DIR}/meanfilter_test_cuda $MAINCMDS \--output=$GPUOUT

      $(projectDir)/tools/mridiff.pl \--results=$(cmdResultFn) \
                                     \--gold=$CPUOUT \
                                     \--test=$GPUOUT \
                                     \--threshold=$(tol)

      testFinish -c $(cmdResultFn) -r $(resultFn) -t $(runtimeFn)
]]