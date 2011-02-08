-- A common file for the MRIconvolve1d tests

module( ..., package.seeall )


require( "tools.freesurfer" )

-- Make access to the FreeSurfer module quicker
local FS = tools.freesurfer



-- Kernel sizes to use
kernels = { 4, 15, 22 }
-- Directions, as specified in mri.h
directions = { x=1, y=0, z=2 }



-- Function to return the list of files to use
function GetFiles( pattern )
   return FS.GetFiles( FS.MRIdir(), pattern )
end

-- Generate a test name from the filename, kernel size and convolution direction
local function testName( input, size, direction )
   local inputItems = FS.split( input, "-" )

   return inputItems[1].."-"..inputItems[2].."-"..size.."-"..direction
end


-- Generates one test configuration for each input, kernel size and convolution direction
function testGen( inputs, sizes, directions, tol )
   local testTable = {}

   for i,input in ipairs(inputs) do
      for s,size in ipairs(sizes) do
	 for dirName,dirVal in pairs(directions) do
	    local tName = testName( input, size, dirName )
	    table.insert( testTable, { id=tName, tol=tol, input=FS.MRIdir()..input, kernel=size, direction=dirVal } )
	 end
      end
   end

   return testTable
end


-- Script to run
script = [[
      $(submit JOBNAME="$(id)", TIME="00:15:00" )

      export MAINCMDS=" --input=$(input) --kernelSize=$(kernel) --direction=$(direction) --repeats=1"
      export CPUOUT="$(outputDir)/$(id).cpu.mgz"
      export GPUOUT="$(outputDir)/$(id).gpu.mgz"

      echo $MAINCMDS

      ${TM_BIN_DIR}/conv1d_test $MAINCMDS \--output=$CPUOUT
      ${TM_BIN_DIR}/conv1d_test_cuda $MAINCMDS \--output=$GPUOUT

      $(projectDir)/tools/mridiff.pl \--results=$(cmdResultFn) \
	                             \--gold=$CPUOUT \
                                     \--test=$GPUOUT \
                                     \--threshold=$(tol)
      
      testFinish -c $(cmdResultFn) -r $(resultFn) -t $(runtimeFn)
]]