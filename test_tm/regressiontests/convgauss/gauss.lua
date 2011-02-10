-- A common file for the MRIconvolveGaussian tests

module( ..., package.seeall )


require( "tools.freesurfer" )

-- Make access to the FreeSurfer module quicker
local FS = tools.freesurfer



-- List of gaussian widths to use
sigmaVals = { 0.5, 0.79788, 1.59597, 2, 6.383, 12.766, 29.2 }

-- Function to return the list of files to use
function GetFiles( pattern )
   return FS.GetFiles( FS.MRIdir(), pattern )
end



-- Function to turn an input filename and gaussian sigma into a test name
local function testName( input, sigma )
   local inputItems = FS.split( input, "-" )

   return inputItems[1].."-"..inputItems[2].."-"..sigma
end



-- Function to generate a list of tests
function testGen( inputs, sigmas, tol )
   local testTable = {}

   -- inputs is a table of input file names
   -- sigmas is a table of gaussian sigmas
   -- tol is the tolerance for mri_diff
   -- Generates nInputs*nSigmas tests

   for i,input in ipairs(inputs) do
      for s,sigma in ipairs(sigmas) do
	 local tName = testName( input, sigma )
	 table.insert( testTable, { id=tName, tol=tol, input=FS.MRIdir()..input, sigma=sigma } )
      end
   end

   return testTable
end



-- Script to run
script = [[
      $(submit JOBNAME="$(id)", TIME="00:15:00" )


      export MAINCMDS=" --input=$(input) --sigma=$(sigma) --repeats=1"
      export CPUOUT="$(outputDir)/$(id).cpu.mgz"
      export GPUOUT="$(outputDir)/$(id).gpu.mgz"

      ${TM_BIN_DIR}/convgaussian_test $MAINCMDS \--output=$CPUOUT
      ${TM_BIN_DIR}/convgaussian_test_cuda $MAINCMDS \--output=$GPUOUT

      $(projectDir)/tools/mridiff.pl \--results=$(cmdResultFn) \
                                     \--gold=$CPUOUT \
                                     \--test=$GPUOUT \
                                     \--threshold=$(tol)

      testFinish -c $(cmdResultFn) -r $(resultFn) -t $(runtimeFn)
]]