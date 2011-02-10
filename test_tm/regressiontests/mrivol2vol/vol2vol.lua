-- -*- lua -*-


-- A common file for the MRImean tests

module( ..., package.seeall )


require( "tools.freesurfer" )

-- Make access to the FreeSurfer module quicker
local FS = tools.freesurfer




-- Common arrays

interpModes = { "trilin", "nearest", }
transformFiles = { "translation1.fsl", "translation2.fsl",
		   "transform1.fsl" }



-- Function to return the list of files to use
function GetFiles( pattern )
   return FS.GetFiles( FS.MRIdir(), pattern )
end



-- Function to turn test parameters into a test name
local function testName( input, transform, precision, mode )
   local inputItems = FS.split( input, "-" )
   local transformItems = FS.split( transform, "%." )

   local myName = inputItems[1].."-"..inputItems[2]
   myName = myName.."-"..transformItems[1]
   myName = myName.."-"..precision.."Out"
   myName = myName.."-"..mode

   return myName
end


-- Function to generate a list of tests
function testGen( inputs, transforms, precisions, interps, threshold )
   local testTable = {}

   -- inputs is a table of input file names
   -- transform is a table of transform matrix names
   -- precisions is a table of output precisions
   -- interps is a table of interpolation modes

   for i,input in ipairs(inputs) do
      for p,precision in ipairs(precisions) do
	 for t,transform in ipairs(transforms) do
	    for im,interp in ipairs(interps) do
	       local tName = testName( input, transform, precision, interp )
	       table.insert( testTable, { id=tName,
					  input=FS.MRIdir()..input,
					  transform=transform,
					  outputPrecision=precision,
					  interp=interp,
					  tol=threshold } )
	    end
	 end
      end
   end

   return testTable
end


-- The script to run
script = [[
      $(submit JOBNAME="$(id)", TIME="01:00:00" )

      export MAINCMDS=" --mov $(input) --targ $(input) --fsl $(testDir)/$(transform) --interp $(interp) --precision $(outputPrecision)"
      export CPUOUT="$(outputDir)/$(id).cpu.mgz"
      export GPUOUT="$(outputDir)/$(id).gpu.mgz"

      ${TM_BIN_DIR}/mri_vol2vol $MAINCMDS \--o $CPUOUT
      ${TM_BIN_DIR}/mri_vol2vol_cuda $MAINCMDS \--o $GPUOUT

      $(projectDir)/tools/mridiff.pl \--results=$(cmdResultFn) \
                                     \--gold=$CPUOUT \
                                     \--test=$GPUOUT \
                                     \--threshold=$(tol)

      testFinish -c $(cmdResultFn) -r $(resultFn) -t $(runtimeFn)
]]