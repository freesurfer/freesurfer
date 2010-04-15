-- A common file for the MRIconvolveGaussian tests

module( ..., package.seeall )


-- List of gaussian widths to use
sigmaVals = { 0.5, 2, 6.383, 12.766, 4, 15, 21, 29.2 }


-- Script to run
script = [[
      export MAINCMDS=" --input=$(projectDir)/inputs/$(input) --sigma=$(sigma) --repeats=1"
      export CPUOUT="$(outputDir)/$(id).cpu.mgz"
      export GPUOUT="$(outputDir)/$(id).gpu.mgz"

      convgaussian_test $MAINCMDS \--output=$CPUOUT
      convgaussian_test_cuda $MAINCMDS \--output=$GPUOUT

      $(projectDir)/tools/mridiff.pl \--results=$(cmdResultFn) \
                                     \--gold=$CPUOUT \
                                     \--test=$GPUOUT \
                                     \--threshold=$(tol)

      testFinish -c $(cmdResultFn) -r $(resultFn) -t $(runtimeFn)
]]