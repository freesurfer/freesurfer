
 1. kvlAtlasMeshAlphaDrawer.h
 2. kvlAtlasMeshVisitCounter.h
 [3. kvlAtlasMeshCollectionPositionCostCalculator.h
     -> only used in kvlAtlasMeshCollectionPositionCostCalculator.cxx, which doesn't seem to be used anywhere
        Should be moved to old; current version is implemented in kvlAtlasMeshCollectionPositionCostCalculator2.cxx
        which doesn't use the rasterizor]
 [4. kvlAtlasMeshDeformationLevenbergMarquardt2.h
      -> only used in kvlAtlasMeshDeformationLevenbergMarquardt2.cxx, which itself doesn't seem to be used anywhere]
 [5. kvlAtlasMeshDeformationLevenbergMarquardt.h
      -> is used in:
           - kvlAtlasMeshDeformationLevenbergMarquardt.cxx -> is being  compiled
           - kvlAtlasMeshDeformationLevenbergMarquardtOptimizer.cxx 
               -> is being compiled
               -> is a kvlAtlasMeshDeformationOptimizer
               -> has a kvlAtlasMeshDeformationLevenbergMarquardt
               -> used in kvlAtlasParameterEstimator.h (in active use for building atlas meshes), but
                  only as an option to select a deformationOptimizer, and could certainly be removed
               -> used in Matlab/kvlDeformOneStep.h as one of the optimizer options, but could certainly be removed
                  (also exposed as Matlab/kvlGetLevenbergMarquardtOptimizer.h, which could also be removed) 
                  and therefore being compiled in Matlab/kvlMatlabRunnerArray.cxx, 
                  and being used in Matlab/preprocessHippoSubfields.m,
                                    Matlab/preprocessHippoSubfields_messedUpWithNewOptimizer.m
                                    Matlab/processHippoSubfields.m
                                    Matlab/segment.m,
                  but this is old stuff and could certainly be removed                
           - kvlAtlasMeshSegmenter.cxx 
                 -> is being compiled 
                 -> included in Matlab/kvlDeformMesh.h (which subclasses it), which is compiled in Matlab/kvlMatlabRunnerArray.cxx, but it's never used and could therefore certainly be removed
                 -> also included in Matlab/kvlSetOptimizerProperties.h,
                    but that seems an error
                 -> only other usage is in kvlAtlasMeshSegmentationDriver, which is being compiled but 
                    only used in Executables/kvlSegmentWithoutGUI.cxx, which is built in Executables,
                    and used in Matlab/kvlPreprocessHippocampalSubfields.sh, but this in turn is the
                    old script that I wrote for hippocampal subfield in FS 5.1! 
      -> double-checked: Eugenio uses LevenbergMarquardt as one of two options in his scripts ("optimizerType") -- the other being Conjugate Gradient. Can certainly be removed (actually, have he's using his current scripts Matlab will complain for not finding the LevenbergMarquardt wrapper if that's selected, so not really urgent)              
      -> it's massively complicated, uses GMM for sparse matrices -> remove]
      
 6. kvlAtlasMeshLabelStatisticsCollector.h
      -> is used in mesh estimation 
 7. kvlAtlasMeshMinLogLikelihoodCalculator.h
      -> is used in AtlasParameterEstimator (function GetMinLogLikelihood)
 8. kvlAtlasMeshMultiAlphaDrawer.h
      -> is used in Matlab/kvlRasterizeAtlasMesh.h
      -> is used in kvlAtlasMeshSmoother.cxx  -> used in Matlab code etc
      -> is used in kvlAtlasMeshToIntensityImagePartialVolumeGradientCalculator.cxx
           - is compiled
           - is used in kvlAtlasMeshSegmenter.h,  which could be removed as already stated above
           - so should be removed
      -> is used in kvlEMSegmenter.cxx
           - is compiled
           - is used in kvlAtlasMeshSegmenter.h, which could be removed as already stated above
           - is used in Executables/kvlEMSegment.cxx, which is being built but could be removed
           - is used in Executables/kvlSamplePositionsFromMeshCollection.cxx, which is never built and could be removed
           - is used in Matlab/kvlBiasFieldCorrectImageBuffer.h -> never used so could be removed
           - is used in Matlab/kvlBiasFieldCorrectImage.h -> is being built into Matlab/kvlMatlabRunnerArray.cxx,
             but never used so should be removed

 9. kvlAtlasMeshPositionGradientCalculator.h
      -> used for building atlas meshes etc
      -> kvlAtlasMeshVertexProcessor2.cxx
          - in function ::CalculateCostAndGradient(), it's used to calculate cost/gradient of prior and data
            part separately (using a lot of hacking). Not sure why that's needed, but can probably be replaced
            by kvlAtlasMeshToIntensityImageGradientCalculator when using the SetSegmentedImage() option (see below)
      -> kvlAtlasMeshVertexProcessor.cxx -> is not actually built, so should be removed (it's complete content,
         as well as that of the header file, is  completely "#if 0"-out
           - is used in kvlAtlasMeshCollectionPositionCostCalculator.h/cxx, but both these files' 
             content is completely "#if 0"-out so that should also be removed
      -> kvlAtlasParameterEstimator.cxx -> but there it could probably replaced by kvlAtlasMeshToIntensityImageGradientCalculator, as detailed below
      
[10. kvlAtlasMeshSummaryDrawer.h
      -> creates a gray-scale weighted-prior-probability map from an atlas (cf. TMI 2009 paper figures etc)
      -> used in kvlEstimateAtlasParametersWithoutGUI.cxx, but that isn't used anymore (should be moved to old)
      -> used in GUI/kvlAtlasMeshViewingConsole.cxx -> is being used as one display option in GUI/kvlViewMeshCollectionWithGUI.cxx
      -> used in GUI/kvlAtlasParameterEstimationConsole.cxx, which is used for executable kvlEstimateAtlasParameters (which actually has a GUI)
      ]
11. kvlAtlasMeshToIntensityImageGradientCalculator.h
      -> main deformation-for-segmentation gradient calculator  
      -> used in kvlAtlasMeshDeformationConjugateGradientOptimizer.h, which is the core optimizer
      -> used in kvlAtlasMeshSegmenter.cxx but that can be removed as already deduced earlier
      -> used in kvlAtlasMeshAveragingConjugateGradientOptimizer.h, which is used in Executables/kvlAverageMeshes.cxx
      -> used in kvlAtlasMeshHamiltonianPositionSampler.cxx
           - not being built so should be removed
           - used in Executables/kvlSamplePositionsFromMeshCollection.cxx, which is never built and should be removed as already deduced earlier

[12. kvlAtlasMeshToIntensityImagePartialVolumeGradientCalculator.h 
      -> could be removed as already stated earlier ]

[13. kvlAtlasMeshToProbabilityImageGradientCalculator.h
      -> only being used in kvlAtlasMeshSegmenter, which can be removed as already deduced earlier ]


-> in terms of gradients, there is 
      kvlAtlasMeshPositionGradientCalculator (for atlas building, but can be replaced by the one for segmentation purposes, as detailed below) 
    and 
      kvlAtlasMeshToIntensityImageGradientCalculator (for segmentation purposes)
    so these need to be reimplemented. Presumably they have a whole lot in common - presumably 
    calls for subclassing
    
-> I'm starting with kvlAtlasMeshToIntensityImageGradientCalculator
     - it has an option SetProbabilityImage, which is propagated in kvlAtlasMeshDeformationConjugateGradientOptimizer and its parent kvlAtlasMeshDeformationOptimizer.h (GetCostAndGradient()). Also in various LevenbergMarquardt stuff and kvlAtlasMeshSegmenter.cxx, but we're going to remove those anyway. 
     SO: this can be removed
     - it also has an option SetSegmentedImage(), but this is actually propagated through the kvlAtlasMeshDeformationOptimizer class into kvlAtlasParameterEstimator (needed for mesh building), where it's actually used. This is
     of course weird, because the functionality presumably is then identical kvlAtlasMeshPositionGradientCalculator, which is also only used for mesh building ??
     Need to check this -> OK I went to check and essentially within the function kvlAtlasParameterEstimator::EstimatePosition(), one of two things happen:
       (1) if Eugenio's flag m_GD = false, then the kvlAtlasMeshDeformationOptimizer framework is called (using the SetSegmentedImage() option). The function CalculateCurrentPositionGradient() is still called for cases where the mesh isn't moving
       (2) otherwise a hand-crafted silly gradient descent optimization routine is used, where the gradients come from 
       the CalculateCurrentPositionGradient() function I mentioned earlier
     Now CalculateCurrentPositionGradient() is what calls the odd AtlasMeshPositionGradientCalculator (although presumably it has the same functionality as the kvlAtlasMeshToIntensityImageGradientCalculator when using the SetSegmentedImage() option) -> YES: the input image is exactly the same (m_LabelImages[ labelImageNumber ]) for both gradient calculator implementations!
     The interface to set the m_GD flag is through:
       - SetModeLM()
       - SetModeCJ()
       - SetModeGD()
       - SetGD()
     By default m_GD=true, meaning the silly gradient descent optimizer is used. Now, kvlAtlasParameterEstimator is
     only used in a few places:
       - kvlAtlasMeshCollectionFastReferencePositionCost.h/kvlAtlasMeshCollectionReferencePositionCost.h
           -> both are being compiled
           -> both are included in kvlAtlasMeshBuilder.cxx:
               * the non-fast version is used in TryToCollapse() [in turn used in AnalyzeEdge() and LoadBalancedAnalyzeEdge()] and OptimizeReferencePosition() [in turn used in TryToRetain()]
               * the fast version is used in TryToCollapseFast() [in turn used in AnalyzeEdgeFast() and LoadBalancedAnalyzeEdgeFast()] and OptimizeReferencePositionFast() [in turn used in TryToRetainFast()]
               * now AnalyzeEdge() is used in ThreaderCallback, LoadBalancedAnalyzeEdge() is never used, and TryToRetain() is used in AnalyzeEdge() and LoadBalancedAnalyzeEdge(). So effectively it's only ever used in ThreaderCallback, which is never used.
               * So all of the non-fast stuff should be removed (both from Builder as from files)
               * There is also a FastThreaderCallback, which is never used either and should be removed
               * The only thing used is LoadBalancedThreaderCallback, which calls  LoadBalancedAnalyzeEdgeFast()
               * Anyway for now this all just means that we have to look at the m_GD flag setting in kvlAtlasMeshCollectionFastReferencePositionCost. There we have:
               
                 switch ( OPTIMIZER_SECTION_2 )
                   {
                   case CONJUGATE_GRADIENT: m_Estimator->SetModeCJ(); break;
                   case GRADIENT_DESCENT: m_Estimator->SetModeGD(); break;
                   case LEVENBERG_MARQUARDT: m_Estimator->SetModeLM(); break;
                   default: break;
                   }

       - kvlMultiResolutionAtlasMesher.h
            -> and there we find:

                switch ( OPTIMIZER_SECTION_1 )
                  {
                  case CONJUGATE_GRADIENT: m_Estimator->SetModeCJ(); break;
                  case GRADIENT_DESCENT: m_Estimator->SetModeGD(); break;
                  case LEVENBERG_MARQUARDT: m_Estimator->SetModeLM(); break;
                  default: break;
                  }
       - OK, now instead of just hardcoding the selected options directly, Eugenio has created a 
         separate file kvlOptimizerChoice.h, which reads:
         
            #define OPTIMIZER_SECTION_1 CONJUGATE_GRADIENT
            #define OPTIMIZER_SECTION_2 GRADIENT_DESCENT
            
        - So now we have it: for estimating deformation of entire meshes in the atlas building process,
          a (proper) conjugate gradient optimizer is used. For subsequent optimization of mini-meshes,
          a very crappy gradient descent is used
       
-> OK, so I'm just never going to implement kvlAtlasMeshPositionGradientCalculator because it's functionality is the same as that of kvlAtlasMeshToIntensityImageGradientCalculator. However, the latter one I'll want to have two subclasses: one for label images, and one for intensity images


--

TODO:

* Fold mixture-models-for-each-class into the likelihood filter; for whole-brain segmentation this will reduce the number of classes from 17+4=21 to 6+4=10, more or less doubling the speed of the gradient-based deformation optimization (4 extra because the baricentric coordinates are always computed).
 
* Disperse functionality in kvlAtlasMeshToIntensityImageGradientCalculatorGPU over a base class 
  (containing everything except the intensity-based stuff) and a derived class, as this will allow us
  to quickly and very cleanly write a class kvlAtlasMeshToSegmentationImageGradientCalculatorGPU for
  registering to a label image (needed in the atlas building code)

* In kvlAtlasMeshAveragingConjugateGradientOptimizerGPU, make the member variable m_Calculator point to
  the base class, and implement the gradient-to-segmentation-image stuff (simply plug the new kvlAtlasMeshToSegmentationImageGradientCalculatorGPU in)
  
* The Conjugate Gradient algorithm implemented in kvlAtlasMeshAveragingConjugateGradientOptimizer(GPU) is
  very sensitive to tiny variations in the gradient computations, exemplified by different optimization paths
  taken when running multi-threaded gradient computations. Have another look a the algorithm to make sure
  it all makes sense.
  
* While removing the non-used Levenberg-Marquardt stuff, make the gradient descent stuff (used only for 
  building atlases) a proper subclass of the optimzer framework.
  


      
     


  