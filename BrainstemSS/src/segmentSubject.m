% Segments the brainstem from the original <MPRAGE
%
% This function is heavily based on the prototypes I wrote for the hippocampal subfields
%
% segmentSubject(subjectName,subjectDir,resolution,atlasMeshFileName,atlasDumpFileName,compressionLUTfileName,K,side,FSdir)
%
% - subjectName: FreeSurfer subject name
% - subjectDir: FreeSurfer subject directory
% - resolution: voxel size at which we want to work (in mm).
% - atlasMeshFileName: the atlas to segment the data
% - atlasDumpFileName: corresponding imageDump.mgz (name *must* be imageDump.mgz)
% - compressionLUTfileName: corresponding compressionLUT.txt
% - K: stiffness of the mesh in the segmentation.
% - optimizerType: 'FixedStepGradientDescent','GradientDescent','ConjugateGradient','L-BFGS'
% - suffix: for output directory, e.g. 'T1based_GGAWLnoSimil','v10',...
% - FSdir: directory with Freesurfer binaries
% - additionalFile: use this file instad of the FS T1 in the analysis. Leave empty if you want to use the T1


function segmentSubject(subjectName,subjectDir,resolution,atlasMeshFileName,atlasDumpFileName,compressionLUTfileName,K,optimizerType,suffix,FSdir,additionalFile)

% Eugenio November 2017: added option for smoother resampling
DEBUG=0;
FAST=0; % set it two one to optimize just a bit (go through code fast) or to 2 to make it super-fast
SMOOTH_LABEL_RESAMPLE=0;
WRITE_POSTERIORS=0;
aux=getenv('WRITE_POSTERIORS');
if ~isempty(aux)
    if str2double(aux)>0
        WRITE_POSTERIORS=1;
    end
end


% sanity check
if nargin<7
    error('Not enough input arguments');
elseif optimizerType(1)~='F' && optimizerType(1)~='G' && optimizerType(1)~='C' && optimizerType(1)~='L'
    error('Optimizer type must be ''FixedStepGradientDescent'',''GradientDescent'',''ConjugateGradient'',''L-BFGS''');
elseif exist([subjectDir '/' subjectName],'dir')==0
    error('Subject directory does not exist');
elseif ~isdeployed && (~isnumeric(resolution))
    error('Resolution must be numeric');
elseif exist(atlasMeshFileName,'file')==0
    error('Provided atlas mesh file does not exist');
elseif exist(atlasDumpFileName,'file')==0
    error('Provided imageDump.mgz does not exist');
elseif exist(compressionLUTfileName,'file')==0
    error('Provided LUT does not exist');
elseif ~isdeployed && (~isnumeric(K))
    error('K must be numeric');
end


% Constants
BRAINSTEM=16;
DE_label_left=28;
DE_label_right=60;

% In case we compiled it...
if isdeployed
    K=str2double(K);
    resolution=str2double(resolution);
else
    addpath([pwd() '/functions']);
    addpath('/usr/local/freesurfer/stable6_0_0/matlab')
    if ismac
        addpath('/autofs/space/panamint_005/users/iglesias/software/freesurfer.GEMS2.MAC/bin')
    elseif isunix
        addpath('/autofs/space/panamint_005/users/iglesias/software/freesurfer.GEMS2/bin')
    else
        error('Neither Linux nor Mac');
    end
end
time_start=clock;

% Clean up KVL memory space
kvlClear;




% Temporary directory: here we have a secret flag: if we are at the Martinos
% Center and we are using the cluster, we want to set USE_SCRATCH to 1 in order
% to avoid massive data flow between the cluster and your machine (assming your
% data is local).
tempdir=[subjectDir '/' subjectName '/tmp/BS_basedOn_'];
aux=getenv('USE_SCRATCH');
if ~isempty(aux)
    if str2double(aux)>0
        if exist('/scratch/','dir')>0
            tempdir=['/scratch/' subjectName '/tmp/BS_basedOn_'];
        end
    end
end
if exist('additionalFile','var')>0
    tempdir=[tempdir 'AddVol/'] ;
    additionalFile=getFullPath(additionalFile); % Eugenio November 2017: before we cd
else
    tempdir=[tempdir 'T1/'] ;
end

if exist(tempdir,'dir')==0
    mkdir(tempdir);
end

cd(tempdir);




% Next: register image dump to automated segmentation
disp('Registering imageDump.mgz to BS mask from ASEG')

% Apply "our template"-to-"FS atlas" transform (pre-computed).
system(['cp ' atlasDumpFileName ' ./imageDump.mgz']);

% Target is masked aseg
targetRegFileName=[tempdir '/BS-DE-binaryMask.mgz'];
targetRegFileNameCropped=[tempdir '/BS-DE-binaryMask_autoCropped.mgz'];
system(['cp ' subjectDir '/' subjectName '/mri/aseg.mgz ./']);
ASEG=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
TARGETREG=ASEG;
TARGETREG.vol=255*double(ASEG.vol==BRAINSTEM);
myMRIwrite(TARGETREG,targetRegFileName,'float',tempdir);

highres=0; if mean(ASEG.volres)<0.99, highres=1; end
if highres==1,
    system([FSdir '/mri_convert ' targetRegFileName ' aux.mgz -odt float -vs 1 1 1 -rt nearest >/dev/null']);
    system(['mv aux.mgz ' targetRegFileName ' >/dev/null']);
end

% cmd=[FSdir '/kvlAutoCrop ' targetRegFileName ' 6'];
% system([cmd ' >/dev/null']);

aux=myMRIread(targetRegFileName,0,tempdir);
[aux.vol,cropping]=cropLabelVol(aux.vol,6);
shift=aux.vox2ras0(1:3,1:3)*[cropping(2)-1; cropping(1)-1; cropping(3)-1];
aux.vox2ras0(1:3,4)=aux.vox2ras0(1:3,4)+shift;
aux.vox2ras1(1:3,4)=aux.vox2ras1(1:3,4)+shift;
aux.vox2ras(1:3,4)=aux.vox2ras(1:3,4)+shift;
aux.tkrvox2ras=[];
myMRIwrite(aux,targetRegFileNameCropped,'float',tempdir);


% Initial affine alignment based just on brainstem
% cmd=[FSdir '/kvlRegister imageDump.mgz ' targetRegFileNameCropped ' 3 2'];
% system(cmd);
% % system([cmd ' >/dev/null']);
% system('mv imageDump_coregistered.mgz imageDump.mgz' );

cmd=[FSdir '/mri_robust_register --mov imageDump.mgz  --dst ' targetRegFileNameCropped ...
    ' -lta trash.lta --mapmovhdr imageDump_coregistered.mgz  --sat 50'];
system(cmd);
% system([cmd ' >/dev/null']);
system('mv imageDump_coregistered.mgz imageDump.mgz' );

cmd=[FSdir '/mri_robust_register --mov imageDump.mgz  --dst ' targetRegFileNameCropped ...
    ' -lta trash.lta --mapmovhdr imageDump_coregistered.mgz --affine --sat 50'];
system(cmd);
% system([cmd ' >/dev/null']);
system('mv imageDump_coregistered.mgz imageDump.mgz' );



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we pretty much copy-paste from preprocessHippoSubfields %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

boundingFileName = [tempdir 'imageDump.mgz']; % Bounding box
meshCollectionFileName = atlasMeshFileName; % The tetrahedral atlas mesh
compressionLookupTableFileName =compressionLUTfileName; % Look-up table belonging to the atlas

%
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );


% Read in aseg and also transform
aux=ASEG;
aux.vol(:)=1;

aux.vol(ASEG.vol==BRAINSTEM | ASEG.vol==7 | ASEG.vol==8 | ASEG.vol==15 | ASEG.vol==28 ...
    | ASEG.vol==46 | ASEG.vol==47 | ASEG.vol==60) = 255;
myMRIwrite(aux,'cheating.mgz','float',tempdir);
[ synIm, transform ] = kvlReadCroppedImage( 'cheating.mgz', boundingFileName );
synImBuffer = kvlGetImageBuffer( synIm );
synSize = size( synImBuffer );
if ~isdeployed && DEBUG>0
    figure
    showImage( synIm )
    title('Synthetic (cheating) Image to segment')
end


% read in collection, set K and apply transform
meshCollection = kvlReadMeshCollection( meshCollectionFileName );
kvlTransformMeshCollection( meshCollection, transform );
kvlSetKOfMeshCollection( meshCollection, K );

% Retrieve the reference mesh, i.e., the mesh representing the average shape.
mesh = kvlGetMesh( meshCollection, -1 );
originalNodePositions = kvlGetMeshNodePositions( mesh );
originalAlphas = kvlGetAlphasInMeshNodes( mesh );


% In this first step, we're cheating in that we're going to segment the
% "fake" synthetic image obtained a binary mask of the cerebellum, 4rd ventricle,
% brainstem and dicencephalon. We group all these because there were some differences
% in the definitions of the boundaries between these strucutures from ASEG to UCSF
% Then, we assign those well-known means (and little variance) to the classes in the mesh
% We assign to each class its corresponding mean
% We also build compressed classes (which group together the hippocampus)
sameGaussianParameters=[];
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'MAC_Medulla','MAC_Pons','MAC_Midbrain','Left-VentralDC','4th-Ventricle','Left-Cerebellum-White-Matter','Left-Cerebellum-Cortex','MAC_Sup_Cerebellum_Ped','Medulla','Pons','SCP','Midbrain'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Caudate','Left-Accumbens-area','Left-Pallidum','3rd-Ventricle','Left-Putamen','Left-Thalamus-Proper','Left-Amygdala','Left-Lateral-Ventricle','Left-choroid-plexus','Left-Hippocampus','Left-Cerebral-White-Matter','Left-Cerebral-Cortex','Background-tissue','Background-CSF','Background'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
cheatingMeans=[255;1];
cheatingVariances=[1;1];



% Compute the "reduced" alphas - those referring to the "super"-structures
[ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
if ( max( abs( sum( reducedAlphas, 2 ) - 1 ) ) > 1e-5 ) % Make sure these vectors really sum to 1
    error( 'The vector of prior probabilities in the mesh nodes must always sum to one over all classes' )
end

% Set the reduced alphas to be the alphas of the mesh
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )



% Eugenio July 2017

% priors = kvlRasterizeAtlasMesh( mesh, synSize );
% MASK=imerode(sum(double(priors)/65535,4)>0.99,createSphericalStrel(5)); % borders are often a bit weird...

for l = 1 : size(reducedAlphas,2)    
    if l==1
       % This is a bit annoying, but the call to kvlRasterize with a single
       % label fills in the voxels outside the cuboid with l=1 (whereas the
       % call with multiple labels does not)
        sillyAlphas=zeros([size(reducedAlphas,1),2],'single');
        sillyAlphas(:,1)=reducedAlphas(:,1);
        sillyAlphas(:,2)=1-sillyAlphas(:,1);
        kvlSetAlphasInMeshNodes( mesh, sillyAlphas )
        prior = kvlRasterizeAtlasMesh( mesh, synSize);
        kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
        sumpriors=prior(:,:,:,1);
    else
        prior = kvlRasterizeAtlasMesh( mesh, synSize, l-1 );
        sumpriors=sumpriors+prior;        
    end
end
MASK=imerode(sum(single( sumpriors / 65535 ),4)>0.99,createSphericalStrel(5));


cheatingImageBuffer=synImBuffer;
cheatingImageBuffer(~MASK)=0;
cheatingImage = kvlCreateImage( cheatingImageBuffer );


if ~isdeployed && DEBUG>0
    priors = kvlRasterizeAtlasMesh( mesh, synSize );
    figure
    for cheatingLabel = 1 : size( reducedAlphas, 2 )
        subplot( 3, 4, cheatingLabel )
        showImage( priors( :, :, :, cheatingLabel ) )
    end
    title('Priors for segmentation of fake intensity image')
end

if ~isdeployed && DEBUG>0
    figure
    showImage( cheatingImage, [], [ 0 110 ] )  % Using dynamic range for display that shows what's going on
    title('Fake intensity image to initialize atlas deformation')
end


% It's good to smooth the mesh, otherwise we get weird compressions of the
% mesh along the boundaries...
meshSmoothingSigma = 3.0; % This is the value Koen had...
fprintf( 'Smoothing mesh with kernel size %f ...', meshSmoothingSigma )
kvlSmoothMeshCollection( meshCollection, meshSmoothingSigma )
fprintf( 'done\n' )
if ~isdeployed && DEBUG>0
    figure
    priors = kvlRasterizeAtlasMesh( mesh, synSize );
    for cheatingLabel = 1 : size( reducedAlphas, 2 )
        subplot( 3, 4, cheatingLabel )
        showImage( priors( :, :, :, cheatingLabel ) )
    end
    title('Smoothed priors')
end

% Eugenio November 2017: GEMS2, more iterations   (note that it uses variances instead of precisions)
% Now the optimization per-se

verbose=0;
maximalDeformationStopCriterion=1e-10;
relativeChangeInCostStopCriterion = 1e-10;
lineSearchMaximalDeformationIntervalStopCriterion=1e-10;
maximumNumberOfDeformationIterations=1000;
BFGSMaximumMemoryLength=12;
maxpuin=300;
if FAST>0,maxpuin=60; end
if FAST>1,maxpuin=5; end

cheatingCalculator = kvlGetCostAndGradientCalculator('AtlasMeshToIntensityImage',...
    cheatingImage, 'Sliding',transform,cheatingMeans,cheatingVariances,ones(size(cheatingMeans)),ones(size(cheatingMeans)));
cheatingOptimizer = kvlGetOptimizer( optimizerType, mesh, cheatingCalculator, ...
    'Verbose', verbose, ...
    'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
    'LineSearchMaximalDeformationIntervalStopCriterion', ...
    lineSearchMaximalDeformationIntervalStopCriterion, ...
    'MaximumNumberOfIterations', maximumNumberOfDeformationIterations, ...
    'BFGS-MaximumMemoryLength', BFGSMaximumMemoryLength );
    

historyOfMinLogLikelihoodTimesPrior = [ 1/eps ];

time_ref_cheat_optimization=clock;
for positionUpdatingIterationNumber = 1 : maxpuin
    disp(['Iteration ' num2str(positionUpdatingIterationNumber)]);
    % Calculate a good step. The first one is very slow because of various set-up issues
    % Eugenio May2018
    maximalDeformation=0;
    try
        tic
        % Eugenio November 2017: GEMS2
        [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( cheatingOptimizer );
        
        elapsedTime = toc;
        disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
        minLogLikelihoodTimesPrior
    end
    if isnan(minLogLikelihoodTimesPrior)
        error('lhood is nan');
    end
    
    historyOfMinLogLikelihoodTimesPrior = [ historyOfMinLogLikelihoodTimesPrior; minLogLikelihoodTimesPrior ];
    
    % Test if we need to stop
    if ( ( maximalDeformation <= maximalDeformationStopCriterion ) | ...
            abs(( ( historyOfMinLogLikelihoodTimesPrior( end-1 ) - historyOfMinLogLikelihoodTimesPrior( end ) ) ...
            / historyOfMinLogLikelihoodTimesPrior( end ))) < relativeChangeInCostStopCriterion   )
        break;
    end
    
    % Show what we have
    if ~isdeployed && DEBUG>0
        subplot( 2, 2, 3 )
        showImage( kvlRasterizeAtlasMesh( mesh, synSize, 1 ) );
        title('Current atlas deformation')
        subplot( 2, 2, 4 )
        plot( historyOfMinLogLikelihoodTimesPrior( 2 : end ) )
        title('History of Log-lhood + Log-prior')
        drawnow
    end
end
% Eugenio November 2017: GEMS2
kvlClear( cheatingOptimizer )
kvlClear( cheatingCalculator )

disp(['Fitting mesh to synthetic image from ASEG took ' num2str(etime(clock,time_ref_cheat_optimization)) ' seconds']);

if positionUpdatingIterationNumber==1
    error('Fitting mesh to synthetic image resulted in no deformation')
end



% OK, we're done. let's modify the mesh atlas in such a way that our computed mesh node positions are
% assigned to what was originally the mesh warp corresponding to the first training subject.
kvlSetAlphasInMeshNodes( mesh, originalAlphas )
updatedNodePositions = kvlGetMeshNodePositions( mesh );
kvlSetMeshCollectionPositions( meshCollection, ... %
    originalNodePositions, ... % reference position (average "shape")
    updatedNodePositions );


% Compare the average shape we started with, with the shape we have computed now in a little movie
% Compare the average shape we started with, with the shape we have computed now in a little movie
if ~isdeployed && DEBUG>0
    figure
    originalPositionColorCodedPriors = ...
        kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( kvlGetMesh( meshCollection, -1 ), synSize ), colors );
    updatedPositionColorCodedPriors = ...
        kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( kvlGetMesh( meshCollection, 0 ), synSize ), colors );
    for i=1:20
        subplot(1,3,1)
        showImage( cheatingImage ), title('Cheating image')
        subplot(1,3,2)
        showImage( cheatingImage ), title('Cheating image')
        subplot(1,3,3)
        showImage( originalPositionColorCodedPriors ), title('Prior before deformation')
        pause(.5)
        subplot(1,3,1)
        showImage( originalPositionColorCodedPriors ), title('Prior before deformation')
        subplot(1,3,2)
        showImage( updatedPositionColorCodedPriors ), title('Prior after deformation')
        subplot(1,3,3)
        showImage( updatedPositionColorCodedPriors ), title('Prior after deformation')
        pause( .5 )
    end
end


% Write the resulting atlas mesh to file IN NATIVE ATLAS SPACE
transformMatrix = double(kvlGetTransformMatrix( transform ));
inverseTransform = kvlCreateTransform( inv( transformMatrix ) );
kvlTransformMeshCollection( meshCollection, inverseTransform );
kvlWriteMeshCollection( meshCollection, 'warpedOriginalMesh.txt' );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we pretty much copy-paste from processHippoSubfields %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clean up the Matlab work space
kvlClear % Clear all the wrapped C++ stuff
close all

% Provide the location of the image to be segmented, as well as the atlas that has been
% pre-registered affinely (i.e., 12 degrees of freedom) to the image.
meshCollectionFileName = 'warpedOriginalMesh.txt.gz'; % The tetrahedral atlas mesh


% Eugenio: we extract a block from norn.mgz and upsample it to the work
% resolution.
% We also mask out non-brain voxels pretty aggressively
imageFileName='T1resampled.mgz';
margin=15; % in mm
marginASEG=5;
if exist('additionalFile','var')>0
    system([FSdir '/mri_convert ' subjectDir '/' subjectName '/mri/norm.mgz ' ...
        tempdir '/normResampled.mgz -odt float -rl ' additionalFile ]);
    system([FSdir '/mri_binarize --min 0.5 --i ' tempdir  '/normResampled.mgz '...
        ' --o ' tempdir '/normMask.mgz']);
    system([FSdir '/mri_mask ' additionalFile ' ' tempdir '/normMask.mgz ' ...
        tempdir '/additionalFileMasked.mgz']);
    system([FSdir '/mri_convert ' subjectDir '/' subjectName '/mri/aseg.mgz ' ...
        tempdir '/asegAddSpace.mgz -odt float -rt nearest -rl ' additionalFile ]);
    A=myMRIread([tempdir '/additionalFileMasked.mgz'],0,tempdir);
    L=myMRIread([tempdir '/asegAddSpace.mgz'],0,tempdir);
else
    A=myMRIread([subjectDir '/' subjectName '/mri/norm.mgz'],0,tempdir);
    L=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
end

[~,cropping]=cropLabelVol(L.vol==BRAINSTEM | L.vol==DE_label_left | L.vol==DE_label_right,round(margin/mean(A.volres)));
Lcrop=applyCropping(L.vol,cropping);
Icrop=applyCropping(A.vol,cropping);
offsetVox=[cropping(2)-1;cropping(1)-1;cropping(3)-1;0];
RAScorner=L.vox2ras0*offsetVox;
vox2ras0New=[ L.vox2ras0(:,1:3)  L.vox2ras0(:,4)+RAScorner];
I=L;
I.vol=Icrop;
I.vox2ras0=vox2ras0New;
myMRIwrite(I,'tempT1.mgz','float',tempdir);
I.vol=Lcrop;
myMRIwrite(I,'tempASEG.mgz','float',tempdir);
cmd=[FSdir '/mri_convert tempT1.mgz '  imageFileName ' -odt float -rt  cubic -vs ' num2str(resolution) ' '  num2str(resolution)  ' '  num2str(resolution)];
system([cmd ' >/dev/null']);
delete tempT1.mgz
cmd=[FSdir '/mri_convert tempASEG.mgz tempASEGresampled.mgz -odt float -rt  nearest -vs ' num2str(resolution) ' '  num2str(resolution)  ' '  num2str(resolution)];
system([cmd ' >/dev/null']);
delete tempASEG.mgz
aux=myMRIread('tempASEGresampled.mgz',0,tempdir);
aux.vol=imdilate(aux.vol>0,createSphericalStrel(round(marginASEG/resolution)));
aux2=myMRIread( imageFileName,0,tempdir);
aux2.vol(aux.vol==0)=0;
myMRIwrite(aux2,imageFileName,'float',tempdir);
delete tempASEGresampled.mgz


% Read the image data from disk. At the same time, construct a 3-D affine transformation (i.e.,
% translation, rotation, scaling, and skewing) as well - this transformation will later be used
% to initially transform the location of the atlas mesh's nodes into the coordinate system of
% the image.
[ image, transform ] = kvlReadCroppedImage( imageFileName, boundingFileName );
imageBuffer = kvlGetImageBuffer( image );
imageSize = size( imageBuffer );
if ~isdeployed && DEBUG>0
    figure
    showImage( imageBuffer ); % Automatically displays middle slices in each direction
    title('Image buffer to segment')
end

% Read the atlas mesh from file, and apply the previously determined transform to the location
% of its nodes.
% If you don't provide an explicit value for K, the value used to construct the atlas will be used (which is
% theoretically the only really valid one...)
meshCollection = kvlReadMeshCollection( meshCollectionFileName );
kvlTransformMeshCollection( meshCollection, transform );
kvlSetKOfMeshCollection( meshCollection, K );


% Retrieve the correct mesh to use from the meshCollection. Above, we pre-deformed
% the average shape hippocampal subfield atlas to match the whole hippocampus segmentation generated by FreeSurfer,
% and wrote it out to the first mesh in the meshCollection we're using now. So make sure to retrieve that one
% and not something else
mesh = kvlGetMesh( meshCollection, 0 );  % Use -1 if you don't want to use the preprocessed one, but really the
% one with average shape


% For reasons that escape me right now, I seem to have somehow decided that areas not covered by the
% mesh have probability 1 of belonging to the class with index 0 (don't ask)
if ~isdeployed && DEBUG>0
    figure
    for labelNumber = 0 : 14  % C-style numbering: 0 corresponds to first element...
        prior = kvlRasterizeAtlasMesh( mesh, imageSize, labelNumber );
        subplot( 4, 4, labelNumber+1 )
        showImage( prior )
        title(['Priors for ' names(labelNumber+1,:)]);
    end
end

% Rather than showing rasterized priors defined by the atlas mesh one by one, as we did above,
% we can color-code them and show everything as one image
if ~isdeployed && DEBUG>0
    priors = kvlRasterizeAtlasMesh( mesh, imageSize ); % Without specifying a specific label, will rasterize all simultaneously
    colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
    figure
    showImage( colorCodedPriors );
    title('Color-coded priors')
    clear colorCodedPriors
end


% Eugenio July 2017: let's avoid rasterizing all priors simultaneously
% (which was faster, but used much more memory)
alphas = kvlGetAlphasInMeshNodes( mesh );
nlabels = size(alphas,2);
for l = 1 : nlabels
    
    if l==1
       % This is a bit annoying, but the call to kvlRasterize with a single
       % label fills in the voxels outside the cuboid with p=1 (whereas the
       % call with multiple labels does not)
        sillyAlphas=zeros([size(originalAlphas,1),2],'single');
        sillyAlphas(:,1)=originalAlphas(:,1);
        sillyAlphas(:,2)=1-sillyAlphas(:,1);
        kvlSetAlphasInMeshNodes( mesh, sillyAlphas )
        prior = kvlRasterizeAtlasMesh( mesh, imageSize);
        kvlSetAlphasInMeshNodes( mesh, originalAlphas )
        sumpriors=prior(:,:,:,1);

    else
       
        prior = kvlRasterizeAtlasMesh( mesh, imageSize, l-1 );
        sumpriors=sumpriors+prior;
        
    end
end

% We're not interested in image areas that fall outside our cuboid ROI where our atlas is defined. Therefore,
% generate a mask of what's inside the ROI. Also, by convention we're skipping all voxels whose intensities
% is exactly zero (which we do to remove image areas far away for FreeSurfer's ASEG hippo segmentation) -
% also include that in the mask
% mask = single( sum( priors, 4 ) / 65535 ) > 0.99;
mask = single( sumpriors / 65535 ) > 0.99;
mask = mask & ( imageBuffer > 0 );
if ~isdeployed && DEBUG>0
    figure
    subplot( 1, 2, 1 )
    showImage( mask )
    title('Mask given by mesh')
end


% Apply the mask to the image we're analyzing by setting the intensity of all voxels not belonging
% to the brain mask to zero. This will automatically discard those voxels in subsequent C++ routines, as
% voxels with intensity zero are simply skipped in the computations.
% Note that it is not sufficient to simply change the intensities in the Matlab matrix imageBuffer holding
% a copy of the image intensities - we also need to explicitly write any modifications to the ITK object!
imageBuffer( find( ~mask ) ) = 0;
kvlSetImageBuffer( image, imageBuffer );
imageBuffer = kvlGetImageBuffer( image );
if ~isdeployed && DEBUG>0
    subplot( 1, 2, 2 )
    showImage( imageBuffer )
    title('Masked image Buffer')
end



% Merge classes
%%%%%%%%%%%%%%%

sameGaussianParameters=[];
% Brainstem structures
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'MAC_Medulla','MAC_Pons','MAC_Midbrain','MAC_Sup_Cerebellum_Ped','Left-VentralDC','Medulla','Pons','SCP','Midbrain'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% CSF structures
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'3rd-Ventricle','Left-Lateral-Ventricle','Background-CSF','4th-Ventricle'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Gray matter structures
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Amygdala','Left-Cerebral-Cortex','Left-Hippocampus'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Caudate
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Caudate'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Accumbens area
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Accumbens-area'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Pallidum
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Pallidum'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Putamen
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Putamen'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Thalamus
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Thalamus-Proper'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Choroid plexus
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-choroid-plexus'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Cerebral white matter
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Cerebral-White-Matter'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% Background: misc tissue
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Background-tissue','Background'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% cerebellum white matter
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Cerebellum-White-Matter'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
% cerebellum cortex
sameGaussianParameters{end+1} = [];
for FreeSurferLabel = {'Left-Cerebellum-Cortex'}
    sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
end
numberOfClasses=length(sameGaussianParameters);



% Compute the "reduced" alphas - those referring to the "super"-structures
[ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
if ( max( abs( sum( reducedAlphas, 2 ) - 1 ) ) > 1e-5 ) % Make sure these vectors really sum to 1
    error( 'The vector of prior probabilities in the mesh nodes must always sum to one over all classes' )
end

% Set the reduced alphas to be the alphas of the mesh
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
if  ~isdeployed && DEBUG>0
    priors = kvlRasterizeAtlasMesh( mesh, imageSize );
    figure
    for reducedLabel = 1 : size( reducedAlphas, 2 )
        subplot( 4, 4, reducedLabel )
        showImage( priors( :, :, :, reducedLabel ) )
        title(['Prior for "reduced" class ' num2str(reducedLabel)]);
    end
    priorsBefore=priors;
    clear priors
end



% Compute hyperparameters for estimation of Gaussian parameters
disp('Computing hyperparameters for estimation of Gaussian parameters')

if exist('additionalFile','var')>0
    DATA=myMRIread([tempdir '/additionalFileMasked.mgz'],0,tempdir);
    ASEG=myMRIread([tempdir '/asegAddSpace.mgz'],0,tempdir);
else
    DATA=myMRIread([subjectDir '/' subjectName '/mri/norm.mgz'],0,tempdir);
    ASEG=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
end


nHyper=zeros(length(sameGaussianParameters),1);
meanHyper=zeros(length(sameGaussianParameters),1);
BGindex=0;
for g=1:length(sameGaussianParameters)
    labels=sameGaussianParameters{g};
    if any(labels==0)
        BGindex=g;
    end
    if any(labels==3 | labels==17 | labels==18)  % gray matter
        listMask=[3 42 17 53 18 54];
    elseif any(labels==2)  % white matter
        listMask=[2 41];
        
    elseif any(labels==178 | labels==34458 | labels==28)  % brainstem + diencephalon
        listMask=[16 28 60];
        
    elseif any(labels==4)  % CSF
        listMask=[4 43 14 15] ;
        
    elseif any(labels==11)  % caudate
        listMask=[11 50];
        
    elseif any(labels==26)  % accumbens
        listMask=[26 58];
        
    elseif any(labels==13)  % pallidum
        listMask=[13 52];
        
    elseif any(labels==12)  % putamen
        listMask=[12 51];
        
    elseif any(labels==10)  % thalamus
        listMask=[10 49];
        
    elseif any(labels==31)  % choroid
        listMask=[31 63];
        
    elseif any(labels==0)  % background
        listMask=[0];
        
    elseif any(labels==7)  % cerebellum WM
        listMask=[7 46];
        
    elseif any(labels==8)  % cerebellum CT
        listMask=[8 47];
    else
        listMask=[];
    end
    if length(listMask)>0
        MASK=zeros(size(DATA.vol));
        for l=1:length(listMask)
            MASK=MASK | ASEG.vol==listMask(l);
        end
        MASK=imerode(MASK,createSphericalStrel(1));
        data=DATA.vol(MASK & DATA.vol>0);
        meanHyper(g)=median(data);
        nHyper(g)=10+0.1*length(data)/prod(DATA.volres);
    end
end
% if any nan, replace by background
ind=find(isnan(meanHyper));
meanHyper(ind)=55;
nHyper(ind)=10;




%
% Multi-resolution scheme
%
% Specify here the size of the standard deviation of the Gaussian kernel used to smooth the priors/mesh. Use
% if you don't want to use multi-resolution

meshSmoothingSigmas = [ 2 1 0 ]';
imageSmoothingSigmas = [0 0 0 ]';
maxItNos=[7 5 3];  % each iteration has 40 deformation steps
if FAST>0
    meshSmoothingSigmas = [ 1 0 ]';
    imageSmoothingSigmas = [0 0 ]';
    maxItNos=[3 3];
end
if FAST>1
    meshSmoothingSigmas = [ 0 ]';
    imageSmoothingSigmas=[0];
    maxItNos=[2];
end


numberOfMultiResolutionLevels = length( meshSmoothingSigmas );

% Now the real work...
if  ~isdeployed && DEBUG>0
    multiResolutionFigure = figure;
    EMResultsFigure = figure;
    deformationMovieFigure = figure;
    costFigure = figure;
end
maskIndices = find( mask );


time_ref_optimization=clock;

imageBufferOrig=imageBuffer;

for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
    
    % Smooth the mesh using a Gaussian kernel.
    kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
    meshSmoothingSigma = meshSmoothingSigmas( multiResolutionLevel );
    fprintf( 'Smoothing mesh collection with kernel size %f ...', meshSmoothingSigma )
    kvlSmoothMeshCollection( meshCollection, meshSmoothingSigma )
    fprintf( 'done\n' )
    
    % Smooth the image using a Gaussian kernel
    imageSigma=imageSmoothingSigmas(multiResolutionLevel);
    if imageSigma>0
        imageBuffer=single(GaussFilt3dMask(imageBufferOrig,imageBufferOrig>0,imageSigma,resolution*ones(1,3)));
        kvlSetImageBuffer(image,imageBuffer);
    end
    
    
    % Show the smoothed atlas
    if  ~isdeployed && DEBUG>0
        figure( multiResolutionFigure )
        priors = kvlRasterizeAtlasMesh( mesh, imageSize );
        for reducedLabel = 1 : size( reducedAlphas, 2 )
            subplot( 3, 4, reducedLabel )
            showImage( priors( :, :, :, reducedLabel ) );
            title(['prior for reduced label ' num2str(reducedLabel)]);
        end
        clear priors
        subplot( 3, 4, 12 )
        showImage( kvlGetImageBuffer( image ) )
        title('image')
    end
    
    % Now with this smoothed atlas, we're ready for the real work. There are essentially two sets of parameters
    % to estimate in our generative model: (1) the mesh node locations (parameters of the prior), and (2) the
    % means and variances of the Gaussian intensity models (parameters of the
    % likelihood function, which is really a hugely simplistic model of the MR imaging process). Let's optimize
    % these two sets alternately until convergence. Optimizing the mesh node locations
    % is provided as a black box type of thing as it's implemented in C++ using complicated code - the other
    % set is much much better to experiment with in Matlab.
    
    
    maximumNumberOfIterations = maxItNos(multiResolutionLevel);  % Maximum number of iterations (includes one imaging model parameter estimation and
    % one deformation optimization; the latter always does 20 steps).
    

    historyOfCost = [ 1/eps ];
    
    % Compute a color coded version of the atlas prior in the atlas's current pose, i.e., *before*
    % we start deforming. We'll use this just for visualization purposes
    if  ~isdeployed && DEBUG>0
        oldColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, imageSize ) );
        figure( deformationMovieFigure )
        showImage( oldColorCodedPriors )
    end
    
    % Let's write this to file. In Unix/Linux systems, you can visualize progress over iterations
    % by doing animate -delay 10 *.png
    if  ~isdeployed && DEBUG>0
        fileName = [ 'colorCodedPrior_multiResolutionLevel' num2str( multiResolutionLevel ) '_iteration' sprintf( '%03d', 0 ) '.png' ];
        tmp = getframe( gcf );
        imwrite( tmp.cdata, fileName );
    end
    
    
    for iterationNumber = 1 : maximumNumberOfIterations
        disp(['Iteration ' num2str(iterationNumber) ' of ' num2str(maximumNumberOfIterations)]);
        %
        % Part I: estimate Gaussian mean and variances using EM
        %
        % See the paper
        %
        %     Automated Model-Based Bias Field Correction of MR Images of the Brain
        %     K. Van Leemput, F. Maes, D. Vandermeulen, P. Suetens
        %    IEEE Transactions on Medical Imaging, vol. 18, no. 10, pp. 885-896, October 1999
        %
        
        % Get the priors as dictated by the current mesh position, as well as the image intensities
        data = double( reshape( kvlGetImageBuffer( image ), [ prod( imageSize ) 1 ] ) ); % Easier to work with vector notation in the computations
        
        % Eugenio July 2017: again, avoid spike of memory use
        priors=zeros([length(maskIndices),numberOfClasses],'uint16');
        for l=1:numberOfClasses
            prior = kvlRasterizeAtlasMesh( mesh, imageSize, l-1 );
            priors(:,l)=prior(maskIndices);
        end
        
%         priors = kvlRasterizeAtlasMesh( mesh, imageSize );
%         priors = reshape( priors, [ prod( imageSize ) numberOfClasses ] ); % Easier to work with vector notation in the computations
%         priors = priors( maskIndices, : );
        
        
        data = data( maskIndices );
        
        % Start EM iterations. Initialize the parameters if this is the
        % first time ever you run this
        posteriors = double( priors ) / 65535;
        EPS=1e-2;
        if ( ( multiResolutionLevel == 1) & ( iterationNumber == 1 ) )
            
            for classNumber = 1 : numberOfClasses
                posterior = posteriors( :, classNumber );
                
                if sum(posterior)>EPS
                    
                    %   mu = data' * posterior / ( sum( posterior ) + eps );
                    %   variance = ( ( data - mu ).^2 )' * posterior / ( sum( posterior ) + eps );
                    
                    mu = (meanHyper(classNumber)*nHyper(classNumber) + data'*posterior) / ( nHyper(classNumber) + sum( posterior ) + EPS );
                    variance = (( ( data - mu ).^2 )' * posterior + nHyper(classNumber)*(mu-meanHyper(classNumber))^2 )/ ( sum( posterior ) + EPS );
                    
                    means( classNumber ) = mu;
                    variances( classNumber ) = variance+EPS;
                    
                else
                    means( classNumber ) = meanHyper(classNumber);
                    variances( classNumber ) = 100;
                    
                end
            end
            variances(variances==0)=100; % added by Eugenio, prevents nans...
            
            
        end % End test need for initialization
        stopCriterionEM = 1e-5;
        historyOfEMCost = [ 1/eps ];
        
        for EMIterationNumber = 1 : 100
            %
            % E-step: compute the posteriors based on the current parameters
            %
            minLogLikelihood = 0;
            for classNumber = 1 : numberOfClasses
                mu = means( classNumber );
                variance = variances( classNumber );
                prior = single( priors( :, classNumber ) ) / 65535;
                posteriors( :, classNumber ) = ( exp( -( data - mu ).^2 / 2 / variance ) .* prior ) ...
                    / sqrt( 2 * pi * variance);
                
                minLogLikelihood = minLogLikelihood+0.5*log(2*pi*variance)-0.5*log(nHyper(classNumber))...
                    +0.5*nHyper(classNumber)/variance*(mu-meanHyper(classNumber))^2; % contribution from prior
                
            end
            normalizer = sum( posteriors, 2 ) + eps;
            % Eugenio July 2017
            % posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfClasses ] );
            posteriors = bsxfun(@rdivide,posteriors, normalizer);
            
            minLogLikelihood =  minLogLikelihood - sum( log( normalizer ) ) % This is what we're optimizing with EM
            if isnan(minLogLikelihood)
                error('lhood is nan');
            end
            
            historyOfEMCost = [ historyOfEMCost; minLogLikelihood ];
            
            
            % Check for convergence
            relativeChangeCost = ( historyOfEMCost(end-1) - historyOfEMCost(end) ) / ...
                historyOfEMCost(end);
            if ( relativeChangeCost < stopCriterionEM )
                % Converged
                disp( 'EM converged!' )
                break;
            end
            
            
            % Show posteriors
            if  ~isdeployed && DEBUG>0
                figure( EMResultsFigure )
                posterior = zeros( imageSize );
                for classNumber = 1 : numberOfClasses
                    subplot( 3, 4, classNumber )
                    posterior( maskIndices ) = posteriors( :, classNumber );
                    showImage( posterior )
                    title(['Posterior for class ' num2str(classNumber)]);
                end
            end
            
            % Also show EM cost function
            if  ~isdeployed && DEBUG>0
                subplot( 3, 4, 12 )
                plot( historyOfEMCost( 2 : end ) )
                title( 'EM cost' )
            end
            
            % Refresh figures
            drawnow
            %pause
            
            
            %
            % M-step: derive parameters from the posteriors
            %
            
            % Update parameters of Gaussian mixture model
            EPS=1e-2;
            for classNumber = 1 : numberOfClasses
                posterior = posteriors( :, classNumber );
                
                if sum(posterior)>EPS
                    
                    %   mu = data' * posterior / ( sum( posterior ) + eps );
                    %   variance = ( ( data - mu ).^2 )' * posterior / ( sum( posterior ) + eps );
                    
                    mu = (meanHyper(classNumber)*nHyper(classNumber) + data'*posterior) / ( nHyper(classNumber) + sum( posterior ) + EPS );
                    variance = (( ( data - mu ).^2 )' * posterior + nHyper(classNumber)*(mu-meanHyper(classNumber))^2 )/ ( sum( posterior ) + EPS );
                    
                    means( classNumber ) = mu;
                    variances( classNumber ) = variance+EPS;
                    
                else
                    means( classNumber ) = meanHyper(classNumber);
                    variances( classNumber ) = 100;
                    
                end
            end
            variances(variances==0)=100; % added by Eugenio, prevents nans...
            
            
        end % End EM iterations
        means'
        (variances').^2
        
        %
        % Part II: update the position of the mesh nodes for the current set of Gaussian parameters
        %
        
        % Do the deformation one step at a time, for maximally positionUpdatingMaximumNumberOfIterations
        % deformation steps or until a step occurs in which the mesh node that moves most moves less than
        % maximalDeformationStopCriterion voxels, whichever comes first. The underlying algorithm is a
        % Levenberg-Marquardt type of algorithm in that it uses the gradient and an approximation of the
        % Hessian to propose a new position. If the new position proposal degrades the cost function (i.e.,
        % the posterior probability of the mesh node positions given the data and the parameters of the
        % imaging model goes down), the Hessian approximation is repeatedly altered by multiplying its diagonal
        % elements with an increasing factor, thereby making the proposal more and more gradient-descent
        % like with smaller-and-smaller step sizes, until a (small) position proposal is obtained that actually
        % improves the cost function. Conversely, every time a good position proposal is obtained, the
        % multiplication of the diagonal elements of the Hessian approximation is decreased the next time
        % around, making the algorithm much more efficient compared to gradient-descent (i.e., take much
        % larger step sizes) whenever it is possible.
        %
        % If no position proposal can be made even when the multiplication factor of the Hessian approximation's
        % diagonal becomes very large, i.e., even when the proposal is a tiny tiny deformation only, the
        % mesh node optimization algorithm gives up and tells you it didn't do anything.
        %
        %
        % NOTE: recall that this procedure is really only one half of a global optimization problem
        % that includes estimating the imaging model parameters (i.e., Gaussian intensity as well.
        % Therefore, it may not make sense to wait 20 minutes to get a really good optimization
        % of the mesh node positions here, as the cost function we're optimizing will change anyway
        % once the imaging model parameters are updated in the next iterations. Since updating the
        % imaging model parameters is very fast compared to updating the mesh nodes, it probably makes
        % sense to re-estimate the imaging model parameters frequently after a partial (not full)
        % optimization of the mesh nodes.
        %
        
        % Eugenio November 2017: GEMS2, more iterations
        
        if ( exist( 'optimizer', 'var' ) == 1 )
            % The optimizer is very memory hungry when run in multithreaded mode.
            % Let's clear any old ones we may have lying around
            try
                kvlClear( optimizer );
                kvlClear( calculator );
            catch ME
            end
        end
        
        haveMoved = false; % Keep track if we've ever moved or not
        positionUpdatingMaximumNumberOfIterations=30;
        %  (note that it uses variances instead of precisions)
        calculator = kvlGetCostAndGradientCalculator('AtlasMeshToIntensityImage',...
            image, 'Sliding',transform,means',variances',ones(size(means')),ones(size(means')));
        
        verbose=0;
        maximalDeformationStopCriterion=1e-10;
        lineSearchMaximalDeformationIntervalStopCriterion=1e-10;
        maximumNumberOfDeformationIterations=1000;
        BFGSMaximumMemoryLength=12;
        
        optimizer = kvlGetOptimizer( optimizerType, mesh, calculator, ...
            'Verbose', verbose, ...
            'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
            'LineSearchMaximalDeformationIntervalStopCriterion', ...
            lineSearchMaximalDeformationIntervalStopCriterion, ...
            'MaximumNumberOfIterations', maximumNumberOfDeformationIterations, ...
            'BFGS-MaximumMemoryLength', BFGSMaximumMemoryLength );


        for positionUpdatingIterationNumber = 1 : positionUpdatingMaximumNumberOfIterations
            % Calculate a good step. The first one is very slow because of various set-up issues
            disp(['Resolution level ' num2str(multiResolutionLevel) ' iteration ' num2str(iterationNumber) ' deformation iterations ' num2str(positionUpdatingIterationNumber)]);
            % Eugenio May2018
            maximalDeformation=0;
            try
                tic
                % Eugenio November 2017: GEMS2
                [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer);
                elapsedTime = toc;
                disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
                minLogLikelihoodTimesPrior
            end
            if isnan(minLogLikelihoodTimesPrior)
                error('lhood is nan');
            end
            if ( maximalDeformation > 0 )
                haveMoved = true;
            end
            
            % Test if we need to stop
            if ( maximalDeformation <= maximalDeformationStopCriterion )
                disp( 'maximalDeformation is too small; stopping' );
                break;
            end
        end
        
        % Show a little movie comparing before and after deformation so far...
        if  ~isdeployed && DEBUG>0
            figure( deformationMovieFigure )
            newColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, imageSize ) );
            for i=1:10
                showImage( oldColorCodedPriors ), title('Previous priors')
                drawnow
                pause( 0.1 )
                showImage( newColorCodedPriors ), title('Updated priors')
                drawnow
                pause( 0.1 )
            end
        end
        
        % Let's write this to file. In Unix/Linux systems, you can visualize progress over iterations
        % by doing animate -delay 10 *.png
        if  ~isdeployed && DEBUG>0
            fileName = [ 'colorCodedPrior_multiResolutionLevel' num2str( multiResolutionLevel ) '_iteration' sprintf( '%03d', iterationNumber ) '.png' ];
            tmp = getframe( gcf );
            imwrite( tmp.cdata, fileName );
        end
        
        
        % Keep track of the cost function we're optimizing
        historyOfCost = [ historyOfCost; minLogLikelihoodTimesPrior ];
        if  ~isdeployed && DEBUG>0
            figure( costFigure )
            plot( historyOfCost( 2 : end ) )
            title( 'Cost' )
        end
        
        % Determine if we should stop the overall iterations over the two set of parameters
        if ( ( ~haveMoved ) || ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / historyOfCost( end ) ) <  1e-6 ) )
            % Converged
            break;
        end
        
        
    end % End looping over global iterations
    
    
    
end % End loop over multiresolution levels

disp(['Fitting mesh to image data took ' num2str(etime(clock,time_ref_optimization)) ' seconds']);


% Restore original image buffer
kvlSetImageBuffer(image,imageBufferOrig);

% OK, now that all the parameters have been estimated, segment the image with all the original
% labels instead of the reduced "super"-structure labels we created.

% Clear some memory
% Eugenio November 2017: GEMS2
kvlClear( optimizer )
kvlClear( calculator )

%
if  ~isdeployed && DEBUG>0
    figure
    subplot( 2, 1, 1 )
    showImage( kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( mesh, imageSize )  ) );
    title('Priors')
    subplot( 2, 1, 2 )
    showImage( imageBuffer )
    title('Image')
end

% Undo the collapsing of several structures into "super"-structures
kvlSetAlphasInMeshNodes( mesh, originalAlphas )
numberOfClasses = size( originalAlphas, 2 );

% Eugenio July 2017 : this type of call is exactly what we're tryint to
% avoid....
data = double( reshape( imageBuffer, [ prod( imageSize ) 1 ] ) ); % Easier to work with vector notation in the computations
data = data( maskIndices );

posteriors=zeros([length(maskIndices),numberOfClasses],'single');
for classNumber = 1 : numberOfClasses
    prior = kvlRasterizeAtlasMesh( mesh, imageSize, classNumber-1 );
    mu = means( reducingLookupTable( classNumber ) );
    variance = variances( reducingLookupTable( classNumber ) );
    posteriors( :, classNumber ) = ( exp( -( data - mu ).^2 / 2 / variance ) ...
        .* (double(prior(maskIndices))/65535) ) / sqrt( 2 * pi * variance);
end
normalizer = sum( posteriors, 2 ) + eps;
posteriors = bsxfun(@rdivide,posteriors, normalizer);
posteriors = uint16( round( posteriors * 65535 ) );

% % Get the priors as dictated by the current mesh position
% data = double( reshape( imageBuffer, [ prod( imageSize ) 1 ] ) ); % Easier to work with vector notation in the computations
% priors = kvlRasterizeAtlasMesh( mesh, imageSize );
% priors = reshape( priors, [ prod( imageSize ) numberOfClasses ] );  % Easier to work with vector notation in the computations
% 
% % Ignore everything that's has zero intensity
% priors = priors( maskIndices, : );
% data = data( maskIndices );
% 
% % Calculate the posteriors
% posteriors = zeros( size( priors ), 'double' );
% for classNumber = 1 : numberOfClasses
%     % Get the parameters from the correct Gaussian
%     mu = means( reducingLookupTable( classNumber ) );
%     variance = variances( reducingLookupTable( classNumber ) );
%     prior = single( priors( :, classNumber ) ) / 65535;
%     
%     posteriors( :, classNumber ) = ( exp( -( data - mu ).^2 / 2 / variance ) .* prior ) ...
%         / sqrt( 2 * pi * variance);
% end
% normalizer = sum( posteriors, 2 ) + eps;
% posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfClasses ] );
% posteriors = uint16( round( posteriors * 65535 ) );

% Display the posteriors
if  ~isdeployed && DEBUG>0
    figure
    subplot( 1, 2, 1 )
    showImage( imageBuffer ); title('Image')
    subplot( 1, 2, 2 )
    tmp = zeros( prod(imageSize), numberOfClasses, 'uint16' );
    tmp( maskIndices, : ) = priors;
    tmp = reshape( tmp, [ imageSize numberOfClasses ] );
    showImage( kvlColorCodeProbabilityImages( tmp, colors ) ); title('Warped priors')
    
    figure
    for classNumber = 1 : numberOfClasses
        
        % Let's create a color image overlaying the posterior on top of gray scale
        overlayColor = [ 0 0 1 ]; % RGB value
        overlay = double( exp( imageBuffer / 1000 )  );
        overlay = overlay - min( overlay(:) );
        overlay = overlay / max( overlay(:) );
        posterior = single( tmp( :, :, :, classNumber ) / 65535 );
        overlay = overlay .* ( 1 - posterior );
        overlay = repmat( overlay, [ 1 1 1 3 ] );
        colorPosterior = reshape( posterior(:) * overlayColor, [ imageSize 3 ] );
        showImage( overlay + colorPosterior );
        
        name = deblank( names( classNumber, : ) );
        title( name, 'interpreter', 'none' )
        
        
        % disp( 'Enter to continue' )
        % pause
        drawnow
        pause( 0.5 )
        
    end
end

% Show the histogram of each structure, where individual voxels are weight according
% to their posterior of belonging to that structure. Overlay on top of it the Gaussian
% intensity model estimated with the EM procedure above - any clear deviations between
% the histogram and the Gaussian model indicate that the model is a bad fit for that
% structure, and that more accurate image modeling may be needed there.
if  ~isdeployed && DEBUG>0
    range = [ min( data(:) ) max( data(:) ) ];
    numberOfBins = 100;
    binEdges = [ 0 : numberOfBins ] / numberOfBins;
    binEdges = binEdges * ( range(2) - range(1) ) + range( 1 );
    binDelta = binEdges(2) - binEdges(1);
    binCenters = binEdges( 1 : numberOfBins ) + binDelta/2;
    [ totalHistogram, binNumbers ] = histc( data, binEdges );
    figure
    bar( binCenters, totalHistogram(1:end-1) )
    figure
    for classNumber = 1 : numberOfClasses
        posterior = double( posteriors( :, classNumber ) / 65535 );
        
        histogram = zeros( numberOfBins, 1 );
        for binNumber = 1 : numberOfBins % Ough this is hugely inefficient in Matlab
            histogram( binNumber ) = sum( posterior( find( binNumbers == binNumber ) ) );
        end
        histogram = histogram / sum( histogram );
        bar( binCenters, histogram )
        hold on
        mu = means( reducingLookupTable( classNumber ) );
        variance = variances( reducingLookupTable( classNumber ) );
        gauss = exp( -( binCenters - mu ).^2 / 2 / variance ) / sqrt( 2 * pi * variance );
        p = plot( binCenters, gauss * binDelta, 'r' );
        hold off
        name = deblank( names( classNumber, : ) );
        title( name, 'interpreter', 'none' )
        
        
        %         disp( 'Enter to continue' )
        %         pause
        drawnow
        pause( 0.5 )
    end
end

% Write the resulting atlas mesh, as well as the image we segmented, to file for future
% reference.
%
kvlWriteMeshCollection( meshCollection, 'warpedMesh.txt' ); % The warped mesh will have index 0, and will be

transformMatrix = double(kvlGetTransformMatrix( transform ));
inverseTransform = kvlCreateTransform( inv( transformMatrix ) );
kvlTransformMeshCollection( meshCollection, inverseTransform );
kvlWriteMeshCollection( meshCollection, 'warpedMeshNoAffine.txt' );



kvlWriteImage( image, 'image.mgz' );
system(['cp ' compressionLookupTableFileName ' .']);


% You can now use the C++ tools distributed with FreeSufer to inspect what we have as
% follows:
%
%    kvlViewMeshCollectionWithGUI warpedMesh.txt.gz 86 113 163 image.mgz
%
% where 86 113 163 are the dimensions of image.mgz, which you need to manually specify (don't ask!)
%


% Eugenio: write discrete labels (MAP)
% Note how we fill in the gaps in the regions outside the FOV with the
% prior!


% First we compute the shift when using kvlReadCroppedImage/kvlWrite
if exist('additionalFile','var')>0
    system([FSdir '/mri_convert asegAddSpace.mgz asmr1.mgz -rl T1resampled.mgz -rt nearest >/dev/null']);
else
    system([FSdir '/mri_convert ' subjectDir '/' subjectName '/mri/aseg.mgz asmr1.mgz -rl T1resampled.mgz -rt nearest >/dev/null']);
end
[ asmr, tasmr ] = kvlReadCroppedImage( 'asmr1.mgz', boundingFileName );
asmrB=kvlGetImageBuffer(asmr);
kvlWriteImage( asmr, 'asmr2.mgz' );
tmp1=myMRIread('asmr1.mgz',0,tempdir);
tmp2=myMRIread('asmr2.mgz',0,tempdir);
tmp3=tmp1; tmp3.vol=tmp2.vol;
myMRIwrite(tmp3,'asmr3.mgz','float',tempdir);
[I,J,K]=ind2sub(size(tmp1.vol),find(tmp1.vol==16));
Ic1=mean(I); Jc1=mean(J); Kc1=mean(K);
[I,J,K]=ind2sub(size(tmp2.vol),find(tmp2.vol==16));
Ic2=mean(I); Jc2=mean(J); Kc2=mean(K);
shift=round([Ic2-Ic1,Jc2-Jc1,Kc2-Kc1]);
shiftNeg=-shift; shiftNeg(shiftNeg<0)=0;
shiftPos=shift; shiftPos(shiftPos<0)=0;

% reorient image.mgz
aux=zeros(size(tmp2.vol)+shiftNeg);
aux2=myMRIread('image.mgz',0,tempdir);
% aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=aux2.vol;
aux(1+shiftNeg(1):shiftNeg(1)+size(aux2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(aux2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(aux2.vol,3))=aux2.vol;
aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
tmp3.vol=aux;
myMRIwrite(tmp3,'image.mgz','float',tempdir);


% Compute posteriors and volumes
% July 2017: make this memory efficient, and add new substructures. We also
% compute the segmentation along the way. 
fid=fopen([tempdir '/volumesBrainstem.txt'],'w');
strOfInterest={'MAC_Medulla','MAC_Pons','MAC_Midbrain','MAC_Sup_Cerebellum_Ped','Medulla','Pons','Midbrain','SCP'};
totVol=0;
found=zeros(1,numberOfClasses);
for i=1:numberOfClasses
    
    if i==1
        sillyAlphas=zeros([size(originalAlphas,1),2],'single');
        sillyAlphas(:,1)=originalAlphas(:,1);
        sillyAlphas(:,2)=1-sillyAlphas(:,1);
        kvlSetAlphasInMeshNodes( mesh, sillyAlphas )
        post = kvlRasterizeAtlasMesh( mesh, imageSize);
        post=post(:,:,:,1);
        kvlSetAlphasInMeshNodes( mesh, originalAlphas );
    else
        post=kvlRasterizeAtlasMesh( mesh, imageSize , i-1);
    end
    post(maskIndices)=posteriors(:,i);
    
    if i==1
        L=ones(size(post));
        MAXP=post;
    else
        M=post>MAXP;
        L(M)=i;
        MAXP(M)=post(M);
    end
    
    found(i)=0;
    str=[];
    vol=0;
    name=names(i,:);
    name=lower(name(name~=' '));
    for j=1:length(strOfInterest)
        if strcmp(name,lower(strOfInterest{j}))>0
            found(i)=j;
        end
    end
    if found(i)>0
        str=strOfInterest{found(i)};
        vol=resolution^3*(sum(double(post(:))/65535));
        fprintf(fid,'%s %f\n',str,vol);
        totVol=totVol+vol;
        if WRITE_POSTERIORS>0
            kk1=double(post)/65535;
            kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
            aux=zeros(size(tmp2.vol)+shiftNeg);
            aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
            aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
            tmp3.vol=aux;
            if exist('additionalFile','var')>0
                myMRIwrite(tmp3,['posterior_' strtrim(lower(names(i,:))) '_' suffix '.additional.mgz'],'float',tempdir);
            else
                myMRIwrite(tmp3,['posterior_' strtrim(lower(names(i,:))) '_' suffix '.mgz'],'float',tempdir);
            end
        end
    end
end
fprintf(fid,'Whole_brainstem %f\n',totVol);
fclose(fid);

% priorsFull = kvlRasterizeAtlasMesh( mesh, imageSize );
% posteriorsFull=priorsFull;
% 
% fid=fopen([tempdir '/volumesBrainstem.txt'],'w');
% strOfInterest={'MAC_Medulla','MAC_Pons','MAC_Midbrain','MAC_Sup_Cerebellum_Ped','Medulla','Pons','Midbrain','SCP'};
% totVol=0;
% found=zeros(1,size(priorsFull,4));
% for i=1:size(priorsFull,4)
%     tmp=posteriorsFull(:,:,:,i);
%     tmp(maskIndices)=posteriors(:,i);
%     posteriorsFull(:,:,:,i)=tmp;
%     found(i)=0;
%     str=[];
%     vol=0;
%     name=names(i,:);
%     name=lower(name(name~=' '));
%     for j=1:length(strOfInterest)
%         if strcmp(name,lower(strOfInterest{j}))>0
%             found(i)=j;
%         end
%     end
%     if found(i)>0
%         str=strOfInterest{found(i)};
%         vol=resolution^3*(sum(sum(sum(double(posteriorsFull(:,:,:,i))/65535))));
%         fprintf(fid,'%s %f\n',str,vol);
%         totVol=totVol+vol;
%         if WRITE_POSTERIORS>0
%             kk1=double(posteriorsFull(:,:,:,i))/65535;
%             kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
%             aux=zeros(size(tmp2.vol)+shiftNeg);
%             aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
%             aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
%             tmp3.vol=aux;
%             if exist('additionalFile','var')>0
%                 myMRIwrite(tmp3,['posterior_' strtrim(lower(names(i,:))) '_' suffix '.additional.mgz'],'float',tempdir);
%             else
%                 myMRIwrite(tmp3,['posterior_' strtrim(lower(names(i,:))) '_' suffix '.mgz'],'float',tempdir);
%             end
%         end
%     end
% end
% fprintf(fid,'Whole_brainstem %f\n',totVol);
% fclose(fid);


% MAP estimates

% Eugenio July 2011
% [~,inds]=max(posteriorsFull,[],4);
inds=L;

kk1=FreeSurferLabels(inds);
kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
aux=zeros(size(tmp2.vol)+shiftNeg);
% aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
aux(1+shiftNeg(1):shiftNeg(1)+size(aux2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(aux2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(aux2.vol,3))=permute(kk1,[2 1 3]);
aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
tmp3.vol=aux;
myMRIwrite(tmp3,'discreteLabels_all.mgz','float',tempdir);
tmp3.vol(tmp3.vol<170)=0;
tmp3Mask=getLargestCC(tmp3.vol>0);
tmp3.vol(~tmp3Mask)=0;
myMRIwrite(tmp3,'discreteLabels.mgz','float',tempdir);




if DEBUG>0
    
    aux=zeros([size(tmp2.vol)+shiftNeg size(priorsFull,4)]);
    aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3),:)=permute(priorsFull,[2 1 3 4]);
    aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end,:);
    tmp3.vol=aux;
    myMRIwrite(tmp3,'priors.mgz','float',tempdir);
    
    aux=zeros([size(tmp2.vol)+shiftNeg size(posteriorsFull,4)]);
    aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3),:)=permute(posteriorsFull,[2 1 3 4]);
    aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end,:);
    tmp3.vol=aux;
    myMRIwrite(tmp3,'posteriors.mgz','float',tempdir);
    
end

if DEBUG>0
    % Interesting for debugging: likelihood terms
    LOGLHOODS=zeros(size(posteriorsFull),'double');
    for i=1:length(means)
        gauss=-.5*log(2*pi*variances(i))-.5*(double(imageBuffer(maskIndices))-means(i)).^2/variances(i);
        kkk=zeros(imageSize);
        kkk(maskIndices)=gauss;
        LOGLHOODS(:,:,:,i)=double(kkk);
    end
    aux=zeros([size(tmp2.vol)+shiftNeg size(LOGLHOODS,4)]);
    aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3),:)=permute(LOGLHOODS,[2 1 3 4]);
    aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end,:);
    tmp3.vol=aux;
    myMRIwrite(tmp3,'loglikelihoods.mgz','float',tempdir);
    tmp3.vol=exp(tmp3.vol);
    myMRIwrite(tmp3,'likelihoods.mgz','float',tempdir);
end


% Convert to 1 mm FreeSurfer Space
if SMOOTH_LABEL_RESAMPLE>0
    refFile=[subjectDir '/' subjectName '/mri/norm.mgz'];
    applyLTAsmoothLabels('discreteLabels.mgz',[],'discreteLabelsResampledT1.mgz',refFile,0,FSpath,tempdir);
else
    system([FSdir '/mri_convert  discreteLabels.mgz  discreteLabelsResampledT1.mgz -rt nearest -odt float ' ...
        ' -rl ' subjectDir '/' subjectName '/mri/norm.mgz']);
end
% Move to MRI directory
if exist('additionalFile','var')>0
    system(['mv discreteLabels.mgz ' subjectDir '/' subjectName '/mri/brainstemSsLabels.' suffix '.additional.mgz']);
    system(['mv discreteLabelsResampledT1.mgz ' subjectDir '/' subjectName '/mri/brainstemSsLabels.' suffix '.additional.FSvoxelSpace.mgz']);
    system(['mv volumesBrainstem.txt ' subjectDir '/' subjectName '/mri/brainstemSsVolumes.' suffix '.additional.txt']);
else
    system(['mv discreteLabels.mgz ' subjectDir '/' subjectName '/mri/brainstemSsLabels.' suffix '.mgz']);
    system(['mv discreteLabelsResampledT1.mgz ' subjectDir '/' subjectName '/mri/brainstemSsLabels.' suffix '.FSvoxelSpace.mgz']);
    system(['mv volumesBrainstem.txt ' subjectDir '/' subjectName '/mri/brainstemSsVolumes.' suffix '.txt']);
end

if WRITE_POSTERIORS>0
    system(['mv posterior_*.mgz ' subjectDir '/' subjectName '/mri/']);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Clean up, I guess... you might wanna skip this with debugging purposes... %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cmd=['rm -r -f ' tempdir];
system([cmd ' >/dev/null']);


disp('Everything done!')
disp(['It took ' num2str(etime(clock,time_start)) ' seconds ']);


if isdeployed
    exit
end


