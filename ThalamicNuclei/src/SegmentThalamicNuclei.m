% Segments the thalamic nuclei from a FS-processed T1 scan (or an
% additional volumes)
%
% This is heavily based on segmentSubjectT1_autoEstimateAlveusML from the
% hippocampal subfield code
%
% SegmentThalamicNuclei(subjectName,subjectDir,resolution,atlasMeshFileName,...
%    atlasDumpFileName,compressionLUTfileName,K,optimizerType,suffix,FSpath,...
%    useTwoComponents,MRFconstant,additionalVol,analysisID,doBFcorrection,BBregisterMode)
%
% - subjectName: FreeSurfer subject name
% - subjectDir: FreeSurfer subject directory
% - resolution: voxel size at which we want to work (in mm).
% - atlasMeshFileName: the atlas to segment the datafbbre
% - atlasDumpFileName: corresponding imageDump.mgz (name *must* be imageDump.mgz)
% - compressionLUTfileName: corresponding compressionLUT.txt
% - K: stiffness of the mesh in the segmentation.
% - optimizerType: 'FixedStepGradientDescent','GradientDescent','ConjugateGradient','L-BFGS'
% - suffix: for output directory, e.g. 'v1.0'
% - FSpath: path to FreeSurfer executables (FS's 'bin' directory)
% - useTwoComponents: use two Gaussians to model thalamus intensities (should always be 1)
% - MRFconstant: to smooth result
% - additionalVol: file name of a volume to replace the T1 in the analysis
% - analysisID: and ID for the analysis. By default, it's T1. But if you
%      use an additional volume, you can use T2, FGATIR, etc.
% - doBFcorrection: set to 1 if you want to bias field correct the additional volaume
% - BBregisterMode: set to 't1' if you have t1-like contrast, 't2' if you have
%     t2-like contrast, and 'none' if you are providing an already registered
%     additional volume
% 
%
function SegmentThalamicNuclei(subjectName,subjectDir,resolution,atlasMeshFileName,...
    atlasDumpFileName,compressionLUTfileName,K,optimizerType,suffix,FSpath,...
    useTwoComponents,MRFconstant,additionalVol,analysisID,doBFcorrection,BBregisterMode)

DEBUG=0;
FAST=0; % set it to one to optimize just a bit (go through code fast)
WRITE_POSTERIORS=0;
WRITE_MESHES=0;
SMOOTH_LABEL_RESAMPLE=0;
THALAMUS_VERBOSE=0;

% March 2022: fix to accommodate 'fs_run_from_mcr'
FSpath = [FSpath '/fs_run_from_mcr ' FSpath '/'];

aux=getenv('WRITE_POSTERIORS');
if ~isempty(aux)
    if str2double(aux)>0
        WRITE_POSTERIORS=1;
    end
end
aux=getenv('WRITE_MESHES');
if ~isempty(aux)
    if str2double(aux)>0
        WRITE_MESHES=1;
    end
end
aux=getenv('THALAMUS_VERBOSE');
if ~isempty(aux)
    if str2double(aux)>0
        THALAMUS_VERBOSE=1;
    end
end

% sanity check
if exist('MRFconstant','var')==0
    MRFconstant=0;
end


if exist('additionalVol','var')==0
    additionalVol=[]; analysisID='T1'; doBFcorrection=0; BBregisterMode='none';
elseif exist(additionalVol,'file')==0
    error('Additional volume does not exist');
elseif exist('analysisID','var')==0
    error('If you provide an additional volume, you also need an analysis ID');
elseif exist('doBFcorrection','var')==0
    error('If you provide an additional volume, you also need to specify if you want to bias field correct it');
elseif exist('BBregisterMode','var')==0
    error('If you provide an additional volume, you also need to specify the contrast for registration with BBregister');
end



if nargin<11
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
elseif ~isdeployed && (~isnumeric(useTwoComponents))
    error('useTwoComponents must be numeric');
elseif ~isdeployed && (~isnumeric(MRFconstant))
    error('MRFconstant must be numeric');
elseif ~isdeployed && (~isnumeric(doBFcorrection))
    error('doBFcorrection must be numeric');
elseif strcmp(lower(BBregisterMode),'none')==0 &&  strcmp(lower(BBregisterMode),'t1')==0  && strcmp(lower(BBregisterMode),'t2')==0 
    error('BBregisterMode must be none, t1 or t2');
end


% Constants
THlabelLeft=10;
THlabelRight=49;
DElabelLeft=28;
DElabelRight=60;

% In case we compiled it...
if isdeployed
    K=str2double(K);
    resolution=str2double(resolution);
    useTwoComponents=str2double(useTwoComponents);
    MRFconstant=str2double(MRFconstant);
    doBFcorrection=str2double(doBFcorrection);
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
tempdir=[subjectDir '/' subjectName '/tmp/THnuclei_' suffix '_useTwoComponents_' num2str(useTwoComponents) '_' analysisID '/'];
aux=getenv('USE_SCRATCH');
if ~isempty(aux)
    if str2double(aux)>0
        if exist('/scratch/','dir')>0
            tempdir=['/scratch/' subjectName '_THnuclei_' suffix '_useTwoComponents_' num2str(useTwoComponents) '_' analysisID '/'];
        end
    end
end
if exist(tempdir,'dir')==0
    mkdir(tempdir);
end


% First thing: registration of additional volume
if ~isempty(additionalVol)
    mov=additionalVol;
    registeredAddVol=[subjectDir '/' subjectName '/mri/' analysisID '.thalamus.' BBregisterMode '.mgz'];
    lta=[tempdir '/additional.lta'];
    
    if strcmp(lower(BBregisterMode),'none')
        disp('Resampling additional volume (skipping registration)');
        system([FSpath '/mri_convert ' mov ' ' registeredAddVol ' -odt float -rl ' subjectDir '/' subjectName '/mri/norm.mgz  >/dev/null']);
    elseif exist(registeredAddVol,'file')
        disp('Additional volume has already been registered/resampled; skipping registration');
    else
        disp('Registering additional volume');
        setenv('SUBJECTS_DIR',subjectDir);
        cmd=[FSpath '/bbregister --s ' subjectName  '  --mov ' mov  ...
            '  --reg ' lta ' --init-header --o ' registeredAddVol ' --' BBregisterMode];
        
        if THALAMUS_VERBOSE, system(cmd); else, system([cmd '  >/dev/null']);  end
        
        % Create animated gif for QC of registration
        auxT1=myMRIread([subjectDir '/' subjectName '/mri/nu.mgz'],0,tempdir);
        auxSEG=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
        auxADD=myMRIread(registeredAddVol,0,tempdir);
        [~,~,auxSl]=ind2sub(size(auxSEG.vol),find(auxSEG.vol==THlabelLeft | auxSEG.vol==THlabelRight));
        Sl=round(mean(auxSl));
        IT1=auxT1.vol(:,:,Sl); IT1=IT1/max(IT1(:))*255; IT1=uint8(IT1);
        ITadd=auxADD.vol(:,:,Sl); ITadd=ITadd/max(ITadd(:))*255; ITadd=uint8(ITadd);
        
        flipFile=[subjectDir '/' subjectName '/mri/transforms/T1-' analysisID '.' suffix '.QC.gif'];
        if exist(flipFile,'file')>0, delete(flipFile); end;
        [imind,cm] = gray2ind(IT1,256);
        imwrite(imind,cm,flipFile,'gif', 'Loopcount',inf);
        [imind,cm] = gray2ind(ITadd,256);
        imwrite(imind,cm,flipFile,'gif','WriteMode','append');
    
        system(['mv ' lta ' ' subjectDir '/' subjectName '/mri/transforms/' analysisID '.' suffix '.lta']);
        
        disp('Registration done!');
    end
end

cd(tempdir);

% Next: register image dump to automated segmentation
disp('Registering imageDump.mgz to mask from ASEG')


% Grab atlas (manually placed in FS atlas coordinate space!)
system(['cp ' atlasDumpFileName ' ./imageDump.mgz']);


% Target is masked aseg (if
targetRegFileName=[tempdir '/THDEbinaryMask.mgz'];
targetRegFileNameCropped=[tempdir '/THDEbinaryMask_autoCropped.mgz'];
ASEG=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
TARGETREG=ASEG;
TARGETREG.vol=255*double(ASEG.vol==THlabelLeft | ASEG.vol==THlabelRight | ASEG.vol==DElabelLeft | ASEG.vol==DElabelRight);
myMRIwrite(TARGETREG,targetRegFileName,'float',tempdir);

highres=0; if mean(ASEG.volres)<0.99, highres=1; end
if highres==1,
    system([FSpath '/mri_convert ' targetRegFileName ' aux.mgz -odt float -vs 1 1 1 -rt nearest >/dev/null']);
    system(['mv aux.mgz ' targetRegFileName ' >/dev/null']);
end


% Replacement of kvlAutoCrop
aux=myMRIread(targetRegFileName,0,tempdir);
[aux.vol,cropping]=cropLabelVol(aux.vol,6);
shift=aux.vox2ras0(1:3,1:3)*[cropping(2)-1; cropping(1)-1; cropping(3)-1];
aux.vox2ras0(1:3,4)=aux.vox2ras0(1:3,4)+shift;
aux.vox2ras1(1:3,4)=aux.vox2ras1(1:3,4)+shift;
aux.vox2ras(1:3,4)=aux.vox2ras(1:3,4)+shift;
aux.tkrvox2ras=[];
myMRIwrite(aux,targetRegFileNameCropped,'float',tempdir);


if 1==1  % This is to use an opened version (closed is better?)
    aux=myMRIread(targetRegFileNameCropped,0,tempdir);
    strel=createSphericalStrel(1);
    % aux.vol=255*double(imdilate(imerode(aux.vol>0,strel),strel));
    aux.vol=255*double(imerode(imdilate(aux.vol>0,strel),strel));
    myMRIwrite(aux,targetRegFileNameCropped,'float',tempdir);
end

% Registration
cmd=[FSpath '/mri_robust_register --mov imageDump.mgz  --dst ' targetRegFileNameCropped ...
    ' -lta trash.lta --mapmovhdr imageDump_coregistered.mgz  --sat 50'];
if THALAMUS_VERBOSE, status=system(cmd); else, status=system([cmd '  >/dev/null']);  end
if status~=0, error('Problem with mri_robust_register'); end
system('mv imageDump_coregistered.mgz imageDump.mgz' );

cmd=[FSpath '/mri_robust_register --mov imageDump.mgz  --dst ' targetRegFileNameCropped ...
    ' -lta trash.lta --mapmovhdr imageDump_coregistered.mgz --affine --sat 50'];
if THALAMUS_VERBOSE, status=system(cmd); else, status=system([cmd '  >/dev/null']);  end
if status~=0, error('Problem with mri_robust_register'); end
system('mv imageDump_coregistered.mgz imageDump.mgz' );




% Now, the idea is to refine the transform based on the thalamus + VDE
% First, we prepare a modifided ASEG that we'll segment

% There's a bunch of labels in the ASEG don't have in our atlas...
ASEGbackup=ASEG;
ASEG.vol(ASEG.vol==5)=4;   % left-inf-lat-vent -> left-lat-vent
ASEG.vol(ASEG.vol==44)=4;  % right-inf-lat-vent -> left-lat-vent
ASEG.vol(ASEG.vol==14)=4;  % 3rd vent -> left-lat-vent
ASEG.vol(ASEG.vol==15)=4;  % 4th vent -> LV (we're killing brainstem anyway...)
ASEG.vol(ASEG.vol==17)=3; % left HP -> left cortex
ASEG.vol(ASEG.vol==53)=3; % right HP -> left cortex
ASEG.vol(ASEG.vol==18)=3; % left amygdala -> left cortex
ASEG.vol(ASEG.vol==54)=3; % right amygdala -> left cortex
ASEG.vol(ASEG.vol==24)=4;  % CSF -> left-lat-vent
ASEG.vol(ASEG.vol==30)=2;  % left-vessel -> left  WM
ASEG.vol(ASEG.vol==62)=2; % right-vessel -> left  WM

% right to left
ASEG.vol(ASEG.vol==41)=2; % WM
ASEG.vol(ASEG.vol==42)=3; % CT
ASEG.vol(ASEG.vol==43)=4; % LV
ASEG.vol(ASEG.vol==46)=7;  % cerebellum WM
ASEG.vol(ASEG.vol==47)=8;  % cerebellum CT
ASEG.vol(ASEG.vol==50)=11; % CA
ASEG.vol(ASEG.vol==51)=12; % PU
ASEG.vol(ASEG.vol==52)=13; % PA
ASEG.vol(ASEG.vol==58)=26; % AA
ASEG.vol(ASEG.vol==63)=31; % CP

ASEG.vol(ASEG.vol==72)=4;  % 5th ventricle -> left-lat-vent
ASEG.vol(ASEG.vol==77)=2;  % WM hippoint -> left WM
ASEG.vol(ASEG.vol==80)=0;  % non-WM hippo -> background
ASEG.vol(ASEG.vol==85)=0;  % optic chiasm -> background

ASEG.vol(ASEG.vol>250)=2;  % CC labels -> left WM

list2kill=[44 62 63 41 42 43 50 51 52 53 54 58];
for k=1:length(list2kill)
    ASEG.vol(ASEG.vol==list2kill(k))=0;
end

ASEG.vol(ASEG.vol==0)=1;



% Write to disk
ASEGmod=ASEG;
myMRIwrite(ASEGmod,'asegMod.mgz','float',tempdir);

ASEG_TH_DE=ASEG;
ASEG_TH_DE.vol(ASEG_TH_DE.vol==DElabelLeft)=THlabelLeft;
ASEG_TH_DE.vol(ASEG_TH_DE.vol==DElabelRight)=THlabelRight;
myMRIwrite(ASEG_TH_DE,'aseg_TH_DE.mgz','float',tempdir);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we pretty much copy-paste from preprocessHippoSubfields %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

asegFileName = [subjectDir '/' subjectName '/mri/aseg.mgz'];  % FreeSurfer's volumetric segmentation results. This are non-probabilistic, "crisp" definitions
boundingFileName = [tempdir 'imageDump.mgz']; % Bounding box
meshCollectionFileName = atlasMeshFileName; % The tetrahedral atlas mesh
compressionLookupTableFileName =compressionLUTfileName; % Look-up table belonging to the atlas

%
[ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );

% Read in aseg and also transform
[ synIm, transform ] = kvlReadCroppedImage( 'aseg_TH_DE.mgz', boundingFileName );
synImBuffer = kvlGetImageBuffer( synIm );
synSize = size( synImBuffer );
if ~isdeployed && DEBUG>0
    figure
    showImage( synIm )
    title('Synthetic Image')
end


% read in collection, set K and apply transform
meshCollection = kvlReadMeshCollection( meshCollectionFileName );
kvlTransformMeshCollection( meshCollection, transform );
kvlSetKOfMeshCollection( meshCollection, K );

% Retrieve the reference mesh, i.e., the mesh representing the average shape.
mesh = kvlGetMesh( meshCollection, -1 );
originalNodePositions = kvlGetMeshNodePositions( mesh );
originalAlphas = kvlGetAlphasInMeshNodes( mesh );

% % Just for illustrative purposes, let's also display this mesh warped onto each
% % of the training subjects, as computed during the group-wise registration during
% % the atlas building
% if ~isdeployed && DEBUG>0
%     for meshNumber = 0 : 14  % C-style indexing
%         pause( .1 )
%         showImage( kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( kvlGetMesh( meshCollection, meshNumber ), synSize ), colors ) )
%         title(['Reference mesh warped back onto subject ' num2str(1+meshNumber)]);
%     end
% end
FreeSurferLabelGroups=[];
FreeSurferLabelGroups{end+1}={'Unknown'};
FreeSurferLabelGroups{end+1}={'Left-Cerebral-White-Matter'};
FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex'};
FreeSurferLabelGroups{end+1}={'Left-Cerebellum-Cortex'};
FreeSurferLabelGroups{end+1}={'Left-Cerebellum-White-Matter'};
FreeSurferLabelGroups{end+1}={'Brain-Stem'};
FreeSurferLabelGroups{end+1}={'Left-Lateral-Ventricle'};
FreeSurferLabelGroups{end+1}={'Left-choroid-plexus'};
FreeSurferLabelGroups{end+1}={'Left-Putamen'};
FreeSurferLabelGroups{end+1}={'Left-Pallidum'};
FreeSurferLabelGroups{end+1}={'Left-Accumbens-area'};
FreeSurferLabelGroups{end+1}={'Left-Caudate'};

% For full model
FreeSurferLabelGroups{end+1}={'Left-L-Sg','Left-LGN','Left-MGN','Left-PuI','Left-PuM','Left-H','Left-PuL',...
    'Left-VPI','Left-PuA','Left-R','Left-MV(Re)','Left-Pf','Left-CM','Left-LP','Left-VLa','Left-VPL','Left-VLp',...
    'Left-MDm','Left-VM','Left-CeM','Left-MDl','Left-Pc','Left-MDv','Left-Pv','Left-CL','Left-VA','Left-VPM',...
    'Left-AV','Left-VAmc','Left-Pt','Left-AD','Left-LD','Left-VentralDC'};
FreeSurferLabelGroups{end+1}={'Right-L-Sg','Right-LGN','Right-MGN','Right-PuI','Right-PuM','Right-H','Right-PuL',...
    'Right-VPI','Right-PuA','Right-R','Right-MV(Re)','Right-Pf','Right-CM','Right-LP','Right-VLa','Right-VPL','Right-VLp',...
    'Right-MDm','Right-VM','Right-CeM','Right-MDl','Right-Pc','Right-MDv','Right-Pv','Right-CL','Right-VA','Right-VPM',...
    'Right-AV','Right-VAmc','Right-Pt','Right-AD','Right-LD','Right-VentralDC'};


sameGaussianParameters=[];
for g=1:length(FreeSurferLabelGroups)
    sameGaussianParameters{end+1} = [];
    for FreeSurferLabel =  FreeSurferLabelGroups{g}
        sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
    end
    if isempty(sameGaussianParameters{end})
        sameGaussianParameters=sameGaussianParameters(1:end-1);
    end
end


cheatingMeans=zeros(length( sameGaussianParameters),1);
cheatingVariances=0.01*ones(length( sameGaussianParameters),1);
for l=1:length(sameGaussianParameters)
    label= sameGaussianParameters{l}(1);
    if label>=8100 && label<8200,  cheatingMeans(l)=THlabelLeft; % left thalamic nuclei + DE-> left TH
    elseif label>=8200,  cheatingMeans(l)=THlabelRight; % right thalamic nuclei + DE -> left TH
    elseif label==0, cheatingMeans(l)=1; % BACKGROUND is 1 instead of 0
    else cheatingMeans(l)=label;
    end
end



% Compute the "reduced" alphas - those referring to the "super"-structures
[ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
if ( max( abs( sum( reducedAlphas, 2 ) - 1 ) ) > 1e-5 ) % Make sure these vectors really sum to 1
    error( 'The vector of prior probabilities in the mesh nodes must always sum to one over all classes' )
end

% Set the reduced alphas to be the alphas of the mesh
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )


% Eugenio July 2017

% priors = kvlRasterizeAtlasMesh( mesh, synSize );
% MASK=imerode(sum(double(priors)/65535,4)>0.99,createSphericalStrel(3));

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
MASK=imerode(sum(single( sumpriors / 65535 ),4)>0.99,createSphericalStrel(3));

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


% We use a multiscale approach here (Koen had a single one with sigma=3)
meshSmoothingSigmas = [ 3.0 2.0]';
% Eugenio November 2017: increased number of iterations by 50%
maxIt=[300,150];


numberOfMultiResolutionLevels = length( meshSmoothingSigmas );

time_ref_cheat_optimization=clock;

historyOfMinLogLikelihoodTimesPrior = [ 1/eps ];

for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
    
    % Smooth the mesh using a Gaussian kernel.
    % It's good to smooth the mesh, otherwise we get weird compressions of the
    % mesh along the boundaries...
    kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
    meshSmoothingSigma = meshSmoothingSigmas( multiResolutionLevel );
    fprintf( 'Smoothing mesh collection with kernel size %f ...', meshSmoothingSigma )
    kvlSmoothMeshCollection( meshCollection, meshSmoothingSigma )
    fprintf( 'done\n' )
    
    
    % Show the smoothed atlas
    if ~isdeployed && DEBUG>0
        figure
        priors = kvlRasterizeAtlasMesh( mesh, synSize );
        for cheatingLabel = 1 : size( reducedAlphas, 2 )
            subplot( 3, 4, cheatingLabel )
            showImage( priors( :, :, :, cheatingLabel ) )
        end
        title('Smoothed priors')
    end
    
    % Eugenio November 2017: GEMS2
    % Set up the black box optimizer for the mesh nodes
    if ( exist( 'cheatingOptimizer', 'var' ) == 1 )
        % The optimizer is very memory hungry when run in multithreaded mode.
        % Let's clear any old ones we may have lying around
        % Eugenio November 2017: GEMS2
        kvlClear( cheatingOptimizer );
        kvlClear( cheatingCalculator );
    end
    

    % Eugenio November 2017: GEMS2  (note that it uses variances instead of precisions)
    % Now the optimization per-se
    cheatingCalculator = kvlGetCostAndGradientCalculator('AtlasMeshToIntensityImage',...
        cheatingImage, 'Sliding',transform,cheatingMeans,cheatingVariances,ones(size(cheatingMeans)),ones(size(cheatingMeans)));
    
    verbose=0;
    maximalDeformationStopCriterion=1e-10;
    lineSearchMaximalDeformationIntervalStopCriterion=1e-10;
    maximumNumberOfDeformationIterations=1000;
    BFGSMaximumMemoryLength=12;
    
    % optimizer = kvlGetOptimizer( optimizerType, mesh, calculator);
    cheatingOptimizer = kvlGetOptimizer( optimizerType, mesh, cheatingCalculator, ...
        'Verbose', verbose, ...
        'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
        'LineSearchMaximalDeformationIntervalStopCriterion', ...
        lineSearchMaximalDeformationIntervalStopCriterion, ...
        'MaximumNumberOfIterations', maximumNumberOfDeformationIterations, ...
        'BFGS-MaximumMemoryLength', BFGSMaximumMemoryLength );
    
    relativeChangeInCostStopCriterion = 1e-10;
    maxpuin=maxIt(multiResolutionLevel);
    if FAST>0
        maxpuin=20;
    end
    
    
    for positionUpdatingIterationNumber = 1 : maxpuin
        if THALAMUS_VERBOSE>0 || mod(positionUpdatingIterationNumber,10)==1
            disp(['Resolution ' num2str(multiResolutionLevel) ', iteration ' num2str(positionUpdatingIterationNumber)]);
        end
        % Calculate a good step. The first one is very slow because of various set-up issues % Eugenio May2018
        maximalDeformation=0;
        try
            tic
            % Eugenio November 2017: GEMS2
            [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( cheatingOptimizer );
            elapsedTime = toc;
            if THALAMUS_VERBOSE
                disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
                minLogLikelihoodTimesPrior
            end
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
end
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
% This is nice because all we need to do is to modify
% imageDump_coregistered with the T1-to-T2 transform to have the warped
% mesh in T2 space :-)
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


% Eugenio: we extract a block from norm.mgz and upsample it to the work
% resolution
% We also mask out non-brain voxels and also the cerebellum,
% brainstem and 3rd/4th ventricles, which can be annoying later on.i
%
imageFileName='T1resampled.mgz';
margin=15; % in mm
if isempty(additionalVol)
    A=myMRIread([subjectDir '/' subjectName '/mri/norm.mgz'],0,tempdir);
elseif doBFcorrection==0
    A=myMRIread(registeredAddVol,0,tempdir);
    aux=myMRIread([subjectDir '/' subjectName '/mri/norm.mgz'],0,tempdir);
    A.vol(aux.vol==0)=0;
else
    
    bfCorrAddVol=[subjectDir '/' subjectName '/mri/' analysisID '.thalamus.' BBregisterMode '.stripped.bfcorr.mgz'];
    if exist(bfCorrAddVol,'file')
        disp('Bias field corrected version of additional volume found; no need to recompute');
        A=myMRIread(bfCorrAddVol,0,tempdir);
    else
        disp('Bias field correction of additional scan');
        
        mri=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
        mri.vol(mri.vol==41)=2; % WM
        mri.vol(mri.vol==42)=3; % CT
        mri.vol(mri.vol==43)=4; % LV
        mri.vol(mri.vol==14)=4; % 3rd VENT
        mri.vol(mri.vol==15)=4; % 4th VENT
        mri.vol(mri.vol==60)=28; % DC
        mri.vol(mri.vol==46)=7; % CWM
        mri.vol(mri.vol==47)=8; % CCT
        mri.vol(mri.vol==49)=10; % LV
        mri.vol(mri.vol==50)=11; % CA
        mri.vol(mri.vol==51)=12; % PU
        mri.vol(mri.vol==52)=13; % PA
        mri.vol(mri.vol==53)=17; % HP
        mri.vol(mri.vol==54)=18; % AM
        mri.vol(mri.vol>=251 & mri.vol<=255)=2; % CC
        mri.vol(imdilate(grad3d(mri.vol)>0,createSphericalStrel(1)))=0;
        llist=[2 3 4 7 8 10 11 12 13 17 18 28];
        M=zeros(size(mri.vol))>1;
        for l=1:length(llist)
            M=M | mri.vol==llist(l);
        end
        mri.vol(~M)=0;
        % we work in log space
        mriAdd=myMRIread(registeredAddVol,0,tempdir);
        biasFieldOrder=4;
        PSI=prepBiasFieldBase(size(mriAdd.vol),biasFieldOrder);
        MASK=mri.vol>0;
        X=log(1+mriAdd.vol(MASK));
        L=mri.vol(MASK);
        PSIv=zeros([numel(X) size(PSI,4)]);
        for i=1:size(PSI,4)
            aux=PSI(:,:,:,i);
            PSIv(:,i)=aux(MASK);
        end
        
        means=zeros([1,length(llist)]);
        vars=zeros([1,length(llist)]);
        for l=1:length(llist)
            M=L==llist(l);
            aux=X(M);
            means(l)=mean(aux);
            vars(l)=var(aux,1);
        end
        ready=0;
        its=0;
        while ready==0
            its=its+1;
            WA=PSIv;
            XX=X;
            for l=1:length(llist)
                M=L==llist(l);
                XX(M)=XX(M)-means(l);
                WA(M,:)=WA(M,:)/vars(l);
            end
            tmp1=WA'*PSIv;
            tmp2=WA'*XX;
            C=tmp1\tmp2;
            Xcorr=X;
            for i=1:length(C)
                Xcorr=Xcorr-C(i)*PSIv(:,i);
            end
            means_old=means;
            vars_old=vars;
            for l=1:length(llist)
                M=L==llist(l);
                aux=Xcorr(M);
                means(l)=mean(aux);
                vars(l)=var(aux,1);
            end
            
            maxdif=max([max(abs(means-means_old)) sqrt(max(abs(vars-vars_old)))]);

            disp(['Bias field correction iteration ' num2str(its) ', maximal difference in parameters ' num2str(maxdif)]);
            
            if maxdif<1e-6 || its>10
                ready=1;
            end
        end
        AddCorr=log(1+mriAdd.vol);
        for i=1:length(C)
            AddCorr=AddCorr-C(i)*PSI(:,:,:,i);
        end
        AddCorr=exp(AddCorr)-1; % back to natural
        AddCorr(mriAdd.vol==0)=0;
        
        aux=mri;
        aux.vol=AddCorr;
        aux2=myMRIread([subjectDir '/' subjectName '/mri/norm.mgz'],0,tempdir);
        aux.vol(aux2.vol==0)=0;
        
        myMRIwrite(aux,bfCorrAddVol,'float',tempdir);
        
        A=aux;
        
    end
end
NORM=A; 
L=ASEG_TH_DE;
[~,cropping]=cropLabelVol(L.vol==THlabelLeft | L.vol==THlabelRight,round(margin/mean(A.volres)));
Lcrop=applyCropping(L.vol,cropping);
Icrop=applyCropping(A.vol,cropping);
offsetVox=[cropping(2)-1;cropping(1)-1;cropping(3)-1;0];
RAScorner=L.vox2ras0*offsetVox;
vox2ras0New=[ L.vox2ras0(:,1:3)  L.vox2ras0(:,4)+RAScorner];
I=L;
I.vol=Icrop;
I.vox2ras0=vox2ras0New;
myMRIwrite(I,'tempT1.mgz','float',tempdir);
cmd=[FSpath '/mri_convert tempT1.mgz '  imageFileName ' -odt float -rt  cubic -vs ' num2str(resolution) ' '  num2str(resolution)  ' '  num2str(resolution)];
system([cmd ' >/dev/null']);
delete tempT1.mgz
system([FSpath '/mri_binarize --i aseg_TH_DE.mgz --min 1.5 --dilate 2 --o asegModBinDilated.mgz >/dev/null']);
system([FSpath '/mri_convert asegModBinDilated.mgz asegModBinDilatedResampled.mgz -odt float -rt nearest -rl T1resampled.mgz >/dev/null']);
system([FSpath '/mri_mask -T 0.5 T1resampled.mgz asegModBinDilatedResampled.mgz T1resampled.mgz >/dev/null']);


% % Eugenio: let's try masking anything that is not close to the thalamus
% dilSize=round(5/mean(A.volres));
% system([FSpath '/mri_binarize --i asegMod.mgz --min 16.5 --max 18.5 --o hippoMask.mgz >/dev/null']);
% system([FSpath '/mri_binarize --i asegMod.mgz --min 16.5 --max 18.5 --o hippoMaskDilated5mm.mgz --dilate ' num2str(dilSize) ' >/dev/null']);
% system([FSpath '/mri_convert hippoMask.mgz hippoMaskResampled.mgz -rt interpolate -rl T1resampled.mgz  -odt float >/dev/null']);
% system([FSpath '/mri_binarize --i  hippoMaskResampled.mgz --min 0.5 --dilate ' num2str(round(3/resolution)) '  --o hippoMaskResampledDilated.mgz  >/dev/null']);
% system([FSpath '/mri_mask -T 0.5 T1resampled.mgz  hippoMaskResampledDilated.mgz T1resampled.mgz >/dev/null']);



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
% mask = imerode(single( sum( priors, 4 ) / 65535 ) > 0.99,createSphericalStrel(5));
mask = imerode(single( sumpriors / 65535 ) > 0.99,createSphericalStrel(5));
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

FreeSurferLabelGroups=[];
FreeSurferLabelGroups{end+1}={'Unknown'};

% FreeSurferLabelGroups{end+1}={'Left-Cerebral-White-Matter'};
FreeSurferLabelGroups{end+1}={'Left-Cerebral-White-Matter','Left-R','Right-R'};
WMind=2;

FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex'};
GMind=3;

FreeSurferLabelGroups{end+1}={'Left-Cerebellum-Cortex'};
FreeSurferLabelGroups{end+1}={'Left-Cerebellum-White-Matter'};
FreeSurferLabelGroups{end+1}={'Brain-Stem'};
FreeSurferLabelGroups{end+1}={'Left-Lateral-Ventricle'};
FreeSurferLabelGroups{end+1}={'Left-choroid-plexus'};
FreeSurferLabelGroups{end+1}={'Left-Putamen'};
FreeSurferLabelGroups{end+1}={'Left-Pallidum'};
FreeSurferLabelGroups{end+1}={'Left-Accumbens-area'};
FreeSurferLabelGroups{end+1}={'Left-Caudate'};
FreeSurferLabelGroups{end+1}={'Left-VentralDC','Right-VentralDC'};

% Note to myself: thalamus must be last component

% FUll model
FreeSurferLabelGroups{end+1}={'Left-L-Sg','Left-LGN','Left-MGN','Left-PuI','Left-PuM','Left-H','Left-PuL',...
    'Left-VPI','Left-PuA','Left-MV(Re)','Left-Pf','Left-CM','Left-LP','Left-VLa','Left-VPL','Left-VLp',...
    'Left-MDm','Left-VM','Left-CeM','Left-MDl','Left-Pc','Left-MDv','Left-Pv','Left-CL','Left-VA','Left-VPM',...
    'Left-AV','Left-VAmc','Left-Pt','Left-AD','Left-LD','Right-L-Sg','Right-LGN','Right-MGN','Right-PuI','Right-PuM','Right-H','Right-PuL',...
    'Right-VPI','Right-PuA','Right-MV(Re)','Right-Pf','Right-CM','Right-LP','Right-VLa','Right-VPL','Right-VLp',...
    'Right-MDm','Right-VM','Right-CeM','Right-MDl','Right-Pc','Right-MDv','Right-Pv','Right-CL','Right-VA','Right-VPM',...
    'Right-AV','Right-VAmc','Right-Pt','Right-AD','Right-LD'};


sameGaussianParameters=[];
for g=1:length(FreeSurferLabelGroups)
    sameGaussianParameters{end+1} = [];
    for FreeSurferLabel =  FreeSurferLabelGroups{g}
        sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
    end
    if isempty(sameGaussianParameters{end})
        sameGaussianParameters=sameGaussianParameters(1:end-1);
    end
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
        subplot( 3, 4, reducedLabel )
        showImage( priors( :, :, :, reducedLabel ) )
        title(['Prior for "reduced" class ' num2str(reducedLabel)]);
    end
    priorsBefore=priors;
    clear priors
end


% Compute hyperparameters for estimation of Gaussian parameters
disp('Computing hyperparameters for estimation of Gaussian parameters')

if isempty(additionalVol)
    DATA=myMRIread([subjectDir '/' subjectName '/mri/norm.mgz'],0,tempdir);
elseif doBFcorrection==0
    DATA=myMRIread(registeredAddVol,0,tempdir);
    aux=myMRIread([subjectDir '/' subjectName '/mri/norm.mgz'],0,tempdir);
    DATA.vol(aux.vol==0)=0;
else
    DATA=myMRIread(bfCorrAddVol,0,tempdir);
end

nHyper=zeros(length(sameGaussianParameters),1);
meanHyper=zeros(length(sameGaussianParameters),1);
BGindex=0;
for g=1:length(sameGaussianParameters)
    labels=sameGaussianParameters{g};
    if any(labels==0)
        BGindex=g;
    end
    if any(labels>8225)  % thalamus
        listMask=[10 49];
    elseif any(labels==28) % VDE
        listMask=[28 60];
    elseif any(labels==0) % background
        listMask=[1];
    else
        listMask=labels;
    end
    
    if length(listMask)>0
        MASK=zeros(size(DATA.vol));
        for l=1:length(listMask)
            MASK=MASK | ASEGmod.vol==listMask(l);
        end
        % MASK=imerode(MASK,createSphericalStrel(1));
        MASK=imerode(MASK,createSphericalStrel(round(1/mean(DATA.volres))));
        data=DATA.vol(MASK & DATA.vol>0);
        meanHyper(g)=median(data);
        if any(labels==28) % special case... VDE is kind of bimodal in FS
            nHyper(g)=10;
        else
            nHyper(g)=10+length(data)*prod(DATA.volres)/resolution^3;
        end
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
if useTwoComponents>0  % iterate a bit more ...
    meshSmoothingSigmas = [ 1.5 1.125 .75 0 ]';
    imageSmoothingSigmas = [0 0 0 0]';
    maxItNos=[7 5 5 3];  % each iteration has 20 deformation steps
else
    meshSmoothingSigmas = [ 1.5 .75 0 ]';
    imageSmoothingSigmas = [0 0 0]';
    maxItNos=[7 5 3];  % each iteration has 20 deformation steps
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
    
    % I added this: if we want to use two Gaussian compoments in the
    % mixture for the thalamus, we do it after the first resolution level
    if useTwoComponents>0 && multiResolutionLevel==2
        
        FreeSurferLabelGroups=[];
        FreeSurferLabelGroups{end+1}={'Unknown'};
        FreeSurferLabelGroups{end+1}={'Left-Cerebral-White-Matter','Left-R','Right-R'};
        FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex'};
        FreeSurferLabelGroups{end+1}={'Left-Cerebellum-Cortex'};
        FreeSurferLabelGroups{end+1}={'Left-Cerebellum-White-Matter'};
        FreeSurferLabelGroups{end+1}={'Brain-Stem'};
        FreeSurferLabelGroups{end+1}={'Left-Lateral-Ventricle'};
        FreeSurferLabelGroups{end+1}={'Left-choroid-plexus'};
        FreeSurferLabelGroups{end+1}={'Left-Putamen'};
        FreeSurferLabelGroups{end+1}={'Left-Pallidum'};
        FreeSurferLabelGroups{end+1}={'Left-Accumbens-area'};
        FreeSurferLabelGroups{end+1}={'Left-Caudate'};
        FreeSurferLabelGroups{end+1}={'Left-VentralDC','Right-VentralDC'};
        
        % Note to myself: medial nuclei MUST be last, and lateral, second to last
   
        
        % For full model
        FreeSurferLabelGroups{end+1}={'Left-L-Sg','Left-LGN','Left-MGN','Left-H',...
            'Left-VPI','Left-MV(Re)','Left-Pf','Left-CM','Left-LP','Left-VLa','Left-VPL','Left-VLp',...
            'Left-VM','Left-CeM','Left-Pc','Left-MDv','Left-Pv','Left-CL','Left-VA','Left-VPM',...
            'Left-AV','Left-VAmc','Left-Pt','Left-AD','Left-LD','Right-L-Sg','Right-LGN','Right-MGN','Right-H',...
            'Right-VPI','Right-MV(Re)','Right-Pf','Right-CM','Right-LP','Right-VLa','Right-VPL','Right-VLp',...
            'Right-VM','Right-CeM','Right-Pc','Right-MDv','Right-Pv','Right-CL','Right-VA','Right-VPM',...
            'Right-AV','Right-VAmc','Right-Pt','Right-AD','Right-LD'};
        FreeSurferLabelGroups{end+1}={'Left-PuA','Left-PuI','Left-PuL','Left-PuM','Left-MDl','Left-MDm',...
            'Right-PuA','Right-PuI','Right-PuL','Right-PuM','Right-MDl','Right-MDm'};
        
        
        
        
        sameGaussianParameters=[];
        for g=1:length(FreeSurferLabelGroups)
            sameGaussianParameters{end+1} = [];
            for FreeSurferLabel =  FreeSurferLabelGroups{g}
                sameGaussianParameters{end} = [ sameGaussianParameters{end} FreeSurferLabels( find( strcmp( FreeSurferLabel, cellstr( names ) ) ) ) ];
            end
            if isempty(sameGaussianParameters{end})
                sameGaussianParameters=sameGaussianParameters(1:end-1);
            end
        end
        numberOfClasses=length(sameGaussianParameters);
        
        % Compute the "reduced" alphas - those referring to the "super"-structures
        [ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
        if ( max( abs( sum( reducedAlphas, 2 ) - 1 ) ) > 1e-5 ) % Make sure these vectors really sum to 1
            error( 'The vector of prior probabilities in the mesh nodes must always sum to one over all classes' )
        end
        % Set the reduced alphas to be the alphas of the mesh
        kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
        
        % Modify hyperparameters
        ThInt=meanHyper(end);
        
        if isempty(additionalVol)
            % lateral, brighter
            nHyper(end)=25;
            meanHyper(end)=ThInt+5;
            % medial, darker
            nHyper(end+1)=25;
            meanHyper(end+1)=ThInt-5;
            
        else
            
            nHyper(end)=25;
            nHyper(end+1)=25;
            
            % lateral, more WM-ish (e.g., darker, in FGATIR)
            meanHyper(end)=ThInt*(0.95+0.1*(meanHyper(WMind)>=meanHyper(GMind)));
            % medial, more GM-ish (e.g., brighter, in FGATIR)
            meanHyper(end+1)=ThInt*(0.95+0.1*(meanHyper(WMind)<meanHyper(GMind)));
            
        end
        
    end  % end of IF changing Gaussian components of thalamus
    
    
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
    
    
    % Iterations
    % Eugenio November 2017: we now do 30 instead of 20, since it's a bit faster
    maximumNumberOfIterations = maxItNos(multiResolutionLevel);  % Maximum number of iterations (includes one imaging model parameter estimation and
    positionUpdatingMaximumNumberOfIterations = 30;
    
    if FAST>0, maximumNumberOfIterations = 3; end  % in case we just wanna cruise throught it :-)
    
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
    
    disp(['Resolution level: ' num2str(multiResolutionLevel) ' of ' num2str(numberOfMultiResolutionLevels)]);
    
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
        %         priors = priors( maskIndices, : );  % Ignore everything that's has zero intensity
        
        data = data( maskIndices );
        
        % Start EM iterations. Initialize the parameters if this is the
        % first time ever you run this
        EPS=1e-2;
        posteriors = double( priors ) / 65535;
        if ( ( iterationNumber == 1 ) &&  ( multiResolutionLevel == 1 || (useTwoComponents>0  && multiResolutionLevel == 2) ) )
            
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
            
            minLogLikelihood =  minLogLikelihood - sum( log( normalizer ) ); % This is what we're optimizing with EM
            if THALAMUS_VERBOSE
                minLogLikelihood
            end
            
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
            for classNumber = 1 : numberOfClasses
                posterior = posteriors( :, classNumber );
                EPS=1e-2;
                
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
        if THALAMUS_VERBOSE
            means'
            (variances').^2
        end
        
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
        
        % Eugenio November 2017: GEMS2
        if ( exist( 'optimizer', 'var' ) == 1 )
            % The optimizer is very memory hungry when run in multithreaded mode.
            % Let's clear any old ones we may have lying around
            try
                % Eugenio November 2017: GEMS2
                kvlClear( optimizer );
                kvlClear( calculator );
            catch ME
            end
        end
        
        
        haveMoved = false; % Keep track if we've ever moved or not
        
         % (note that it uses variances instead of precisions)
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
            if THALAMUS_VERBOSE>0 || mod(positionUpdatingIterationNumber,10)==1
                disp(['Resolution level ' num2str(multiResolutionLevel) ' iteration ' num2str(iterationNumber) ' deformation iterations ' num2str(positionUpdatingIterationNumber)]);
            end
            % Eugenio May2018
            maximalDeformation=0;
            try
                tic
                % Eugenio November 2017: GEMS2
                [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer );
                elapsedTime = toc;
                if THALAMUS_VERBOSE
                    disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
                    minLogLikelihoodTimesPrior
                end
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

disp(['Fitting mesh to image took ' num2str(etime(clock,time_ref_optimization)) ' seconds']);

% Restore original image buffer
kvlSetImageBuffer(image,imageBufferOrig);

% OK, now that all the parameters have been estimated, segment the image with all the original
% labels instead of the reduced "super"-structure labels we created.

% Clear some memory
% Eugenio November 2017: GEMS2
try
    kvlClear( optimizer )
    kvlClear( calculator )
end

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

% 
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
    
    priors = kvlRasterizeAtlasMesh( mesh, imageSize );
    priors = reshape( priors, [ prod( imageSize ) numberOfClasses ] );  % Easier to work with vector notation in the computations
    priors = priors( maskIndices, : );
    
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
system([FSpath '/mri_convert asegMod.mgz asmr1.mgz -rl T1resampled.mgz -rt nearest -odt float >/dev/null']);
[ asmr, tasmr ] = kvlReadCroppedImage( 'asmr1.mgz', boundingFileName );
asmrB=kvlGetImageBuffer(asmr);
kvlWriteImage( asmr, 'asmr2.mgz' );
tmp1=myMRIread('asmr1.mgz',0,tempdir);
tmp2=myMRIread('asmr2.mgz',0,tempdir);
tmp3=tmp1; tmp3.vol=tmp2.vol;
myMRIwrite(tmp3,'asmr3.mgz','float',tempdir);
[I,J,K]=ind2sub(size(tmp1.vol),find(tmp1.vol==10));
Ic1=mean(I); Jc1=mean(J); Kc1=mean(K);
[I,J,K]=ind2sub(size(tmp2.vol),find(tmp2.vol==10));
Ic2=mean(I); Jc2=mean(J); Kc2=mean(K);
shift=round([Ic2-Ic1,Jc2-Jc1,Kc2-Kc1]);
shiftNeg=-shift; shiftNeg(shiftNeg<0)=0;
shiftPos=shift; shiftPos(shiftPos<0)=0;

% reorient image.mgz
aux=zeros(size(tmp2.vol)+shiftNeg);
aux2=myMRIread('image.mgz',0,tempdir);
aux(1+shiftNeg(1):shiftNeg(1)+size(aux2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(aux2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(aux2.vol,3))=aux2.vol;
aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
tmp3.vol=aux;
myMRIwrite(tmp3,'image.mgz','float',tempdir);




% Memory efficient version for posteriors and MAP segmentation
strInt={'L-Sg','LGN','MGN','PuI','PuM','H','PuL',...
    'VPI','PuA','MV(Re)','Pf','CM','LP','VLa','VPL','VLp',...
    'MDm','VM','CeM','MDl','Pc','MDv','Pv','CL','VA','VPM',...
    'AV','VAmc','Pt','AD','LD'}; % ,'R'  I now leave reticular out
fid=fopen([tempdir '/volumesThalamus.txt'],'w');
totVolL=0;
totVolR=0;
foundL=zeros([1,numberOfClasses]);
foundR=zeros([1,numberOfClasses]);
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
    
    foundL(i)=0;
    foundR(i)=0;
    
    name=names(i,:);
    name=lower(name(name~=' '));
    f=find(name=='-');
    if ~isempty(f) && strcmp(name(1:f(1)-1),'left')
        side='left';        
    elseif ~isempty(f) && strcmp(name(1:f(1)-1),'right')
        side='right';
    else
        side=[];
    end
    
    for j=1:length(strInt)
        if strcmp(name,lower(['left-' strInt{j}]))>0
            foundL(i)=j;
        elseif strcmp(name,lower(['right-' strInt{j}]))>0
            foundR(i)=j;
        end
    end
    
    if foundL(i)>0 || foundR(i)>0
        vol=resolution^3*(sum(double(post(:))/65535));
        if foundL(i)>0
            str=strInt{foundL(i)};
            fprintf(fid,'Left-%s %f\n',str,vol);
            if isempty(strfind(lower(names(i,:)),'hippocampal-fissure'))  % don't count the fissure towards the total volume
                totVolL=totVolL+vol;
            end
        else
            str=strInt{foundR(i)};
            fprintf(fid,'Right-%s %f\n',str,vol);
            totVolR=totVolR+vol;
        end
        
        if WRITE_POSTERIORS>0 
            if useTwoComponents>0
                suffix2=[suffix '.' analysisID]; % default
            else
                suffix2=[suffix '.OneGaussComp.' analysisID];
            end
            kk1=double(post)/65535;
            kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
            aux=zeros(size(tmp2.vol)+shiftNeg);
            aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
            aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
            tmp3.vol=aux;
            myMRIwrite(tmp3,['posterior_' side '_' strtrim(lower(names(i,:))) '.' suffix2 '.mgz'],'float',tempdir);
        end
        
    end
end
fprintf(fid,'Left-Whole_thalamus %f\n',totVolL);
fprintf(fid,'Right-Whole_thalamus %f\n',totVolR);
fclose(fid);




% MAP estimates
% [~,inds]=max(posteriorsFull,[],4);
inds=L;
kk1=FreeSurferLabels(inds);
kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
aux=zeros(size(tmp2.vol)+shiftNeg);
aux(1+shiftNeg(1):shiftNeg(1)+size(kk1,2),1+shiftNeg(2):shiftNeg(2)+size(kk1,1),1+shiftNeg(3):shiftNeg(3)+size(kk1,3))=permute(kk1,[2 1 3]);
aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
tmp3.vol=aux;
myMRIwrite(tmp3,'discreteLabels_all.mgz','float',tempdir);
tmp3.vol(tmp3.vol<100 & tmp3.vol~=10 & tmp3.vol~=49 )=0; 

% kill reticular
LEFT_RET=FreeSurferLabels(find( strcmp( 'Left-R', cellstr( names ) ) ));
RIGHT_RET=FreeSurferLabels(find( strcmp( 'Right-R', cellstr( names ) ) ));
tmp3.vol(tmp3.vol==LEFT_RET)=0;
tmp3.vol(tmp3.vol==RIGHT_RET)=0;

% Get only connectec components
% (Sometimes the two thalami are not connected...)
m1=getLargestCC(tmp3.vol<8200 & (tmp3.vol>100 | tmp3.vol==THlabelLeft) );
m2=getLargestCC(tmp3.vol>8200 | tmp3.vol==THlabelRight)  ;
tmp3Mask =m1 | m2;
tmp3.vol(~tmp3Mask)=0;
myMRIwrite(tmp3,'discreteLabels.mgz','float',tempdir);



% Eugenio July 2017
% I disabled this for now ...
if  MRFconstant>0

    % Eugenio July 2017
    error('MRF smoothing disabled for now');
    
%     EPS=1e-12;
%     [~,inds]=max(posteriorsFull,[],4);
%     tmp=FreeSurferLabels(inds);
%     kk=zeros(size(tmp)); kk(maskIndices)=1; tmp=tmp.*kk;
%     tmp(tmp<100)=0; tmp(tmp>300)=0;
%     [~,cropping]=cropLabelVol(tmp);
%     Ct=zeros([cropping(4)-cropping(1)+1,cropping(5)-cropping(2)+1,cropping(6)-cropping(3)+1,numberOfClasses]);
%     for c=1:numberOfClasses
%         Ct(:,:,:,c)=-log(EPS+double(posteriorsFull(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6),c))/65535);
%     end
%     factor=-256/log(EPS);
%     Ct=int32(round(Ct*factor));
%     unaryTermWeight=int32(round(MRFconstant*factor));
%     
%     siz=[size(Ct,1) size(Ct,2) size(Ct,3)];
%     h = GCO_Create(prod(siz),numberOfClasses);
%     DC = zeros([numberOfClasses,prod(siz)],'int32');
%     for c=1:numberOfClasses
%         aux=Ct(:,:,:,c);
%         DC(c,:)=aux(:);
%     end
%     GCO_SetDataCost(h,DC);
%     aux=int32(double(unaryTermWeight)*(ones(numberOfClasses)-eye(numberOfClasses)));
%     GCO_SetSmoothCost(h,aux);
%     
%     row=zeros([prod(siz)*3,1]);
%     col=zeros([prod(siz)*3,1]);
%     t=1;
%     
%     Ifrom=1:siz(1)-1;
%     Ito=2:siz(1);
%     inc=length(Ito);
%     for j=1:siz(2)
%         J=j*ones(size(Ifrom));
%         for k=1:siz(3)
%             K=k*ones(size(Ifrom));
%             row(t:t+inc-1)=sub2ind(siz,Ifrom,J,K);
%             col(t:t+inc-1)=sub2ind(siz,Ito,J,K);
%             t=t+inc;
%         end
%     end
%     
%     Jfrom=1:siz(2)-1;
%     Jto=2:siz(2);
%     inc=length(Jto);
%     for i=1:siz(1)
%         I=i*ones(size(Jfrom));
%         for k=1:siz(3)
%             K=k*ones(size(Jfrom));
%             row(t:t+inc-1)=sub2ind(siz,I,Jfrom,K);
%             col(t:t+inc-1)=sub2ind(siz,I,Jto,K);
%             t=t+inc;
%         end
%     end
%     
%     Kfrom=1:siz(3)-1;
%     Kto=2:siz(3);
%     inc=length(Kto);
%     for i=1:siz(1)
%         I=i*ones(size(Kfrom));
%         for j=1:siz(2)
%             J=j*ones(size(Kfrom));
%             row(t:t+inc-1)=sub2ind(siz,I,J,Kfrom);
%             col(t:t+inc-1)=sub2ind(siz,I,J,Kto);
%             t=t+inc;
%         end
%     end
%     
%     row=row(1:t-1);
%     col=col(1:t-1);
%     
%     NEIGH=sparse(row,col,ones(size(row)),prod(siz),prod(siz));
%     GCO_SetNeighbors(h,NEIGH);
%     
%     
%     GCO_Expansion(h);      % Compute optimal labeling via alpha-expansion
%     ind=reshape(GCO_GetLabeling(h),siz);
%     
%     SEG=FreeSurferLabels(ind);
%     SEG(SEG<100)=0;
%     SEG(SEG>300)=0;
%     
%     data=zeros(size(inds));
%     data(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6))=SEG;
%     aux=zeros(size(tmp2.vol)+shiftNeg);
%     aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(data,[2 1 3]);
%     aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
%     tmp3.vol=aux;
%     tmp3Mask=getLargestCC(tmp3.vol>0);
%     tmp3.vol(~tmp3Mask)=0;
%     myMRIwrite(tmp3,'discreteLabels_MRF.mgz','float',tempdir);
    
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
% Eugenio November 2017: smarter/smoother resampling
if SMOOTH_LABEL_RESAMPLE>0
    refFile=[subjectDir '/' subjectName '/mri/norm.mgz'];
    applyLTAsmoothLabels('discreteLabels.mgz',[],'discreteLabelsResampledT1.mgz',refFile,0,FSpath,tempdir);
else
    system([FSpath '/mri_convert  discreteLabels.mgz  discreteLabelsResampledT1.mgz -rt nearest -odt float ' ...
    ' -rl ' subjectDir '/' subjectName '/mri/norm.mgz >/dev/null']);
end


% Move to MRI directory (modify suffix to reflect number of components)
if useTwoComponents>0
    suffix2=[suffix '.' analysisID]; % default
else
    suffix2=[suffix '.OneGaussComp.' analysisID];
end

% Move to MRI directory
system(['mv discreteLabels.mgz ' subjectDir '/' subjectName '/mri/ThalamicNuclei.' suffix2 '.mgz']);
system(['mv discreteLabelsResampledT1.mgz ' subjectDir '/' subjectName '/mri/ThalamicNuclei.' suffix2 '.FSvoxelSpace.mgz']);
system(['mv volumesThalamus.txt ' subjectDir '/' subjectName '/mri/ThalamicNuclei.' suffix2 '.volumes.txt']);


if WRITE_POSTERIORS>0
    system(['mv posterior*.mgz '  subjectDir '/' subjectName '/mri/']);
end

% Eugenio November 2017
if WRITE_MESHES>0
    system(['mv warpedMesh.txt.gz   ' subjectDir '/' subjectName '/mri/ThalamicNuclei.' suffix2 '.txt.gz']);
    system(['mv warpedMeshNoAffine.txt.gz   ' subjectDir '/' subjectName '/mri/ThalamicNucleiMeshAtlasSpace.' suffix2 '.txt.gz']);
    system(['mv image.mgz   ' subjectDir '/' subjectName '/mri/ThalamicImageForMesh.' suffix2 '.mgz']);    
    fid=fopen([subjectDir '/' subjectName '/mri/ThalamicAffineTransformMesh.' suffix2 '.txt'],'w');
    fprintf(fid,'%f %f %f %f \n',transformMatrix(1,1),transformMatrix(1,2),transformMatrix(1,3),transformMatrix(1,4));
    fprintf(fid,'%f %f %f %f \n',transformMatrix(2,1),transformMatrix(2,2),transformMatrix(2,3),transformMatrix(2,4));  
    fprintf(fid,'%f %f %f %f \n',transformMatrix(3,1),transformMatrix(3,2),transformMatrix(3,3),transformMatrix(3,4));
    fprintf(fid,'%f %f %f %f \n',0,0,0,1);
    fclose(fid);
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
