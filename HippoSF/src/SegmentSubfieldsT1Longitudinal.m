% segments the subfields (of multiple time points) in a longitudinal fashion.
% It uses ASEG the base as initialization for the subject-specific atlas.
% It is based on  segmentSubjectT1_autoEstimateAlveusM
%
% This version uses the same initialization for all time points (computed
% from the base)
%
% SegmentSubfieldsT1Longitudinal(subjectDir,resolution,atlasMeshFileName,atlasDumpFileName,compressionLUTfileName,K,side,optimizerType,suffix,MRFconstant,subjectBase,subjectTP1,subjectTP2,...)
%
% - subjectName: FreeSurfer subject name
% - subjectDir: FreeSurfer subject directory
% - resolution: voxel size at which we want to work (in mm).
% - atlasMeshFileName: the atlas to segment the data
% - atlasDumpFileName: corresponding imageDump.mgz (name *must* be imageDump.mgz)
% - compressionLUTfileName: corresponding compressionLUT.txt
% - K: stiffness of the mesh in the segmentation.
% - side: 'left' or 'right'
% - optimizerType: 'FixedStepGradientDescent','GradientDescent','ConjugateGradient','L-BFGS'
% - suffix: for output directory, e.g. 'T1based_GGAWLnoSimil'
% - FSpath: path to FreeSurfer executables
% - MRFconstant (optional): make it >0 for MRF cleanup (5 is reasonable, larger is smoother)
%           It does NOT affect volumes, which are computed from soft posteriors anyway
% - subjectBase: subject name of base
% - subjectTP1: subjectname of time point 1
% - subjectTP2: subjectname of time point 2
% - ...
% - subjectTPn: subjectname of time point n


function SegmentSubfieldsT1Longitudinal(varargin)


if nargin<14
    error('Requires at least 14 arguments (one time point)')
end

subjectDir=varargin{1};
resolution=varargin{2};
atlasMeshFileName=varargin{3};
atlasDumpFileName=varargin{4};
compressionLUTfileName=varargin{5};
Katl=varargin{6};
Ktp=varargin{7};
side=varargin{8};
optimizerType=varargin{9};
suffix=varargin{10};
FSpath=varargin{11};
MRFconstant=varargin{12};
subjectBase=varargin{13};
nTP=nargin-13;
subjectTPs=cell([1,nTP]);
for t=1:nTP
    subjectTPs{t}=varargin{13+t};
end


% In case we compiled it...
if isdeployed
    Katl=str2double(Katl);
    Ktp=str2double(Ktp);
    resolution=str2double(resolution);
    MRFconstant=str2double(MRFconstant);
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

% Sanity check
if strcmp(side,'left')==0 && strcmp(side,'right')==0
    error('Side must be ''left'' or ''right''');
elseif optimizerType(1)~='F' && optimizerType(1)~='G' && optimizerType(1)~='C' && optimizerType(1)~='L'
    error('Optimizer type must be ''FixedStepGradientDescent'',''GradientDescent'',''ConjugateGradient'',''L-BFGS''');
elseif exist([subjectDir '/' subjectBase],'dir')==0
    error('Subject directory for base does not exist');
elseif ~isdeployed && (~isnumeric(resolution))
    error('Resolution must be numeric');
elseif exist(atlasMeshFileName,'file')==0
    error('Provided atlas mesh file does not exist');
elseif exist(atlasDumpFileName,'file')==0
    error('Provided imageDump.mgz does not exist');
elseif exist(compressionLUTfileName,'file')==0
    error('Provided LUT does not exist');
elseif ~isdeployed && (~isnumeric(Katl))
    error('Katl must be numeric');
elseif ~isdeployed && (~isnumeric(Ktp))
    error('Ktp must be numeric');
elseif ~isdeployed && (~isnumeric(MRFconstant))
    error('MRFconstant must be numeric');
end
for t=1:nTP
    if exist([subjectDir '/' subjectTPs{t}],'dir')==0
        error(['Subject directory for time point ' num2str(t) ' does not exist']);
    end
end


% Constants
HippoLabelLeft=17;
HippoLabelRight=53;
MAX_GLOBAL_ITS=5;
DEBUG=0;
FAST=0; % set it to one to optimize just a bit (go through code fast)
WRITE_POSTERIORS=0;
% Eugenio November 2017: added option to write meshes and smoother resampling
WRITE_MESHES=0;
SMOOTH_LABEL_RESAMPLE=0;
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

% March 2021: fix to accommodate 'fs_run_from_mcr'
FSpath = [FSpath '/fs_run_from_mcr ' FSpath '/'];

% File name for temporary directory
tempdir=[subjectDir '/' subjectBase '/tmp/hippoSF_T1_' suffix '_' side '/'];
aux=getenv('USE_SCRATCH');
if ~isempty(aux)
    if str2double(aux)>0
        if exist('/scratch/','dir')>0
            tempdir=['/scratch/' subjectBase '_hippoSF_T1_' suffix '_' side '/'];
        end
    end
end

% Temporary directory: create if necessary
if exist(tempdir,'dir')==0
    mkdir(tempdir);
end
cd(tempdir);

% Produce error if resolutions are not all the same
aux=myMRIread([subjectDir '/' subjectBase '/mri/aseg.mgz'],0,tempdir);
volresPrev=aux.volres;
for t=1:nTP
    if any (volresPrev ~= aux.volres)
        error('Time points/base do not all have the same resolution');
    else
        volresPrev = aux.volres;
    end
end
highres=0; if mean(aux.volres)<0.99, highres=1; end


%%%%%%%%%%%%%%%%%%%%%%%%
% Stopwatch: let's go! %
%%%%%%%%%%%%%%%%%%%%%%%%
time_start=clock;

% Clean up KVL memory space
kvlClear;


% Next: register image dump to automated segmentation of base.
% This will provide the global affine coordinate frame for the mesh
disp('Registering imageDump.mgz to hippocampal mask from ASEG of Base Point')

% Grab atlas (manually placed in FS atlas coordinate space!)
system(['cp ' atlasDumpFileName ' ./imageDump.mgz']);

% flip LR if right side - we only rotate along LR axis not to bias left vs right hippo segmentation
if strcmp(side,'right')>0
    aux=myMRIread('imageDump.mgz',0,tempdir);
    aux.vox2ras0(1,:)=-aux.vox2ras0(1,:);
    aux.vox2ras1(1,:)=-aux.vox2ras1(1,:);
    aux.vox2ras(1,:)=-aux.vox2ras(1,:);
    myMRIwrite(aux,'imageDump.mgz','float',tempdir);
end


% Target is masked aseg (hippocampus and amygdala)
targetRegFileName=[tempdir '/hippoAmygBinaryMask.mgz'];
targetRegFileNameCropped=[tempdir '/hippoAmygBinaryMask_autoCropped.mgz'];
targetRegFileNameCroppedOpened=[tempdir '/hippoAmygBinaryMask_autoCropped_opened.mgz'];
ASEG=myMRIread([subjectDir '/' subjectBase '/mri/aseg.mgz'],0,tempdir);
TARGETREG=ASEG;
if strcmp(side,'left')>0
    TARGETREG.vol=255*double(ASEG.vol==HippoLabelLeft | ASEG.vol==HippoLabelLeft+1);
else
    TARGETREG.vol=255*double(ASEG.vol==HippoLabelRight | ASEG.vol==HippoLabelRight+1);
end
myMRIwrite(TARGETREG,targetRegFileName,'float',tempdir);
if highres==1,
    system([FSpath '/mri_convert ' targetRegFileName ' aux.mgz -odt float -vs 1 1 1 -rt nearest >/dev/null']);
    system(['mv aux.mgz ' targetRegFileName ' >/dev/null']);
end


% cmd=[FSpath '/kvlAutoCrop ' targetRegFileName ' 6'];
% system([cmd ' >/dev/null']);

aux=myMRIread(targetRegFileName,0,tempdir);
[aux.vol,cropping]=cropLabelVol(aux.vol,6);
shift=aux.vox2ras0(1:3,1:3)*[cropping(2)-1; cropping(1)-1; cropping(3)-1];
aux.vox2ras0(1:3,4)=aux.vox2ras0(1:3,4)+shift;
aux.vox2ras1(1:3,4)=aux.vox2ras1(1:3,4)+shift;
aux.vox2ras(1:3,4)=aux.vox2ras(1:3,4)+shift;
aux.tkrvox2ras=[];
myMRIwrite(aux,targetRegFileNameCropped,'float',tempdir);



% Initial affine alignment based just on (opened version of) hippocampus
aux=myMRIread(targetRegFileNameCropped,0,tempdir);
strel=createSphericalStrel(1);
aux.vol=255*double(imdilate(imerode(aux.vol>0,strel),strel));
myMRIwrite(aux,targetRegFileNameCroppedOpened,'float',tempdir);

% cmd=[FSpath '/kvlRegister imageDump.mgz ' targetRegFileNameCroppedOpened ' 3 2'];
% system(cmd);
% system('mv imageDump_coregistered.mgz imageDump.mgz' );

cmd=[FSpath '/mri_robust_register --mov imageDump.mgz  --dst ' targetRegFileNameCroppedOpened ...
    ' -lta trash.lta --mapmovhdr imageDump_coregistered.mgz  --sat 50'];
system(cmd);
system('cp imageDump_coregistered.mgz imageDump_rigidly_coregistered.mgz' );
system('mv imageDump_coregistered.mgz imageDump.mgz' );

cmd=[FSpath '/mri_robust_register --mov imageDump.mgz  --dst ' targetRegFileNameCroppedOpened ...
    ' -lta rigidToAffineBase.lta --mapmovhdr imageDump_coregistered.mgz --affine --sat 50'];
system(cmd);
system('mv imageDump_coregistered.mgz imageDump.mgz' );


% Also, we compute an affine aligment of the time points to the base based
% solely of the hippocampal segmentation from aseg
% We need to remember the determinant of the transform - we'll need to
% divide the final volume estimates by it
%
disp('Registering hippocampal masks of time points to mask of base')
volumeFactors=zeros([1,nTP]);

fixed=targetRegFileNameCropped;
for t=1:nTP
    aux=myMRIread([subjectDir '/' subjectTPs{t} '/mri/aseg.mgz'],0,tempdir);
    if strcmp(side,'left')>0
        aux.vol=255*double(aux.vol==HippoLabelLeft | aux.vol==HippoLabelLeft+1);
    else
        aux.vol=255*double(aux.vol==HippoLabelRight | aux.vol==HippoLabelRight+1);
    end
    myMRIwrite(aux,['hippoAmygBinaryMask_tp_' num2str(t) '.mgz'],'float',tempdir);
    
    if highres==1,
        system([FSpath '/mri_convert hippoAmygBinaryMask_tp_' num2str(t) '.mgz aux.mgz -odt float -vs 1 1 1 -rt nearest >/dev/null']);
        system(['mv aux.mgz hippoAmygBinaryMask_tp_' num2str(t) '.mgz >/dev/null']);
    end
    
    %     cmd=[FSpath '/kvlAutoCrop hippoAmygBinaryMask_tp_' num2str(t) '.mgz 6'];
    %     system([cmd ' >/dev/null']);
    
    aux=myMRIread(['hippoAmygBinaryMask_tp_' num2str(t) '.mgz'],0,tempdir);
    [aux.vol,cropping]=cropLabelVol(aux.vol,6);
    shift=aux.vox2ras0(1:3,1:3)*[cropping(2)-1; cropping(1)-1; cropping(3)-1];
    aux.vox2ras0(1:3,4)=aux.vox2ras0(1:3,4)+shift;
    aux.vox2ras1(1:3,4)=aux.vox2ras1(1:3,4)+shift;
    aux.vox2ras(1:3,4)=aux.vox2ras(1:3,4)+shift;
    aux.tkrvox2ras=[];
    myMRIwrite(aux,['hippoAmygBinaryMask_tp_' num2str(t) '_autoCropped.mgz'],'float',tempdir);
    
    
    
    % Initial affine alignment based just on  hippocampus
    moving=['hippoAmygBinaryMask_tp_' num2str(t) '_autoCropped.mgz'];
    
    %     AFFINE
    cmd=[FSpath '/mri_robust_register --mov ' moving '  --dst ' fixed ...
        ' -lta tp_' num2str(t) '_to_base.lta  --mapmovhdr kk2.mgz --affine   --sat 50'];
    system(cmd);
    
    % %     % RIGID
    %     cmd=[FSpath '/mri_robust_register --mov ' moving '  --dst ' fixed ...
    %         ' -lta tp_' num2str(t) '_to_base.lta  --mapmovhdr kk2.mgz    --sat 50'];
    %     system(cmd);
    
    aux=my_lta_read(['tp_' num2str(t) '_to_base.lta']);
    volumeFactors(t)=det(aux(1:3,1:3));  %  when >1, it means that hippo is bigger in base (so we'll divide by these at the end)
end
system('rm refTP*.lta');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we initialize the mesh deformations of the subject-specific     %
% atlas using the ASEGs from the base and the longitudinal FS stream. %
% For the time points, we just use the initialization of the base     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% we pretty much copy-paste from Koen's preprocessHippoSubfields.m
for t=0:nTP
    
    kvlClear % release memory
    
    if t==0
        asegFileName = [subjectDir '/' subjectBase '/mri/aseg.mgz'];  % FreeSurfer's volumetric segmentation results. This are non-probabilistic, "crisp" definitions
        meshCollectionFileName = atlasMeshFileName; % The tetrahedral atlas mesh
        warpedMeshCollectionFileName='warpedOriginalMeshBase.txt';
    else
        asegFileName = ['./asegRotated_tp_' num2str(t) '.mgz'];
        %         system([FSpath '/mri_vol2vol --mov ' subjectDir '/' subjectTPs{t} '/mri/aseg.mgz --o ' asegFileName ...
        %             ' --no-resample --targ ' subjectDir '/' subjectBase '/mri/aseg.mgz --lta tp_' num2str(t) '_to_base.lta']);
        system([FSpath '/mri_convert  ' subjectDir '/' subjectTPs{t} '/mri/aseg.mgz  ' asegFileName ...
            ' -rt nearest -odt float -rl ' subjectDir '/' subjectBase '/mri/aseg.mgz -at tp_' num2str(t) '_to_base.lta']);
        meshCollectionFileName = 'warpedOriginalMeshBase.txt';
        warpedMeshCollectionFileName=['warpedOriginalMeshTP' num2str(t) '.txt'];
    end
    boundingFileName = [tempdir 'imageDump.mgz']; % Bounding box
    compressionLookupTableFileName =compressionLUTfileName; % Look-up table belonging to the atlas
    [ FreeSurferLabels, names, colors ] = kvlReadCompressionLookupTable( compressionLookupTableFileName );
    
    % Prepare modified ASEG we will segment
    ASEG=myMRIread(asegFileName,0,tempdir);
    ASEG.vol(ASEG.vol==15)=0;  % 4th vent -> background (we're killing brainstem anyway...)
    ASEG.vol(ASEG.vol==16)=0;  % get rid of brainstem
    ASEG.vol(ASEG.vol==7)=0;   % get rid of left cerebellum WM ...
    ASEG.vol(ASEG.vol==8)=0;   % ... and of left cerebellum CT
    ASEG.vol(ASEG.vol==46)=0;  % get rid of right cerebellum WM ...
    ASEG.vol(ASEG.vol==47)=0;  % ... and of right cerebellum CT
    ASEG.vol(ASEG.vol==80)=0;  % non-WM hippo -> background
    ASEG.vol(ASEG.vol==85)=0;  % optic chiasm -> background
    ASEG.vol(ASEG.vol==72)=4;  % 5th ventricle -> left-lat-vent
    
    if strcmp(side,'left')>0
        ASEG.vol(ASEG.vol==5)=4;   % left-inf-lat-vent -> left-lat-vent
        ASEG.vol(ASEG.vol==30)=2;  % left-vessel -> left  WM
        ASEG.vol(ASEG.vol==14)=4;  % 3rd vent -> left-lat-vent
        ASEG.vol(ASEG.vol==24)=4;  % CSF -> left-lat-vent
        ASEG.vol(ASEG.vol==77)=2;  % WM hippoint -> left WM
        ASEG.vol(ASEG.vol>250)=2;  % CC labels -> left WM
        list2kill=[44 62 63 41 42 43 49 50 51 52 53 54 58 60];
        for k=1:length(list2kill)
            ASEG.vol(ASEG.vol==list2kill(k))=0;
        end
    else
        BU=ASEG.vol;
        ASEG.vol(:)=0;
        ASEG.vol(BU==44)=4; % right-inf-lat-vent -> left-lat-vent
        ASEG.vol(BU==62)=2; % right-vessel -> left  WM
        ASEG.vol(BU==14)=4;  % 3rd vent -> left-lat-vent
        ASEG.vol(BU==24)=4;  % CSF -> left-lat-vent
        ASEG.vol(BU==77)=2;  % WM hippoint -> left WM
        ASEG.vol(BU>250)=2;  % CC labels -> left WM
        % left to right
        ASEG.vol(BU==41)=2; % WM
        ASEG.vol(BU==42)=3; % CT
        ASEG.vol(BU==43)=4; % LV
        ASEG.vol(BU==49)=10; % TH
        ASEG.vol(BU==50)=11; % CA
        ASEG.vol(BU==51)=12; % PU
        ASEG.vol(BU==52)=13; % PA
        ASEG.vol(BU==53)=17; % HP
        ASEG.vol(BU==54)=18; % AM
        ASEG.vol(BU==58)=26; % AA
        ASEG.vol(BU==60)=28; % DC
        ASEG.vol(BU==63)=31; % CP
    end
    ASEG.vol(ASEG.vol==0)=1; % KVL code ignores zeros and we don't want that here
    myMRIwrite(ASEG,['asegMod_tp_' num2str(t) '.mgz'],'float',tempdir);
    
    % Eugenio July 2017
    % We now merge hippo, amygdala and cortex in cheating image
    ASEGcha=ASEG;
    ASEGcha.vol(ASEGcha.vol==17)=3;
    ASEGcha.vol(ASEGcha.vol==18)=3;
    % Eugenio, February 2022
    % myMRIwrite(ASEG,['asegModCHA_tp_' num2str(t) '.mgz'],'float',tempdir);
    myMRIwrite(ASEGcha,['asegModCHA_tp_' num2str(t) '.mgz'],'float',tempdir);
    
    
    % Read in modified aseg and also transform
    % Eugenio July 2017
    % [ synIm, transform ] = kvlReadCroppedImage( ['asegMod_tp_' num2str(t) '.mgz'], boundingFileName );
    [ synIm, transform ] = kvlReadCroppedImage( ['asegModCHA_tp_' num2str(t) '.mgz'], boundingFileName );
    
    synImBuffer = kvlGetImageBuffer( synIm );
    synSize = size( synImBuffer );
    if ~isdeployed && DEBUG>0
        figure
        showImage( synIm )
        title(['Synthetic Image: time point' num2str(t)])
    end
    
    % read in collection, set K and apply transform
    meshCollection = kvlReadMeshCollection( meshCollectionFileName );
    kvlTransformMeshCollection( meshCollection, transform );
    kvlSetKOfMeshCollection( meshCollection, Katl ); % liberal stiffness here
    
    % For the base, retrieve the reference mesh, for the time points, the deformed
    if t==0
        mesh = kvlGetMesh( meshCollection, -1 );
    else
        mesh = kvlGetMesh( meshCollection, 0 );
    end
    originalNodePositions = kvlGetMeshNodePositions( mesh );
    originalAlphas = kvlGetAlphasInMeshNodes( mesh );
    
    % Eugenio July 2017: head / body, and also merged hippos/amygdala/cortex, and moved fissure with background
    FreeSurferLabelGroups=[];
    FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex',... % cortex
        'Left-Hippocampus','alveus','subiculum-body','subiculum-head','Hippocampal_tail',... % hippo
        'molecular_layer_HP-body','molecular_layer_HP-head','GC-ML-DG-body','GC-ML-DG-head',...
        'CA4-body','CA4-head','CA1-body','CA1-head','CA3-body','CA3-head','HATA','fimbria',...
        'presubiculum-body','presubiculum-head','parasubiculum','Left-hippocampus-intensity-abnormality', ...
        'Left-Amygdala','Lateral-nucleus','Paralaminar-nucleus','Basal-nucleus',... % amygdala
        'Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
        'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus',...
        'Anterior-amygdaloid-area-AAA'};
    FreeSurferLabelGroups{end+1}={'Left-Cerebral-White-Matter'};
    FreeSurferLabelGroups{end+1}={'Left-Lateral-Ventricle'};
    FreeSurferLabelGroups{end+1}={'Left-choroid-plexus'};
    FreeSurferLabelGroups{end+1}={'Background','hippocampal-fissure','Background-CSF','Background-vessels','Background-tissue','Unknown'};
    FreeSurferLabelGroups{end+1}={'Left-VentralDC'};
    FreeSurferLabelGroups{end+1}={'Left-Putamen'};
    FreeSurferLabelGroups{end+1}={'Left-Pallidum'};
    FreeSurferLabelGroups{end+1}={'Left-Thalamus-Proper'};
    FreeSurferLabelGroups{end+1}={'Left-Accumbens-area'};
    FreeSurferLabelGroups{end+1}={'Left-Caudate'};
    FreeSurferLabelGroups{end+1}={'SUSPICIOUS'};
    
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
    
    % Eugenio July 2017 
    % In February 2022, I added:  & labels~=215
    cheatingMeans=zeros(length( sameGaussianParameters),1);
    cheatingVariances=0.01*ones(length( sameGaussianParameters),1);
    for l=1:length(sameGaussianParameters)
        labels= sameGaussianParameters{l};
        if any(labels>=200 & labels<=226 & labels~=215),  cheatingMeans(l)=3; %  cheatingMeans(l)=17; % HIPPO SF -> HIPPO
        elseif any(labels>=7000),  cheatingMeans(l)=3;  % cheatingMeans(l)=18; % AMYGDALOID SUBNUCLEI -> AMYGDALA
        elseif any(labels==0), cheatingMeans(l)=1; % BACKGROUND is 1 instead of 0
        elseif any(labels==999), cheatingMeans(l)=55; cheatingVariances(l)=55^2; % This is the generic, "suspicious" label we use for cysts...
        else cheatingMeans(l)=labels(1);
        end
    end
    
    % Compute the "reduced" alphas - those referring to the "super"-structures
    [ reducedAlphas, reducingLookupTable ] = kvlReduceAlphas( originalAlphas, compressionLookupTableFileName, sameGaussianParameters );
    if ( max( abs( sum( reducedAlphas, 2 ) - 1 ) ) > 1e-5 ) % Make sure these vectors really sum to 1
        error( 'The vector of prior probabilities in the mesh nodes must always sum to one over all classes' )
    end
    
    % Set the reduced alphas to be the alphas of the mesh
    kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
    
    % Create mask for valid area
    
    
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
    
    
    if ~isdeployed && DEBUG>0
        priors = kvlRasterizeAtlasMesh( mesh, synSize );
        figure
        for cheatingLabel = 1 : size( reducedAlphas, 2 )
            subplot( 3, 5, cheatingLabel )
            showImage( priors( :, :, :, cheatingLabel ) )
        end
        title(['Time point ' num2str(t) ': priors for segmentation of fake intensity image']);
    end
    
    % Create a mask with the valid image area
    
    cheatingImageBuffer=synImBuffer;
    cheatingImageBuffer(~MASK)=0;
    cheatingImage = kvlCreateImage( cheatingImageBuffer );
    
    if ~isdeployed && DEBUG>0
        figure
        showImage( cheatingImage, [], [ 0 110 ] )  % Using dynamic range for display that shows what's going on
        title(['Time point ' num2str(t) ': fake intensity image to initialize atlas deformation']);
    end
    
    % We use a multiscale approach for the base, only the finer resolution for
    % the time points
    if t==0
        meshSmoothingSigmas = [ 3.0 2.0]';
    else
        % meshSmoothingSigmas = [ 2.0]';
        meshSmoothingSigmas = [ ]';
    end
    
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
                subplot( 3, 5, cheatingLabel )
                showImage( priors( :, :, :, cheatingLabel ) )
            end
            title(['Time point ' num2str(t) ': Smoothed priors']);
        end
        
        % Eugenio November 2017: GEMS2
        % Set up the black box optimizer for the mesh nodes
        if ( exist( 'cheatingOptimizer', 'var' ) == 1 )
            % The optimizer is very memory hungry when run in multithreaded mode.
            % Let's clear any old ones we may have lying around
            % Eugenio November 2017: GEMS2
            try
                kvlClear( cheatingOptimizer );
                kvlClear( cheatingCalculator );
            end
        end
        
        
        % Eugenio November 2017: GEMS2   (note that it uses variances instead of precisions)
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
        
        
        % Main optimization loop
        for positionUpdatingIterationNumber = 1 : maxpuin
            disp(['Time point ' num2str(t) ' resolution ' num2str(multiResolutionLevel) ', iteration ' num2str(positionUpdatingIterationNumber)]);
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
                title(['Time point ' num2str(t) ': current atlas deformation'])
                subplot( 2, 2, 4 )
                plot( historyOfMinLogLikelihoodTimesPrior( 2 : end ) )
                title(['Time point ' num2str(t) ': history of Log-lhood + Log-prior'])
                drawnow
            end
        end
    end
    try  % Eugenio November 2017
        kvlClear( cheatingOptimizer );
        kvlClear( cheatingCalculator );
    end
    
    
    disp(['Time point ' num2str(t) ': fitting mesh to synthetic image from ASEG took ' num2str(etime(clock,time_ref_cheat_optimization)) ' seconds']);
    
    if positionUpdatingIterationNumber==1
        error('Fitting mesh to synthetic image resulted in no deformation')
    end
    
    % OK, we're done. let's modify the mesh atlas in such a way that our computed mesh node positions are
    % assigned to what was originally the mesh warp corresponding to the first training subject.
    % Also, the reference for the time points is set to the deformed base
    kvlSetAlphasInMeshNodes( mesh, originalAlphas )
    updatedNodePositions = kvlGetMeshNodePositions( mesh );
    kvlSetMeshCollectionPositions( meshCollection, originalNodePositions,  updatedNodePositions );
    
    
    % Compare the average shape we started with, with the shape we have computed now in a little movie
    if ~isdeployed && DEBUG>0
        figure
        originalPositionColorCodedPriors = ...
            kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( kvlGetMesh( meshCollection, -1 ), synSize ), colors );
        updatedPositionColorCodedPriors = ...
            kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( kvlGetMesh( meshCollection, 0 ), synSize ), colors );
        for i=1:20
            subplot(1,3,1)
            showImage( cheatingImage ), title(['Time point ' num2str(t) ': cheating image'])
            subplot(1,3,2)
            showImage( cheatingImage ), title(['Time point ' num2str(t) ': cheating image'])
            subplot(1,3,3)
            showImage( originalPositionColorCodedPriors ), title(['Time point ' num2str(t) ': prior before deformation'])
            pause(.5)
            subplot(1,3,1)
            showImage( originalPositionColorCodedPriors ), title(['Time point ' num2str(t) ': prior before deformation'])
            subplot(1,3,2)
            showImage( updatedPositionColorCodedPriors ), title(['Time point ' num2str(t) ': prior after deformation'])
            subplot(1,3,3)
            showImage( updatedPositionColorCodedPriors ), title(['Time point ' num2str(t) ': prior after deformation'])
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
    kvlWriteMeshCollection( meshCollection, warpedMeshCollectionFileName );
    
end  %% end of loop over time points

disp('Done with initializations');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we do the iterative segmentation of the time points %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


kvlClear % Clear all the wrapped C++ stuff
close all


% Get data ready: image at work resolution, masking out non-brain voxels and also the cerebellum
% We also mask out brainstem and 3rd/4th ventricles, which can be annoying later on
disp('Extracting blocks from norm.mgz and upsample to work resolution');
L=myMRIread([subjectDir '/' subjectBase '/mri/aseg.mgz'],0,tempdir);
margin=15; % in mm
if strcmp(side,'left')
    [~,cropping]=cropLabelVol(L.vol==HippoLabelLeft,round(margin/mean(L.volres)));
else
    [~,cropping]=cropLabelVol(L.vol==HippoLabelRight,round(margin/mean(L.volres)));
end
Lcrop=applyCropping(L.vol,cropping);
offsetVox=[cropping(2)-1;cropping(1)-1;cropping(3)-1;0];
RAScorner=L.vox2ras0*offsetVox;
vox2ras0New=[ L.vox2ras0(:,1:3)  L.vox2ras0(:,4)+RAScorner];
A=myMRIread([subjectDir '/' subjectBase '/mri/norm.mgz'],0,tempdir);
Icrop=applyCropping(A.vol,cropping);
I=L;
I.vol=Icrop;
I.vox2ras0=vox2ras0New;
myMRIwrite(I,'tempT1.mgz','float',tempdir);
cmd=[FSpath '/mri_convert tempT1.mgz T1resampled_0.mgz  -odt float -rt  cubic -vs ' num2str(resolution) ' '  num2str(resolution)  ' '  num2str(resolution)];
system([cmd ' >/dev/null']);
delete tempT1.mgz
system([FSpath '/mri_binarize --i asegMod_tp_0.mgz --min 1.5 --dilate 2 --o asegModBinDilated.mgz >/dev/null']);
system([FSpath '/mri_convert asegModBinDilated.mgz asegModBinDilatedResampled.mgz -odt float -rt nearest -rl T1resampled_0.mgz >/dev/null']);
system([FSpath '/mri_mask -T 0.5 T1resampled_0.mgz asegModBinDilatedResampled.mgz T1resampled_0.mgz >/dev/null']);

% Eugenio: let's try masking anything that is not close to the hippo,
% brainstem and 3rd/4th ventricles, which can be annoying later on.
dilSize=round(5/mean(A.volres));
system([FSpath '/mri_binarize --i asegMod_tp_0.mgz --min 16.5 --max 18.5 --o hippoMask.mgz >/dev/null']);
system([FSpath '/mri_binarize --i asegMod_tp_0.mgz --min 16.5 --max 18.5 --o hippoMaskDilated.mgz --dilate ' num2str(dilSize) ' >/dev/null']);
system([FSpath '/mri_convert hippoMask.mgz hippoMaskResampled.mgz -rt interpolate -rl T1resampled_0.mgz  -odt float >/dev/null']);
system([FSpath '/mri_binarize --i  hippoMaskResampled.mgz --min 0.5 --dilate ' num2str(round(3/resolution)) '  --o hippoMaskResampledDilated.mgz  >/dev/null']);
system([FSpath '/mri_mask -T 0.5 T1resampled_0.mgz  hippoMaskResampledDilated.mgz T1resampled_0.mgz >/dev/null']);

[ imageBase, transform ] = kvlReadCroppedImage( 'T1resampled_0.mgz', boundingFileName );
imageBaseBuffer = kvlGetImageBuffer( imageBase );
imageSize = size( imageBaseBuffer );

images=cell([1,nTP]);
imageBuffers=cell([1,nTP]);
for t=1:nTP
    normFileName=[subjectDir '/' subjectTPs{t} '/mri/norm.mgz'];
    resampledFileName=['T1resampled_' num2str(t) '.mgz'];
    
    system([FSpath '/mri_convert ' normFileName ' ' resampledFileName ' -odt float -rl T1resampled_0.mgz -rt cubic -at tp_' num2str(t) '_to_base.lta  >/dev/null']);
    system([FSpath '/mri_mask -T 0.5 ' resampledFileName '  hippoMaskResampledDilated.mgz ' resampledFileName ' >/dev/null']);
    [ images{t}, ~ ] = kvlReadCroppedImage( resampledFileName, boundingFileName );
    imageBuffers{t} = kvlGetImageBuffer( images{t} );
end
if ~isdeployed && DEBUG>0
    figure, showImage(imageBaseBuffer), title('Base image');
    for t=1:nTP
        figure, showImage( imageBuffers{t} ); title(['Buffer to segment: time point ' num2str(t)]);
    end
end


% We're not interested in image areas that fall outside our cuboid ROI where our atlas is defined. Therefore,
% generate a mask of what's inside the ROI. Also, by convention we're skipping all voxels whose intensities
% is exactly zero
disp('Further masking of image buffers with mesh ROI');
meshCollection = kvlReadMeshCollection( 'warpedOriginalMeshBase.txt.gz' );
kvlTransformMeshCollection( meshCollection, transform );
mesh = kvlGetMesh( meshCollection, 0 );

% Eugenio July 2017
% priors = kvlRasterizeAtlasMesh( mesh, imageSize ); % Without specifying a specific label, will rasterize all simultaneously
% mask = imerode(single( sum( priors, 4 ) / 65535 ) > 0.99,createSphericalStrel(5));
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
mask = imerode(single( sumpriors / 65535 ) > 0.99,createSphericalStrel(5));


mask = mask & ( imageBaseBuffer > 0 );
maskIndices = find( mask );
kvlClear(meshCollection); kvlClear(mesh);

imageBaseBuffer( find( ~mask ) ) = 0;
kvlSetImageBuffer( imageBase, imageBaseBuffer );
imageBaseBuffer = kvlGetImageBuffer( imageBase );
for t=1:nTP
    imageBuffers{t}( find( ~mask ) ) = 0;
    kvlSetImageBuffer( images{t}, imageBuffers{t} );
    imageBuffers{t} = kvlGetImageBuffer( images{t} );
end

if ~isdeployed && DEBUG>0
    colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
    figure
    showImage( colorCodedPriors );
    title('Color-coded priors of base')
    clear colorCodedPriors
end
if ~isdeployed && DEBUG>0
    for t=1:nTP
        figure
        subplot( 1, 2, 1 )
        showImage( mask )
        title('Mask given by mesh of base')
        subplot(1,2,2)
        showImage(imageBuffers{t})
        title(['Time point ' num2str(t) ': masked image buffer']);
    end
end


% Merge classes with similar intensity properties
disp('Merging classes with similar intensity properties ("reduced alphas")');

% Eugenio July 2017
FreeSurferLabelGroups=[];
if highres==0
    FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex','Left-Hippocampus','Left-Amygdala','subiculum-head','subiculum-body','Hippocampal_tail','GC-ML-DG-head','GC-ML-DG-body','CA4-head','CA4-body','presubiculum-head','presubiculum-body',...
        'CA1-head','CA1-body','parasubiculum','CA3-head','CA3-body','HATA','Lateral-nucleus','Paralaminar-nucleus',...
        'Basal-nucleus','Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
        'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus',...
        'Anterior-amygdaloid-area-AAA','molecular_layer_HP-body','molecular_layer_HP-head'};
else
    FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex','Left-Hippocampus','Left-Amygdala','subiculum-head','subiculum-body','Hippocampal_tail','GC-ML-DG-head','GC-ML-DG-body','CA4-head','CA4-body','presubiculum-head','presubiculum-body',...
        'CA1-head','CA1-body','parasubiculum','CA3-head','CA3-body','HATA','Lateral-nucleus','Paralaminar-nucleus',...
        'Basal-nucleus','Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
        'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus',...
        'Anterior-amygdaloid-area-AAA'};
    FreeSurferLabelGroups{end+1}={'molecular_layer_HP-body','molecular_layer_HP-head'};
end
FreeSurferLabelGroups{end+1}={'Left-Cerebral-White-Matter','fimbria'};
FreeSurferLabelGroups{end+1}={'alveus'};
FreeSurferLabelGroups{end+1}={'Left-Lateral-Ventricle','Background-CSF','SUSPICIOUS','Left-hippocampus-intensity-abnormality'};
FreeSurferLabelGroups{end+1}={'hippocampal-fissure'};
FreeSurferLabelGroups{end+1}={'Left-Pallidum'};
FreeSurferLabelGroups{end+1}={'Left-Putamen'};
FreeSurferLabelGroups{end+1}={'Left-Caudate'};
FreeSurferLabelGroups{end+1}={'Left-Thalamus-Proper'};
FreeSurferLabelGroups{end+1}={'Left-choroid-plexus'};
FreeSurferLabelGroups{end+1}={'Left-VentralDC'};
FreeSurferLabelGroups{end+1}={'Left-Accumbens-area'};
FreeSurferLabelGroups{end+1}={'Unknown','Background-tissue'};



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

% Compute hyperparameters for estimation of Gaussian parameters
disp('Computing hyperparameters for estimation of Gaussian parameters')
nHyper=cell([1,t]);
meanHyper=cell([1,t]);
for t=1:nTP
    DATA=myMRIread([subjectDir '/' subjectTPs{t} '/mri/norm.mgz'],0,tempdir);
    WMPARC=myMRIread([subjectDir '/' subjectTPs{t} '/mri/wmparc.mgz'],0,tempdir);
    aux=myMRIread('hippoMaskDilated.mgz',0,tempdir);
    WMPARC.vol(WMPARC.vol==0 & aux.vol==0)=-1;
    
    nHyper{t}=zeros(length(sameGaussianParameters),1);
    meanHyper{t}=zeros(length(sameGaussianParameters),1);
    BGindex=0;
    for g=1:length(sameGaussianParameters)
        labels=sameGaussianParameters{g};
        if any(labels==0)
            BGindex=g;
        end
        if any(labels==3 | labels==17 | labels==18 | labels > 7000 | labels==226)  % gray matter + alveus/ML
            if strcmp(side,'left')>0, listMask=17; else, listMask=53; end
            % listMask=[3 42 17 53 18 54];
        elseif any(labels==2)  % white matter
            if strcmp(side,'left')>0, listMask=[3006 3007 3016]; else, listMask=[4006 4007 4016]; end
            % listMask=[2 41];
        elseif any(labels==26)    % accumbens area
            if strcmp(side,'left')>0, listMask=26; else, listMask=58; end
            % listMask=[26 58];
        elseif any(labels==4)  % CSF
            if strcmp(side,'left')>0, listMask=4; else, listMask=43; end
            % listMask=[4 43 14 15] ;
        elseif any(labels==0) % Background
            listMask=0;
        elseif any(labels==13)
            if strcmp(side,'left')>0, listMask=13; else, listMask=52; end
            % listMask=[13 52];
        elseif any(labels==12)
            if strcmp(side,'left')>0, listMask=12; else, listMask=51; end
            % listMask=[12 51];
        elseif any(labels==11)
            if strcmp(side,'left')>0, listMask=11; else, listMask=50; end
            % listMask=[11 50];
        elseif any(labels==10)
            if strcmp(side,'left')>0, listMask=10; else, listMask=49; end
            % listMask=[10 49];
        elseif any(labels==31)
            if strcmp(side,'left')>0, listMask=31; else, listMask=63; end
            % listMask=[31 63];
        elseif any(labels==28)
            if strcmp(side,'left')>0, listMask=28; else, listMask=60; end
            % listMask=[28 60];
        else
            listMask=[];
        end
        if length(listMask)>0
            MASK=zeros(size(DATA.vol));
            for l=1:length(listMask)
                MASK=MASK | WMPARC.vol==listMask(l);
            end
            MASK=imerode(MASK,createSphericalStrel(round(1/mean(DATA.volres))));
            data=DATA.vol(MASK & DATA.vol>0);
            meanHyper{t}(g)=median(data);
            nHyper{t}(g)=10+length(data)*prod(DATA.volres)/resolution^3;
        end
    end
    % if any nan, replace by background
    ind=find(isnan(meanHyper{t}));
    meanHyper{t}(ind)=55;
    nHyper{t}(ind)=10;
end


% Here's the part where we simulate partial voluming!
% disp('Estimating typical intensities of molecular layer and alveus')
disp('Estimating typical intensities of  alveus and fissure')
WMind=-1;
GMind=-1;
ALind=-1;
MLind=-1;
FISSind=-1;
CSFind=-1;
for g=1:length(sameGaussianParameters)
    labels=sameGaussianParameters{g};
    if any(labels==2)
        WMind=g;
    end
    if any(labels==3)
        GMind=g;
    end
    if any(labels==201)
        ALind=g;
    end
    if any(labels==245) && highres>0  % Eugenio July 2017 (changed 214 by 245)
        MLind=g;
    end
    if any(labels==215)
        FISSind=g;
    end
    if any(labels==4)
        CSFind=g;
    end
end

meshCollection = kvlReadMeshCollection( 'warpedOriginalMeshBase.txt.gz' );
kvlTransformMeshCollection( meshCollection, transform );
mesh = kvlGetMesh( meshCollection, 0 );
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )

% Eugenio July 2017: again, rasterize priors one at the time
for l = 1 : size(reducedAlphas,2)
    if l==1
        % This is a bit annoying, but the call to kvlRasterize with a single
        % label fills in the voxels outside the cuboid with p=1 (whereas the
        % call with multiple labels does not)
        sillyAlphas=zeros([size(reducedAlphas,1),2],'single');
        sillyAlphas(:,1)=reducedAlphas(:,1);
        sillyAlphas(:,2)=1-reducedAlphas(:,1);
        kvlSetAlphasInMeshNodes( mesh, sillyAlphas )
        prior = kvlRasterizeAtlasMesh( mesh, imageSize);
        kvlSetAlphasInMeshNodes( mesh, reducedAlphas )
        prior=prior(:,:,:,1);
        sumpriors=prior;
        L=ones(size(prior));
        PMAX=prior;
    else
        prior = kvlRasterizeAtlasMesh( mesh, imageSize, l-1 );
        sumpriors=sumpriors+prior;
        M=prior>PMAX;
        L(M)=l;
        PMAX(M)=prior(M);
    end
end
suma=single(sumpriors)/65535;
maskPriors=suma>.97;
kvlClear(meshCollection); kvlClear(mesh);

% priors = kvlRasterizeAtlasMesh( mesh, imageSize ); % Without specifying a specific label, will rasterize all simultaneously
% kvlClear(meshCollection); kvlClear(mesh);
% priors=double(priors)/63535;
% suma=sum(priors,4);
% maskPriors=suma>.97;
% priors=priors./(eps+repmat(suma,[1 1 1 numberOfClasses]));
% [~,L]=max(priors,[],4);




for t=1:nTP
    I=zeros(size(L));
    for l=1:numberOfClasses
        if l==ALind || l==MLind
            I(L==l)=meanHyper{t}(WMind);
        elseif l==FISSind
            I(L==l)=meanHyper{t}(CSFind);
        else
            I(L==l)=meanHyper{t}(l);
        end
    end
    I(~maskPriors)=0;
    I_PV=GaussFilt3d(I,mean(DATA.volres)/(2.355*resolution));
    
    if ALind~=-1
        data=I_PV(L==ALind); % it's multimodal, so median won't cut it...
        [density,v]=ksdensity(data);
        [trash,idx]=max(density);
        meanHyper{t}(ALind)=median(v(idx));
        nHyper{t}(ALind)=(nHyper{t}(GMind)+nHyper{t}(WMind))/2;
    end
    
    if highres>0
        data=I_PV(L==MLind);
        meanHyper{t}(MLind)=median(data);
        nHyper{t}(MLind)=(nHyper{t}(WMind)+nHyper{t}(GMind))/2;
    end
    
    if FISSind~=-1
        data=I_PV(L==FISSind);
        meanHyper{t}(FISSind)=median(data);
        nHyper{t}(FISSind)=(nHyper{t}(CSFind)+nHyper{t}(GMind))/2;
    end
    
    
end


% Read in meshes: base, time points, and auxiliary (read from base)
meshCollectionBase = kvlReadMeshCollection( 'warpedOriginalMeshBase.txt.gz' );
kvlTransformMeshCollection( meshCollectionBase, transform );
kvlSetKOfMeshCollection( meshCollectionBase, Katl );
meshBase = kvlGetMesh( meshCollectionBase, 0 );
subjectAtlasPositions = kvlGetMeshNodePositions( meshBase );
kvlSetAlphasInMeshNodes( meshBase, reducedAlphas )
meshOrig = kvlGetMesh( meshCollectionBase,-1);
atlasPositions = kvlGetMeshNodePositions(meshOrig);
kvlClear(meshOrig);


meshCollections=cell([1,nTP]);
meshes=cell([1,nTP]);
subjectTPpositions=cell([1,nTP]);
for t=1:nTP
    meshCollections{t} = kvlReadMeshCollection(['warpedOriginalMeshTP' num2str(t) '.txt.gz']);
    kvlTransformMeshCollection( meshCollections{t}, transform );
    kvlSetKOfMeshCollection( meshCollections{t}, Ktp );
    meshes{t} = kvlGetMesh( meshCollections{t}, 0 );
    subjectTPpositions{t} = kvlGetMeshNodePositions( meshes{t} );
    kvlSetAlphasInMeshNodes( meshes{t}, reducedAlphas )
end

% This will be the mesh collection we'll use for temporary data
meshCollection = kvlReadMeshCollection( 'warpedOriginalMeshBase.txt.gz' );
kvlTransformMeshCollection( meshCollection, transform );
kvlSetKOfMeshCollection( meshCollection, Katl );
mesh = kvlGetMesh( meshCollection, 0 );
nodePositions = kvlGetMeshNodePositions( mesh );
kvlSetAlphasInMeshNodes( mesh, reducedAlphas )


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL LOOP: ATLAS ESTIMATION - SEGMENTATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

globalIt=0;

globalReady=0;

if FAST>0, MAX_GLOBAL_ITS=2; end

while globalReady==0
    
    globalIt = globalIt+1;
    
    disp(' ');
    disp('*******************');
    disp([' Global iteration ' num2str(globalIt) ]);
    disp('*******************');
    disp(' ');
    disp('   ***** Estimating subject-specific atlas *****');
    
    
    % Eugenio: November 2017: switched from kvlAverageMeshes to Matlab
    % interface 
    kvlSetMeshCollectionPositions(meshCollection,atlasPositions,subjectAtlasPositions); 
    meshSA=kvlGetMesh(meshCollection,0);
        
    cmd='kvlSetMeshCollectionPositions(meshCollection,atlasPositions';
    for t=1:nTP
        cmd=[cmd ',subjectTPpositions{' num2str(t) '}'];
    end
    cmd=[cmd ');'];
    eval(cmd);
    
    if ( exist( 'optimizer', 'var' ) == 1 )
        % The optimizer is very memory hungry when run in multithreaded mode.
        % Let's clear any old ones we may have lying around
        try
            kvlClear( optimizer );
            kvlClear( calculator );
            clear optimizer
            clear calculator
        catch ME
        end
    end
    
    calculator = kvlGetAverageAtlasMeshPositionCostAndGradientCalculator(meshCollection,Katl,Ktp,transform);
    optimizer = kvlGetOptimizer( optimizerType, meshSA, calculator, ...
        'Verbose', verbose, ...
        'MaximalDeformationStopCriterion', maximalDeformationStopCriterion, ...
        'LineSearchMaximalDeformationIntervalStopCriterion', ...
        lineSearchMaximalDeformationIntervalStopCriterion, ...
        'MaximumNumberOfIterations', maximumNumberOfDeformationIterations, ...
        'BFGS-MaximumMemoryLength', BFGSMaximumMemoryLength );
    

    hoc=[];
    maxitsatlasupdate=400;
    for positionUpdatingIterationNumber = 1 : maxitsatlasupdate
        disp(['Subject atlas deformation, step ' num2str(positionUpdatingIterationNumber) ' of ' num2str(maxitsatlasupdate)]);
        
        % Eugenio May2018
        maximalDeformation=0;
        try
            tic
            [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer );
            elapsedTime = toc;
            disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
            minLogLikelihoodTimesPrior
            hoc(end+1)=minLogLikelihoodTimesPrior;
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
              
    subjectAtlasPositions = kvlGetMeshNodePositions(meshSA);
    kvlSetMeshCollectionPositions(meshCollectionBase, atlasPositions, subjectAtlasPositions);
    kvlSetKOfMeshCollection( meshCollectionBase, Katl );
    for t=1:nTP
        kvlSetMeshCollectionPositions(meshCollections{t}, subjectAtlasPositions, subjectTPpositions{t} );
        kvlSetKOfMeshCollection( meshCollections{t}, Ktp );
    end
   
    
    
%     % Prepare mesh collection with atlas and current deformations (input to
%     % kvlAverageMeshes)
%     
%     cmd='kvlSetMeshCollectionPositions(meshCollection,atlasPositions';
%     for t=1:nTP
%         cmd=[cmd ',subjectTPpositions{' num2str(t) '}'];
%     end
%     cmd=[cmd ');'];
%     eval(cmd);
%     
%     kvlWriteMeshCollection(meshCollection,'input.txt');
%     
%     
%     
%     kvlSetMeshCollectionPositions(meshCollection,atlasPositions,subjectAtlasPositions);
%     kvlWriteMeshCollection(meshCollection,'initialization.txt');
%     
%     % Call kvlAverageMeshes to get the subject-specific atlas
%     cmd=[FSpath '/kvlAverageMeshes'];
%     if exist(cmd,'file')==0
%         cmd='/cluster/koen/eugenio/GEMS-Release-linux/bin/kvlAverageMeshes';
%         if exist(cmd,'file')==0
%             error('Command "kvlAverageMeshes" not found')
%         end
%     end
%     
%     
%     cmd=[cmd ' input.txt initialization.txt ' ...
%         num2str(Katl) ' ' num2str(Ktp) ' SubjectAtlas.txt ' num2str(nTP)];
%     system(cmd);
%     % system([cmd ' >/dev/null']);
%     
%     % Read output, and update the reference positions of the time points
%     mc = kvlReadMeshCollection('SubjectAtlas.txt.gz');
%     m = kvlGetMesh(mc,0);
%     subjectAtlasPositions = kvlGetMeshNodePositions( m );
%     kvlClear(mc); kvlClear(m);
%     kvlSetMeshCollectionPositions(meshCollectionBase, atlasPositions, subjectAtlasPositions);
%     kvlSetKOfMeshCollection( meshCollectionBase, Katl );
%     for t=1:nTP
%         kvlSetMeshCollectionPositions(meshCollections{t}, subjectAtlasPositions, subjectTPpositions{t} );
%         kvlSetKOfMeshCollection( meshCollections{t}, Ktp );
%     end
    
    
    disp(' ');
    disp('   ***** Done with estimation of subject-specific atlas *****');
    disp( ' ');
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Time to deform meshes %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    disp(' ');
    disp('   ***** Deformation of meshes according to image data *****');
    disp(' ');
    
    for t = 1:nTP  % loop around time points
        
        disp(' ');
        disp(['  Deforming time point ' num2str(t) ' of ' num2str(nTP)]);
        disp(' ');
        
        % Multi-resolution scheme
        if globalIt <= 2
            meshSmoothingSigmas = [ 1.5 .75 ]';
            imageSmoothingSigmas = [0 0]';
            maxItNos=[6 3];  % each iteration has 40 deformation steps
        else
            meshSmoothingSigmas = [ .75 0 ]';
            imageSmoothingSigmas = [0 0]';
            maxItNos=[2 1];  % each iteration has 40 deformation steps
        end
        
        numberOfMultiResolutionLevels = length( meshSmoothingSigmas );
        
        time_ref_optimization=clock;
        
        imageBufferOrig=imageBuffers{t};
        
        for multiResolutionLevel = 1 : numberOfMultiResolutionLevels
            
            % Smooth the mesh using a Gaussian kernel.
            meshes{t} = kvlGetMesh( meshCollections{t}, 0 );
            kvlSetAlphasInMeshNodes( meshes{t}, reducedAlphas )
            meshSmoothingSigma = meshSmoothingSigmas( multiResolutionLevel );
            fprintf( 'Smoothing mesh collection with kernel size %f ...', meshSmoothingSigma )
            kvlSmoothMeshCollection( meshCollections{t}, meshSmoothingSigma )
            fprintf( 'done\n' )
            
            % Smooth the image using a Gaussian kernel
            imageSigma=imageSmoothingSigmas(multiResolutionLevel);
            if imageSigma>0
                imageBuffers{t}=single(GaussFilt3dMask(imageBufferOrig,imageBufferOrig>0,imageSigma,resolution*ones(1,3)));
                kvlSetImageBuffer(images{t},imageBuffer);
            end
            
            
            % Now with this smoothed atlas, we're ready for the real work. There are essentially two sets of parameters
            % to estimate in our generative model: (1) the mesh node locations (parameters of the prior), and (2) the
            % means and variances of the Gaussian intensity models (parameters of the
            % likelihood function, which is really a hugely simplistic model of the MR imaging process). Let's optimize
            % these two sets alternately until convergence. Optimizing the mesh node locations
            % is provided as a black box type of thing as it's implemented in C++ using complicated code - the other
            % set is much much better to experiment with in Matlab.
            
            
            
            % Maximum number of iterations (includes one imaging model parameter estimation and
            % one deformation optimization; the latter always does 20 steps).
            maximumNumberOfIterations = maxItNos(multiResolutionLevel);
            
            if FAST>0, maximumNumberOfIterations = 1; end  % in case we just wanna cruise throught it :-)
            
            historyOfCost = [ 1/eps ];
            
            % Compute a color coded version of the atlas prior in the atlas's current pose, i.e., *before*
            % we start deforming. We'll use this just for visualization purposes
            if  ~isdeployed && DEBUG>0
                deformationMovieFigure = 25;
                oldColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( meshes{t}, imageSize ) );
                figure( deformationMovieFigure )
                showImage( oldColorCodedPriors )
            end
            
            
            for iterationNumber = 1 : maximumNumberOfIterations
                disp(['Iteration ' num2str(iterationNumber) ' of ' num2str(maximumNumberOfIterations)]);
                
                %
                % Part I: estimate Gaussian mean and variances using EM
                %
                
                % Get the priors as dictated by the current mesh position, as well as the image intensities
                data = double( reshape( kvlGetImageBuffer( images{t} ), [ prod( imageSize ) 1 ] ) ); % Easier to work with vector notation in the computations
                
                % Eugenio July 2017: again, avoid spike of memory use
                priors=zeros([length(maskIndices),numberOfClasses],'uint16');
                for l=1:numberOfClasses
                    prior = kvlRasterizeAtlasMesh( meshes{t}, imageSize, l-1 );
                    priors(:,l)=prior(maskIndices);
                end
                
                %                 priors = kvlRasterizeAtlasMesh( meshes{t}, imageSize );
                %                 priors = reshape( priors, [ prod( imageSize ) numberOfClasses ] ); % Easier to work with vector notation in the computations
                %                 priors = priors( maskIndices, : );
                
                
                data = data( maskIndices );
                
                % Start EM iterations. Initialize the parameters if this is the
                % first time ever you run this
                EPS=1e-2;
                posteriors = double( priors ) / 65535;
                if ( ( multiResolutionLevel == 1) & ( iterationNumber == 1 ) & (globalIt==1) )
                    
                    for classNumber = 1 : numberOfClasses
                        
                        posterior = posteriors( :, classNumber );
                        
                        if sum(posterior)>EPS
                            
                            mu = (meanHyper{t}(classNumber)*nHyper{t}(classNumber) + data'*posterior) / ( nHyper{t}(classNumber) + sum( posterior ) + EPS );
                            variance = (( ( data - mu ).^2 )' * posterior + nHyper{t}(classNumber)*(mu-meanHyper{t}(classNumber))^2 )/ ( sum( posterior ) + EPS );
                            
                            means{t}( classNumber ) = mu;
                            variances{t}( classNumber ) = variance+EPS;
                            
                        else
                            means{t}( classNumber ) = meanHyper{t}(classNumber);
                            variances{t}( classNumber ) = 100;
                            
                        end
                        
                    end
                    variances{t}(variances{t}==0)=100; % added by Eugenio, prevents nans...
                    
                end % End test need for initialization
                
                stopCriterionEM = 1e-5;
                historyOfEMCost = [ 1/eps ];
                
                for EMIterationNumber = 1 : 100
                    %
                    % E-step: compute the posteriors based on the current parameters
                    %
                    minLogLikelihood = 0;
                    for classNumber = 1 : numberOfClasses
                        mu = means{t}( classNumber );
                        variance = variances{t}( classNumber );
                        prior = single( priors( :, classNumber ) ) / 65535;
                        posteriors( :, classNumber ) = ( exp( -( data - mu ).^2 / 2 / variance ) .* prior ) ...
                            / sqrt( 2 * pi * variance);
                        
                        minLogLikelihood = minLogLikelihood+0.5*log(2*pi*variance)-0.5*log(nHyper{t}(classNumber))...
                            +0.5*nHyper{t}(classNumber)/variance*(mu-meanHyper{t}(classNumber))^2; % contribution from prior
                        
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
                    
                    
                    
                    % Also show EM cost function
                    if  ~isdeployed && DEBUG>0
                        EMcostHistoryFigure = 30;
                        figure(EMcostHistoryFigure);
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
                            
                            mu = (meanHyper{t}(classNumber)*nHyper{t}(classNumber) + data'*posterior) / ( nHyper{t}(classNumber) + sum( posterior ) + EPS );
                            variance = (( ( data - mu ).^2 )' * posterior + nHyper{t}(classNumber)*(mu-meanHyper{t}(classNumber))^2 )/ ( sum( posterior ) + EPS );
                            
                            means{t}( classNumber ) = mu;
                            variances{t}( classNumber ) = variance+EPS;
                            
                        else
                            means{t}( classNumber ) = meanHyper{t}(classNumber);
                            variances{t}( classNumber ) = 100;
                            
                        end
                    end
                    variances{t}(variances{t}==0)=100; % added by Eugenio, prevents nans...
                    
                end % End EM iterations
                means{t}
                (variances{t}).^2
                
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
              
                
                
                % Eugenio November 2011
                if ( exist( 'optimizer', 'var' ) == 1 )
                    % The optimizer is very memory hungry when run in multithreaded mode.
                    % Let's clear any old ones we may have lying around
                    try
                        kvlClear( optimizer );
                        kvlClear( calculator );
                        clear optimizer
                        clear calculator
                    catch ME
                    end
                end
                
                % For some reason, clearing the optimizer also kills the meshes
                meshes{t} = kvlGetMesh( meshCollections{t}, 0 );
                
                haveMoved = false; % Keep track if we've ever moved or not
                
                % Eugenio November 2017: GEMS2, more iterations
                calculator = kvlGetCostAndGradientCalculator('AtlasMeshToIntensityImage',...
                    images{t}, 'Sliding',transform,means{t}',variances{t}',ones(size(means{t}')),ones(size(means{t}')));
                
                verbose=0;
                maximalDeformationStopCriterion=1e-10;
                lineSearchMaximalDeformationIntervalStopCriterion=1e-10;
                maximumNumberOfDeformationIterations=1000;
                BFGSMaximumMemoryLength=12;
                positionUpdatingMaximumNumberOfIterations=30;
                
                optimizer = kvlGetOptimizer( optimizerType, meshes{t}, calculator, ...
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
                        % Eugenio November 2011
                        [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlStepOptimizer( optimizer );
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
                    newColorCodedPriors = kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( meshes{t}, imageSize ) );
                    for i=1:10
                        showImage( oldColorCodedPriors ), title('Previous priors')
                        drawnow
                        pause( 0.1 )
                        showImage( newColorCodedPriors ), title('Updated priors')
                        drawnow
                        pause( 0.1 )
                    end
                end
                
                
                % Keep track of the cost function we're optimizing
                historyOfCost = [ historyOfCost; minLogLikelihoodTimesPrior ];
                if  ~isdeployed && DEBUG>0
                    costFigure = 30;
                    figure( costFigure )
                    plot( historyOfCost( 2 : end ) )
                    title( 'Cost' )
                end
                
                % Determine if we should stop the overall iterations over the two set of parameters
                if ( ( ~haveMoved ) || ( ( ( historyOfCost( end-1 ) - historyOfCost( end ) ) / historyOfCost( end ) ) <  1e-6 ) )
                    % Converged
                    break;
                end
                
                
            end % End looping over EM-def iterations
            
        end % End loop over multiresolution levels
        
        disp(['Fitting mesh to image data mask took, at this global iteration and time point, ' num2str(etime(clock,time_ref_optimization)) ' seconds']);
        
        % Restore original image buffer
        kvlSetImageBuffer(images{t},imageBufferOrig);
        
        % Clear some memory
        % Eugenio November 2017
        if ( exist( 'optimizer', 'var' ) == 1 )
            try
                kvlClear( optimizer )
                kvlClear( calculator )
                clear optimizer
                clear calculator
            end
        end
        
    end % End of loop over time points
    
    % Update positions
    for t=1:nTP
        meshes{t} = kvlGetMesh( meshCollections{t}, 0 );
        subjectTPpositions{t} = kvlGetMeshNodePositions( meshes{t} );
    end
    
    if globalIt >= MAX_GLOBAL_ITS
        globalReady=1;
    end
    
end  % end of loop over global iterations

%
if  ~isdeployed && DEBUG>0
    figure
    for t=1:nTP
        subplot( 2, nTP, t )
        showImage( kvlColorCodeProbabilityImages( kvlRasterizeAtlasMesh( meshes{t}, imageSize )  ) );
        title(['Time point ' num2str(t) ': priors'])
        subplot( 2, nTP, t+nTP )
        showImage( imageBuffers{t} )
        title(['Time point ' num2str(t) ': image'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final part: compute volumes and segmentations at each time point %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for t=1:nTP
    
    % Undo the collapsing of several structures into "super"-structures
    kvlSetAlphasInMeshNodes( meshes{t}, originalAlphas )
    numberOfClasses = size( originalAlphas, 2 );
    
    
    % Eugenio July 2017 : this type of call is exactly what we're trying to
    % avoid....
    data = double( reshape( imageBuffers{t}, [ prod( imageSize ) 1 ] ) ); % Easier to work with vector notation in the computations
    data = data( maskIndices );
    
    posteriors=zeros([length(maskIndices),numberOfClasses],'single');
    for classNumber = 1 : numberOfClasses
        prior = kvlRasterizeAtlasMesh( meshes{t}, imageSize, classNumber-1 );
        mu = means{t}( reducingLookupTable( classNumber ) );
        variance = variances{t}( reducingLookupTable( classNumber ) );
        posteriors( :, classNumber ) = ( exp( -( data - mu ).^2 / 2 / variance ) ...
            .* (double(prior(maskIndices))/65535) ) / sqrt( 2 * pi * variance);
    end
    normalizer = sum( posteriors, 2 ) + eps;
    posteriors = bsxfun(@rdivide,posteriors, normalizer);
    posteriors = uint16( round( posteriors * 65535 ) );
    
    %     % Get the priors as dictated by the current mesh position
    %     data = double( reshape( imageBuffers{t}, [ prod( imageSize ) 1 ] ) ); % Easier to work with vector notation in the computations
    %     priors = kvlRasterizeAtlasMesh( meshes{t}, imageSize );
    %     priors = reshape( priors, [ prod( imageSize ) numberOfClasses ] );  % Easier to work with vector notation in the computations
    %
    %     % Ignore everything that's has zero intensity
    %     priors = priors( maskIndices, : );
    %     data = data( maskIndices );
    %
    %     % Calculate the posteriors
    %     posteriors = zeros( size( priors ), 'double' );
    %     for classNumber = 1 : numberOfClasses
    %         % Get the parameters from the correct Gaussian
    %         mu = means{t}( reducingLookupTable( classNumber ) );
    %         variance = variances{t}( reducingLookupTable( classNumber ) );
    %         prior = single( priors( :, classNumber ) ) / 65535;
    %
    %         posteriors( :, classNumber ) = ( exp( -( data - mu ).^2 / 2 / variance ) .* prior ) ...
    %             / sqrt( 2 * pi * variance);
    %     end
    %     normalizer = sum( posteriors, 2 ) + eps;
    %     posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfClasses ] );
    %     posteriors = uint16( round( posteriors * 65535 ) );
    
    % Display the posteriors
    if  ~isdeployed && DEBUG>0
        figure
        subplot( 1, 2, 1 )
        showImage( imageBuffers{t} ); title(['Time point ' num2str(t) ': image'])
        subplot( 1, 2, 2 )
        tmp = zeros( prod(imageSize), numberOfClasses, 'uint16' );
        tmp( maskIndices, : ) = priors;
        tmp = reshape( tmp, [ imageSize numberOfClasses ] );
        showImage( kvlColorCodeProbabilityImages( tmp, colors ) ); title(['Time point ' num2str(t) ': warped priors'])
        
        figure
        for classNumber = 1 : numberOfClasses
            
            % Let's create a color image overlaying the posterior on top of gray scale
            overlayColor = [ 0 0 1 ]; % RGB value
            overlay = double( exp( imageBuffers{t} / 1000 )  );
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
            mu = means{t}( reducingLookupTable( classNumber ) );
            variance = variances{t}( reducingLookupTable( classNumber ) );
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
    kvlWriteMeshCollection( meshCollections{t}, ['warpedMesh_tp_' num2str(t) '.txt'] ); % The warped mesh will have index 0, and will be
    
    transformMatrix = double(kvlGetTransformMatrix( transform ));
    inverseTransform = kvlCreateTransform( inv( transformMatrix ) );
    kvlTransformMeshCollection( meshCollections{t}, inverseTransform );
    kvlWriteMeshCollection( meshCollections{t}, ['warpedMeshNoAffine_tp_' num2str(t) '.txt'] );
    
    kvlWriteImage( images{t}, ['image_tp_' num2str(t) '.mgz'] );
    if t==1, system(['cp ' compressionLookupTableFileName ' .']); end
    
    
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
    % We only do this for the first time point: the others have the same shift
    if t==1
        disp('Computing shift introduced by kvlRead')
        system([FSpath '/mri_convert asegMod_tp_1.mgz asmr1.mgz -rl T1resampled_1.mgz -rt nearest -odt float >/dev/null']);
        [ asmr, tasmr ] = kvlReadCroppedImage( 'asmr1.mgz', boundingFileName );
        asmrB=kvlGetImageBuffer(asmr);
        kvlWriteImage( asmr, 'asmr2.mgz' );
        tmp1=myMRIread('asmr1.mgz',0,tempdir);
        tmp2=myMRIread('asmr2.mgz',0,tempdir);
        tmp3=tmp1; tmp3.vol=tmp2.vol;
        myMRIwrite(tmp3,'asmr3.mgz','float',tempdir);
        [I,J,K]=ind2sub(size(tmp1.vol),find(tmp1.vol==17));
        Ic1=mean(I); Jc1=mean(J); Kc1=mean(K);
        [I,J,K]=ind2sub(size(tmp2.vol),find(tmp2.vol==17));
        Ic2=mean(I); Jc2=mean(J); Kc2=mean(K);
        shift=round([Ic2-Ic1,Jc2-Jc1,Kc2-Kc1]);
        shiftNeg=-shift; shiftNeg(shiftNeg<0)=0;
        shiftPos=shift; shiftPos(shiftPos<0)=0;
    end
    
    % reorient image.mgz
    aux=zeros(size(tmp2.vol)+shiftNeg);
    aux2=myMRIread(['image_tp_' num2str(t) '.mgz'],0,tempdir);
    aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=aux2.vol;
    aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
    tmp3.vol=aux;
    myMRIwrite(tmp3,['image_tp_' num2str(t) '.mgz'],'float',tempdir);
    
    
    % July 2017: make this memory efficient, and add new substructures. We also
    % compute the segmentation along the way. Also, we do hippo and amygdala
    % together. Ah! And we also do volumes for head and body
    fidHP=fopen([tempdir '/volumesHippo_tp_' num2str(t) '.txt'],'w');
    strOfInterestHP={'subiculum-body','subiculum-head','Hippocampal_tail','molecular_layer_HP-body','molecular_layer_HP-head','hippocampal-fissure','GC-ML-DG-body','GC-ML-DG-head','CA4-body','CA4-head','presubiculum-body','presubiculum-head','CA1-body','CA1-head','parasubiculum','fimbria','CA3-body','CA3-head','HATA'};
    totVolHP=0;
    foundHP=zeros([1,numberOfClasses]);
    HPbodyList={'subiculum-body','CA1-body','presubiculum-body','molecular_layer_HP-body','CA3-body','GC-ML-DG-body','CA4-body','fimbria'};
    HPheadList={'subiculum-head','presubiculum-head','CA1-head','parasubiculum','molecular_layer_HP-head','GC-ML-DG-head','CA4-head','CA3-head','HATA'};
    totVolHPbody=0;
    totVolHPhead=0;
    
    fidAM=fopen([tempdir '/volumesAmygdala_tp_' num2str(t) '.txt'],'w');
    strOfInterestAM={'Left-Amygdala','Lateral-nucleus','Paralaminar-nucleus',...
        'Basal-nucleus','Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
        'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus','Anterior-amygdaloid-area-AAA'};
    totVolAM=0;
    foundAM=zeros([1,numberOfClasses]);
    
    for i=1:numberOfClasses
        
        if i==1
            sillyAlphas=zeros([size(originalAlphas,1),2],'single');
            sillyAlphas(:,1)=originalAlphas(:,1);
            sillyAlphas(:,2)=1-sillyAlphas(:,1);
            kvlSetAlphasInMeshNodes( meshes{t}, sillyAlphas )
            post = kvlRasterizeAtlasMesh( meshes{t}, imageSize);
            post=post(:,:,:,1);
            kvlSetAlphasInMeshNodes( meshes{t}, originalAlphas );
        else
            post=kvlRasterizeAtlasMesh( meshes{t}, imageSize , i-1);
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
        
        foundAM(i)=0;
        foundHP(i)=0;
        
        name=names(i,:);
        name=lower(name(name~=' '));
        
        for j=1:length(strOfInterestHP)
            if strcmp(name,lower(strOfInterestHP{j}))>0
                foundHP(i)=j;
            end
        end
        
        for j=1:length(strOfInterestAM)
            if strcmp(name,lower(strOfInterestAM{j}))>0
                foundAM(i)=j;
            end
        end
        
        if foundHP(i)>0 || foundAM(i)>0
            % Fix by Eugenio in February 2022
            % vol=resolution^3*(sum(double(post(:))/65535));
            vol=resolution^3*(sum(double(post(:))/65535)) / volumeFactors(t);
            
            if foundHP(i)>0
                str=strOfInterestHP{foundHP(i)};
                fprintf(fidHP,'%s %f\n',str,vol);
                if isempty(strfind(lower(names(i,:)),'hippocampal-fissure'))  % don't count the fissure towards the total volume
                    totVolHP=totVolHP+vol;
                end
            else
                str=strOfInterestAM{foundAM(i)};
                fprintf(fidAM,'%s %f\n',str,vol);
                totVolAM=totVolAM+vol;
            end
            
            for j=1:length(HPbodyList)
                if strcmp(name,lower(HPbodyList{j}))>0
                    totVolHPbody=totVolHPbody+vol;
                end
            end
            for j=1:length(HPheadList)
                if strcmp(name,lower(HPheadList{j}))>0
                    totVolHPhead=totVolHPhead+vol;
                end
            end
            
            if WRITE_POSTERIORS>0
                kk1=double(post)/65535;
                kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
                aux=zeros(size(tmp2.vol)+shiftNeg);
                aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
                aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
                tmp3.vol=aux;
                myMRIwrite(tmp3,['posterior_' side '_' strtrim(lower(names(i,:))) '_tp' num2str(t) '_T1_' suffix '.mgz'],'float',tempdir);
            end
            
        end
    end
    if sum(foundHP>0)>1
        fprintf(fidHP,'Whole_hippocampal_body %f\n',totVolHPbody);
        fprintf(fidHP,'Whole_hippocampal_head %f\n',totVolHPhead);
        fprintf(fidHP,'Whole_hippocampus %f\n',totVolHP);
        fclose(fidHP);
    else
        fclose(fidHP);
        delete([tempdir '/volumesHippo.txt']);
    end
    if sum(foundAM>0)>1
        fprintf(fidAM,'Whole_amygdala %f\n',totVolAM);
        fclose(fidAM);
    else
        fclose(fidAM);
        delete([tempdir '/volumesAmygdala.txt']);
    end
    
    
    %     % Compute posteriors and volumes - and write volumes to text files
    %     priorsFull = kvlRasterizeAtlasMesh( meshes{t}, imageSize );
    %     posteriorsFull=priorsFull;
    %
    %     fid=fopen([tempdir '/volumesHippo_tp_' num2str(t) '.txt'],'w');
    %     % strOfInterest={'alveus','subiculum','Hippocampal_tail','molecular_layer_HP','hippocampal-fissure','GC-ML-DG','CA4','presubiculum','CA1','parasubiculum','fimbria','CA3','HATA'};
    %     % no alveus
    %     strOfInterest={'subiculum','Hippocampal_tail','molecular_layer_HP','hippocampal-fissure','GC-ML-DG','CA4','presubiculum','CA1','parasubiculum','fimbria','CA3','HATA'};
    %     totVol=0;
    %     found=zeros(1,size(priorsFull,4));
    %     for i=1:size(priorsFull,4)
    %         tmp=posteriorsFull(:,:,:,i);
    %         tmp(maskIndices)=posteriors(:,i);
    %         posteriorsFull(:,:,:,i)=tmp;
    %         found(i)=0;
    %         str=[];
    %         vol=0;
    %         name=names(i,:);
    %         name=lower(name(name~=' '));
    %         for j=1:length(strOfInterest)
    %             if strcmp(name,lower(strOfInterest{j}))>0
    %                 found(i)=j;
    %             end
    %         end
    %         if found(i)>0
    %             str=strOfInterest{found(i)};
    %             vol=resolution^3*(sum(sum(sum(double(posteriorsFull(:,:,:,i))/65535))))/volumeFactors(t);
    %             fprintf(fid,'%s %f\n',str,vol);
    %             if isempty(strfind(lower(names(i,:)),'hippocampal-fissure'))  % don't count the fissure towards the total volume
    %                 totVol=totVol+vol;
    %
    %                 if WRITE_POSTERIORS>0
    %                     kk1=double(posteriorsFull(:,:,:,i))/65535;
    %                     kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
    %                     aux=zeros(size(tmp2.vol)+shiftNeg);
    %                     aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
    %                     aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
    %                     tmp3.vol=aux;
    %                     myMRIwrite(tmp3,['posterior_' side '_' strtrim(lower(names(i,:))) '_tp' num2str(t) '_T1_' suffix '.mgz'],'float',tempdir);
    %                 end
    %
    %             end
    %         end
    %     end
    %     fprintf(fid,'Whole_hippocampus %f\n',totVol);
    %     fclose(fid);
    %
    %
    %     fid=fopen([tempdir '/volumesAmygdala_tp_' num2str(t) '.txt'],'w');
    %     strOfInterest={'Left-Amygdala','Lateral-nucleus','Paralaminar-nucleus',...
    %         'Basal-nucleus','Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
    %         'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus','Anterior-amygdaloid-area-AAA'};
    %     totVol=0;
    %     found=zeros(1,size(priorsFull,4));
    %     for i=1:size(priorsFull,4)
    %         tmp=posteriorsFull(:,:,:,i);
    %         tmp(maskIndices)=posteriors(:,i);
    %         posteriorsFull(:,:,:,i)=tmp;
    %         found(i)=0;
    %         str=[];
    %         vol=0;
    %         name=names(i,:);
    %         name=lower(name(name~=' '));
    %         for j=1:length(strOfInterest)
    %             if strcmp(name,lower(strOfInterest{j}))>0
    %                 found(i)=j;
    %             end
    %         end
    %         if found(i)>0
    %             str=strOfInterest{found(i)};
    %             vol=resolution^3*(sum(sum(sum(double(posteriorsFull(:,:,:,i))/65535))))/volumeFactors(t);
    %             fprintf(fid,'%s %f\n',str,vol);
    %             totVol=totVol+vol;
    %
    %             if WRITE_POSTERIORS>0
    %                     kk1=double(posteriorsFull(:,:,:,i))/65535;
    %                     kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
    %                     aux=zeros(size(tmp2.vol)+shiftNeg);
    %                     aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
    %                     aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
    %                     tmp3.vol=aux;
    %                     myMRIwrite(tmp3,['posterior_' side '_' strtrim(lower(names(i,:))) '_tp' num2str(t) '_T1_' suffix '.mgz'],'float',tempdir);
    %             end
    %
    %         end
    %     end
    %     if sum(found>0)>1
    %         fprintf(fid,'Whole_amygdala %f\n',totVol);
    %         fclose(fid);
    %     else
    %         fclose(fid);
    %         delete([tempdir '/volumesAmygdala_tp_' num2str(t) '.txt']);
    %     end
    
    
    % MAP estimates
    
    % Eugenio July 2011
    % [~,inds]=max(posteriorsFull,[],4);
    inds=L;
    
    
    kk1=FreeSurferLabels(inds);
    kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
    aux=zeros(size(tmp2.vol)+shiftNeg);
    aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
    aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
    tmp3.vol=aux;
    myMRIwrite(tmp3,['discreteLabels_all_tp_' num2str(t) '.mgz'],'float',tempdir);
    tmp3.vol(tmp3.vol<200)=0;
    
    % Eugenio July 2017
    % tmp3.vol(tmp3.vol>226 & tmp3.vol<7000)=0;
    tmp3.vol(tmp3.vol>246 & tmp3.vol<7000)=0;
    
    tmp3.vol(tmp3.vol==201)=0; % alveus
    tmp3Mask=getLargestCC(tmp3.vol>0);
    tmp3.vol(~tmp3Mask)=0;
    myMRIwrite(tmp3,['discreteLabels_tp_' num2str(t) '.mgz'],'float',tempdir);
    
    
    % Eugenio July 2011
    % Write merged versions to disk as well
    % First: tail / body /head
    
    HippoBodyLabel=231;
    HippoHeadLabel=232;
    
    tmp4=tmp3;
    for c=1:numberOfClasses
        name=names(c,:);
        name=lower(name(name~=' '));
        
        found=0;
        for j=1:length(HPbodyList)
            if strcmp(name,lower(HPbodyList{j}))>0
                found=1;
            end
        end
        if found==1
            tmp4.vol(tmp4.vol==FreeSurferLabels(c))=HippoBodyLabel;
        end
        
        found=0;
        for j=1:length(HPheadList)
            if strcmp(name,lower(HPheadList{j}))>0
                found=1;
            end
        end
        if found==1
            tmp4.vol(tmp4.vol==FreeSurferLabels(c))=HippoHeadLabel;
        end
    end
    tmp4.vol(tmp4.vol==215)=0; % kill the fissure
    myMRIwrite(tmp4,['discreteLabelsWholeBodyHead_tp_' num2str(t) '.mgz'],'float',tempdir);
    
    % Second: head+body of each subfield
    tmp4=tmp3;
    tmp4.vol(tmp3.vol==233 | tmp3.vol==234)=204; % presubiculum
    tmp4.vol(tmp3.vol==235 | tmp3.vol==236)=205; % subiculum
    tmp4.vol(tmp3.vol==237 | tmp3.vol==238)=206; % CA1
    tmp4.vol(tmp3.vol==239 | tmp3.vol==240)=208; % CA3
    tmp4.vol(tmp3.vol==241 | tmp3.vol==242)=209; % CA4
    tmp4.vol(tmp3.vol==243 | tmp3.vol==244)=210; % GC-DG
    tmp4.vol(tmp3.vol==245 | tmp3.vol==246)=214; % ML
    
    myMRIwrite(tmp4,['discreteLabelsMergedBodyHead_tp_' num2str(t) '.mgz'],'float',tempdir);
    
    % Third: same as above, but getting rid of internal labels
    tmp4.vol(tmp4.vol==210)=209;  % GC-DG -> CA4
    
    % Molecular layer: replace by nearest label that is not background or
    % fissure
    [~,cropping]=cropLabelVol(tmp4.vol==214,2);
    VOL=applyCropping(tmp4.vol,cropping);
    llist=unique(VOL);
    llist=llist(llist~=0 & llist~=215 & llist~=214);
    mask=VOL==214;
    for l=1:length(llist)
        label=llist(l);
        dmap=bwdist(VOL==label);
        if l==1
            mini=dmap(mask);
            seg=label;
        else
            dist=dmap(mask);
            m=dist<mini;
            mini(m)=dist(m);
            seg(m)=label;
        end
    end
    VOL(mask)=seg;
    tmp4.vol(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6))=VOL;
    myMRIwrite(tmp4,['discreteLabelsMergedBodyHeadNoMLorGCDG_tp_' num2str(t) '.mgz'],'float',tempdir);
    
    
    % Eugenio July 2017
    % I disabled this for now ...
    if  MRFconstant>0
        
        % Eugenio July 2017
        error('MRF smoothing disabled for now');
        
        %         EPS=1e-12;
        %         [~,inds]=max(posteriorsFull,[],4);
        %         tmp=FreeSurferLabels(inds);
        %         kk=zeros(size(tmp)); kk(maskIndices)=1; tmp=tmp.*kk;
        %         tmp(tmp<200)=0; tmp(tmp>226 & tmp<7000)=0;
        %         [~,cropping]=cropLabelVol(tmp);
        %         Ct=zeros([cropping(4)-cropping(1)+1,cropping(5)-cropping(2)+1,cropping(6)-cropping(3)+1,numberOfClasses]);
        %         for c=1:numberOfClasses
        %             Ct(:,:,:,c)=-log(EPS+double(posteriorsFull(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6),c))/65535);
        %         end
        %         factor=-256/log(EPS);
        %         Ct=int32(round(Ct*factor));
        %         unaryTermWeight=int32(round(MRFconstant*factor));
        %
        %         siz=[size(Ct,1) size(Ct,2) size(Ct,3)];
        %         h = GCO_Create(prod(siz),numberOfClasses);
        %         DC = zeros([numberOfClasses,prod(siz)],'int32');
        %         for c=1:numberOfClasses
        %             aux=Ct(:,:,:,c);
        %             DC(c,:)=aux(:);
        %         end
        %         GCO_SetDataCost(h,DC);
        %         aux=int32(double(unaryTermWeight)*(ones(numberOfClasses)-eye(numberOfClasses)));
        %         GCO_SetSmoothCost(h,aux);
        %
        %         row=zeros([prod(siz)*3,1]);
        %         col=zeros([prod(siz)*3,1]);
        %         t=1;
        %
        %         Ifrom=1:siz(1)-1;
        %         Ito=2:siz(1);
        %         inc=length(Ito);
        %         for j=1:siz(2)
        %             J=j*ones(size(Ifrom));
        %             for k=1:siz(3)
        %                 K=k*ones(size(Ifrom));
        %                 row(t:t+inc-1)=sub2ind(siz,Ifrom,J,K);
        %                 col(t:t+inc-1)=sub2ind(siz,Ito,J,K);
        %                 t=t+inc;
        %             end
        %         end
        %
        %         Jfrom=1:siz(2)-1;
        %         Jto=2:siz(2);
        %         inc=length(Jto);
        %         for i=1:siz(1)
        %             I=i*ones(size(Jfrom));
        %             for k=1:siz(3)
        %                 K=k*ones(size(Jfrom));
        %                 row(t:t+inc-1)=sub2ind(siz,I,Jfrom,K);
        %                 col(t:t+inc-1)=sub2ind(siz,I,Jto,K);
        %                 t=t+inc;
        %             end
        %         end
        %
        %         Kfrom=1:siz(3)-1;
        %         Kto=2:siz(3);
        %         inc=length(Kto);
        %         for i=1:siz(1)
        %             I=i*ones(size(Kfrom));
        %             for j=1:siz(2)
        %                 J=j*ones(size(Kfrom));
        %                 row(t:t+inc-1)=sub2ind(siz,I,J,Kfrom);
        %                 col(t:t+inc-1)=sub2ind(siz,I,J,Kto);
        %                 t=t+inc;
        %             end
        %         end
        %
        %         row=row(1:t-1);
        %         col=col(1:t-1);
        %
        %         NEIGH=sparse(row,col,ones(size(row)),prod(siz),prod(siz));
        %         GCO_SetNeighbors(h,NEIGH);
        %
        %
        %         GCO_Expansion(h);      % Compute optimal labeling via alpha-expansion
        %         ind=reshape(GCO_GetLabeling(h),siz);
        %
        %         SEG=FreeSurferLabels(ind);
        %         SEG(SEG>226 & SEG<7000)=0; SEG(SEG<200)=0;  SEG(SEG==201)=0;
        %
        %         data=zeros(size(inds));
        %         data(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6))=SEG;
        %         aux=zeros(size(tmp2.vol)+shiftNeg);
        %         aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(data,[2 1 3]);
        %         aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
        %         tmp3.vol=aux;
        %         tmp3Mask=getLargestCC(tmp3.vol>0);
        %         tmp3.vol(~tmp3Mask)=0;
        %         myMRIwrite(tmp3,['discreteLabels_MRF_tp_ ' num2str(t) '.mgz'],'float',tempdir);
        %
    end
    
    
    % Convert to 1 mm FreeSurfer Space
    % Eugenio July 2017: add new segmentation maps
    % Eugenio November 2017: alternative resampling strategy
    if SMOOTH_LABEL_RESAMPLE>0
        
        refFile=[subjectDir '/' subjectTPs{t} '/mri/norm.mgz'];
        tfile=['tp_' num2str(t) '_to_base.lta'];
        
        applyLTAsmoothLabels(['discreteLabels_tp_' num2str(t) '.mgz'],tfile,['discreteLabelsResampledT1_tp_' num2str(t) '.mgz'],refFile,1,FSpath,tempdir);
        applyLTAsmoothLabels(['discreteLabelsWholeBodyHead_tp_' num2str(t) '.mgz'],tfile,['discreteLabelsWholeBodyHeadResampledT1_tp_' num2str(t) '.mgz'],refFile,1,FSpath,tempdir);
        applyLTAsmoothLabels(['discreteLabelsMergedBodyHead_tp_' num2str(t) '.mgz'],tfile,['discreteLabelsMergedBodyHeadResampledT1_tp_' num2str(t) '.mgz'],refFile,1,FSpath,tempdir);
        applyLTAsmoothLabels(['discreteLabelsMergedBodyHeadNoMLorGCDG_tp_' num2str(t) '.mgz'],tfile,['discreteLabelsMergedBodyHeadNoMLorGCDGResampledT1_tp_' num2str(t) '.mgz'],refFile,1,FSpath,tempdir);
        
    else
        
        system([FSpath '/mri_convert  discreteLabels_tp_' num2str(t) '.mgz  discreteLabelsResampledT1_tp_' num2str(t) '.mgz -rt nearest -odt float ' ...
            ' -rl ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz -ait tp_' num2str(t) '_to_base.lta']);
        system([FSpath '/mri_convert  discreteLabelsWholeBodyHead_tp_' num2str(t) '.mgz  discreteLabelsWholeBodyHeadResampledT1_tp_' num2str(t) '.mgz -rt nearest -odt float ' ...
            ' -rl ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz -ait tp_' num2str(t) '_to_base.lta']);
        system([FSpath '/mri_convert  discreteLabelsMergedBodyHead_tp_' num2str(t) '.mgz  discreteLabelsMergedBodyHeadResampledT1_tp_' num2str(t) '.mgz -rt nearest -odt float ' ...
            ' -rl ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz -ait tp_' num2str(t) '_to_base.lta']);
        system([FSpath '/mri_convert  discreteLabelsMergedBodyHeadNoMLorGCDG_tp_' num2str(t) '.mgz  discreteLabelsMergedBodyHeadNoMLorGCDGResampledT1_tp_' num2str(t) '.mgz -rt nearest -odt float ' ...
            ' -rl ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz -ait tp_' num2str(t) '_to_base.lta']);
    end
    
end  % end of loop over time points


% We need to warp back to original space
for t=1:nTP
    % Eugenio July 2017: go over all segmentations
    system([FSpath '/mri_vol2vol --mov discreteLabels_tp_' num2str(t) '.mgz  --o discreteLabels_tp_' num2str(t) '_origSpace.mgz --no-resample ' ...
        ' --targ ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz --lta-inv tp_' num2str(t) '_to_base.lta']);
    system([FSpath '/mri_vol2vol --mov discreteLabelsWholeBodyHead_tp_' num2str(t) '.mgz  --o discreteLabelsWholeBodyHead_tp_' num2str(t) '_origSpace.mgz --no-resample ' ...
        ' --targ ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz --lta-inv tp_' num2str(t) '_to_base.lta']);
    system([FSpath '/mri_vol2vol --mov discreteLabelsMergedBodyHead_tp_' num2str(t) '.mgz  --o discreteLabelsMergedBodyHead_tp_' num2str(t) '_origSpace.mgz --no-resample ' ...
        ' --targ ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz --lta-inv tp_' num2str(t) '_to_base.lta']);
    system([FSpath '/mri_vol2vol --mov discreteLabelsMergedBodyHeadNoMLorGCDG_tp_' num2str(t) '.mgz  --o discreteLabelsMergedBodyHeadNoMLorGCDG_tp_' num2str(t) '_origSpace.mgz --no-resample ' ...
        ' --targ ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz --lta-inv tp_' num2str(t) '_to_base.lta']);
    
    if WRITE_POSTERIORS>0
        d=dir(['posterior_' side '_*_tp' num2str(t) '_T1_' suffix '.mgz']);
        for i=1:length(d)
            if isempty(strfind(d(i).name,'left-amygdala'))
                system([FSpath '/mri_vol2vol --mov ' d(i).name '  --o ' d(i).name  ' --no-resample ' ...
                    ' --targ ' subjectDir '/' subjectTPs{t} '/mri/norm.mgz --lta-inv tp_' num2str(t) '_to_base.lta']);
                
            end
        end
    end
    
end

% Finally, move results to  MRI directory
for t=1:nTP
    
    % Eugenio July 2017: add new segmentation maps, and simplified code
    % (left/right)
    
    system(['mv discreteLabels_tp_' num2str(t) '_origSpace.mgz ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoAmygLabels-T1.' suffix '.mgz']);
    system(['mv discreteLabelsResampledT1_tp_' num2str(t) '.mgz ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoAmygLabels-T1.' suffix '.FSvoxelSpace.mgz']);
    system(['mv discreteLabelsWholeBodyHead_tp_' num2str(t) '_origSpace.mgz ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoAmygLabels-T1.' suffix '.HBT.mgz']);
    system(['mv discreteLabelsWholeBodyHeadResampledT1_tp_' num2str(t) '.mgz ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoAmygLabels-T1.' suffix '.HBT.FSvoxelSpace.mgz']);
    system(['mv discreteLabelsMergedBodyHead_tp_' num2str(t) '_origSpace.mgz ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoAmygLabels-T1.' suffix '.FS60.mgz']);
    system(['mv discreteLabelsMergedBodyHeadResampledT1_tp_' num2str(t) '.mgz ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoAmygLabels-T1.' suffix '.FS60.FSvoxelSpace.mgz']);
    system(['mv discreteLabelsMergedBodyHeadNoMLorGCDG_tp_' num2str(t) '_origSpace.mgz ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoAmygLabels-T1.' suffix '.CA.mgz']);
    system(['mv discreteLabelsMergedBodyHeadNoMLorGCDGResampledT1_tp_' num2str(t) '.mgz ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoAmygLabels-T1.' suffix '.CA.FSvoxelSpace.mgz']);
    
    system(['mv volumesHippo_tp_' num2str(t) '.txt ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.hippoSfVolumes-T1.' suffix '.txt']);
    system(['mv volumesAmygdala_tp_' num2str(t) '.txt ' subjectDir '/' subjectTPs{t} '/mri/' side(1) 'h.amygNucVolumes-T1.' suffix '.txt']);
    
    
    if WRITE_POSTERIORS>0
        d=dir(['posterior_' side '_*_tp' num2str(t) '_T1_' suffix '.mgz']);
        for i=1:length(d)
            if isempty(strfind(d(i).name,'left-amygdala'))
                str=['_tp' num2str(t)];
                p=strfind(d(i).name,str);
                p=p(1);
                tname=[d(i).name(1:p-1) d(i).name(p+length(str):end)];
                system(['mv  ' d(i).name ' ' subjectDir '/' subjectTPs{t} '/mri/' tname]);
            end
        end
    end
    
    if WRITE_MESHES>0
        system(['mv warpedMesh_tp_' num2str(t) '.txt.gz   ' subjectDir '/'  subjectTPs{t} '/mri/' side(1) 'h.hippoAmygMesh-T1.' suffix '.txt.gz']);
        system(['mv warpedMeshNoAffine_tp_' num2str(t) '.txt.gz   ' subjectDir '/'  subjectTPs{t} '/mri/' side(1) 'h.hippoAmygMeshAtlasSpace-T1.' suffix '.txt.gz']);
        system(['mv image_tp_' num2str(t) '.mgz   ' subjectDir '/'  subjectTPs{t} '/mri/' side(1) 'h.imageForMesh-T1.' suffix '.mgz']);
        fid=fopen([subjectDir '/'  subjectTPs{t} '/mri/' side(1) 'h.affineTransformMesh-T1.' suffix '.txt'],'w');
        fprintf(fid,'%f %f %f %f \n',transformMatrix(1,1),transformMatrix(1,2),transformMatrix(1,3),transformMatrix(1,4));
        fprintf(fid,'%f %f %f %f \n',transformMatrix(2,1),transformMatrix(2,2),transformMatrix(2,3),transformMatrix(2,4));
        fprintf(fid,'%f %f %f %f \n',transformMatrix(3,1),transformMatrix(3,2),transformMatrix(3,3),transformMatrix(3,4));
        fprintf(fid,'%f %f %f %f \n',0,0,0,1);
        fclose(fid);
    end
    
    
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
