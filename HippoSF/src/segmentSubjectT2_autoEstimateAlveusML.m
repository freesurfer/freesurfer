
% Segments the subfields from a T2 scan; requires the FS recon from the original MPRAGE

%
% segmentSubjectT2(subjectName,subjectDir,T2volumeFileName,resolution,atlasMeshFileName,atlasDumpFileName,compressionLUTfileName,K,side,optimizerType,suffix,MRFconstant,ByPassBF,UseWholeBrainInHyperPar)
%
% - subjectName: FreeSurfer subject name
% - subjectDir: FreeSurfer subject directory
% - T2volumeFileName: T2 image
% - resolution: voxel size at which we want to work (in mm).
% - atlasMeshFileName: the atlas to segment the datas
% - atlasDumpFileName: corresponding imageDump.mgz (name *must* be imageDump.mgz)
% - compressionLUTfileName: corresponding compressionLUT.txt
% - K: stiffness of the mesh in the segmentation (around 0.05)
% - side: 'left' or 'right'
% - optimizerType: 'FixedStepGradientDescent','GradientDescent','ConjugateGradient','L-BFGS'
% - suffix: for output directory, e.g. 'T1based_GGAWLnoSimil','v10',...
% - suffixUser: for the user to identify outputs from different T2s
% - FSpath: path to FreeSurfer executables
% - MRFconstant (optional): make it >0 for MRF cleanup (5 is reasonable, larger is smoother)
%           It does NOT affect volumes, which are computed from soft posteriors anyway


function segmentSubjectT2_autoEstimateAlveusML(subjectName,subjectDir,T2volumeFileName,resolution,atlasMeshFileName,atlasDumpFileName,compressionLUTfileName,K,side,optimizerType,suffix,suffixUser,FSpath,MRFconstant,ByPassBF,UseWholeBrainInHyperPar)


% Eugenio November 2017: added option to write meshes and smoother resampling
DEBUG=0;
FAST=0; % set it to one to optimize just a bit (go through code fast)
WRITE_POSTERIORS=0;
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


% sanity check
if exist('MRFconstant','var')==0
    MRFconstant=0;
end

if nargin<12
    error('Not enough input arguments');
elseif strcmp(side,'left')==0 && strcmp(side,'right')==0
    error('Side must be ''left'' or ''right''');
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
elseif ~isdeployed && (~isnumeric(MRFconstant))
    error('MRFconstant must be numeric');
end


% Constants
HippoLabelLeft=17;
HippoLabelRight=53;

% In case we compiled it...
if isdeployed
    K=str2double(K);
    resolution=str2double(resolution);
    if exist('MRFconstant','var')==0
        MRFconstant=0;
    else
        MRFconstant=str2double(MRFconstant);
    end
    if exist('ByPassBF','var')==0
        ByPassBF=0;
    else
        ByPassBF=str2double(ByPassBF);
    end
    if exist('UseWholeBrainInHyperPar','var')==0
        UseWholeBrainInHyperPar=0;
    else
        UseWholeBrainInHyperPar=str2double(UseWholeBrainInHyperPar);
    end
else
    if exist('MRFconstant','var')==0
        MRFconstant=0;
    end
    if exist('ByPassBF','var')==0
        ByPassBF=0;
    end
    if exist('UseWholeBrainInHyperPar','var')==0
        UseWholeBrainInHyperPar=0;
    end
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
tempdir=[subjectDir '/' subjectName '/tmp/hippoSF_T2_' suffix '_' suffixUser '_' side '/'];
aux=getenv('USE_SCRATCH');
if ~isempty(aux)
    if str2double(aux)>0
        if exist('/scratch/','dir')>0
            tempdir=['/scratch/' subjectName '_hippoSF_T2_' suffix '_' suffixUser '_' side '/'];
        end
    end
end
if exist(tempdir,'dir')==0
    mkdir(tempdir);
end

T2volumeFileName=getFullPath(T2volumeFileName); % Eugenio November 2017: before we cd

cd(tempdir);

f=find(T2volumeFileName=='/');
if isempty(f)
    T2name=T2volumeFileName;
else
    T2name=T2volumeFileName(f(end)+1:end);
end
T2transform=[subjectDir '/' subjectName '/mri/transforms/T1_to_' suffixUser '.' suffix '.lta'];
T2transformInfo=[subjectDir '/' subjectName '/mri/transforms/T1_to_' suffixUser '.' suffix '.info'];


% First thing: are T1 and T2 registered?
if exist(T2transform,'file')>0 && exist(T2transformInfo,'file')>0
    disp('Found transform from T1 to additional volume; no need to register')
    fid=fopen(T2transformInfo);
    tline=fgetl(fid); tline=fgetl(fid);
    tline=fgetl(fid); tline=fgetl(fid);
    if strcmp(tline,T2name)==0
        error('Found T1-T2 transform does not correspond to T2 input. If you are using a different T2, please change the ID of the analysis');
    end
else
    
    fid=fopen(T2transformInfo,'w');
    fprintf(fid,'The transform:\n');
    fprintf(fid,'%s\n',T2transform);
    fprintf(fid,'was created with additional image file:\n');
    fprintf(fid,'%s\n',T2name);
    fclose(fid);
    
    disp('Registering norm.mgz to additional volume');
    
    WHBRmaskFile=[tempdir '/wholeBrainMask.mgz'];
    
    % Read in ASEG and NU and create whole brain mask
    % ASEG=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
    NORM=myMRIread([subjectDir '/' subjectName '/mri/norm.mgz'],0,tempdir);
    WHBRmask=NORM;
    WHBRmask.vol=NORM.vol>0;
    myMRIwrite(WHBRmask,WHBRmaskFile,'float',tempdir);
    % Register norm.mgz to T2!
    cmd=[FSpath '/mri_robust_register --mov ' subjectDir '/' subjectName '/mri/norm.mgz --maskmov ' WHBRmaskFile ' --dst ' T2volumeFileName ' --lta ' T2transform ' --noinit --cost NMI -nosym'];
    system([cmd ' >/dev/null']);
    delete(WHBRmaskFile);
    
    fid=fopen(T2transformInfo,'w');
    fprintf(fid,'The transform:\n');
    fprintf(fid,'%s\n',T2transform);
    fprintf(fid,'was created with additional image file:\n');
    fprintf(fid,'%s\n',T2name);
    fclose(fid);
    
    % Create animated gif for QC of registration
    auxT1=myMRIread([subjectDir '/' subjectName '/mri/nu.mgz'],0,tempdir);
    auxSEG=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
    cmd=[FSpath '/mri_convert ' T2volumeFileName ' ' tempdir '/T2resampledLikeT1.mgz -rl ' subjectDir '/' subjectName '/mri/nu.mgz -odt float -ait ' T2transform];
    system([cmd ' >/dev/null']);
    auxT2=myMRIread([tempdir '/T2resampledLikeT1.mgz'],0,tempdir);
    [~,~,auxSl]=ind2sub(size(auxSEG.vol),find(auxSEG.vol==HippoLabelLeft | auxSEG.vol==HippoLabelRight));
    Sl=round(mean(auxSl));
    IT1=auxT1.vol(:,:,Sl); IT1=IT1/max(IT1(:))*255; IT1=uint8(IT1);
    IT2=auxT2.vol(:,:,Sl); IT2=IT2/max(IT2(:))*255; IT2=uint8(IT2);
    
    flipFile=[subjectDir '/' subjectName '/mri/transforms/T1_to_' suffixUser '.' suffix '.QC.gif'];
    if exist(flipFile,'file')>0, delete(flipFile); end;
    [imind,cm] = gray2ind(IT1,256);
    imwrite(imind,cm,flipFile,'gif', 'Loopcount',inf);
    [imind,cm] = gray2ind(IT2,256);
    imwrite(imind,cm,flipFile,'gif','WriteMode','append');
    
    disp('Done with registration!')
end



% Next: register image dump to automated segmentation
disp('Registering imageDump.mgz to hippocampal mask from ASEG')

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


% Target is masked aseg
targetRegFileName=[tempdir '/hippoAmygBinaryMask.mgz'];
targetRegFileNameCropped=[tempdir '/hippoAmygBinaryMask_autoCropped.mgz'];
ASEG=myMRIread([subjectDir '/' subjectName '/mri/aseg.mgz'],0,tempdir);
TARGETREG=ASEG;
if strcmp(side,'left')>0
    TARGETREG.vol=255*double(ASEG.vol==HippoLabelLeft | ASEG.vol==HippoLabelLeft+1);
else
    TARGETREG.vol=255*double(ASEG.vol==HippoLabelRight | ASEG.vol==HippoLabelRight+1);
end
myMRIwrite(TARGETREG,targetRegFileName,'float',tempdir);

highres=0; if mean(ASEG.volres)<0.99, highres=1; end
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



% Initial affine alignment based just on hippocampus

if 0==1  % This is to use a softened version (opening the mask first)
    aux=myMRIread(targetRegFileNameCropped,0,tempdir);
    strel=createSphericalStrel(1);
    aux.vol=imdilate(imerode(aux.vol>0,strel),strel);
    tmp1=bwdist(aux.vol);
    tmp2=bwdist(~aux.vol);
    D=zeros(size(aux.vol));
    D(~aux.vol)=-tmp1(~aux.vol);
    D(aux.vol)=tmp2(aux.vol)-1;
    rho=1;
    P=exp(D)./(exp(-D)+exp(D));
    P=round(255*P);
    aux.vol=P;
    myMRIwrite(aux,targetRegFileNameCropped,'float',tempdir);
end

if 1==1  % This is to use an opened version
    aux=myMRIread(targetRegFileNameCropped,0,tempdir);
    strel=createSphericalStrel(1);
    aux.vol=255*double(imdilate(imerode(aux.vol>0,strel),strel));
    myMRIwrite(aux,targetRegFileNameCropped,'float',tempdir);
end

% cmd=[FSpath 'kvlRegister imageDump.mgz ' targetRegFileNameCropped ' 3 2'];
% system(cmd);
% % system([cmd ' >/dev/null']);
% system('mv imageDump_coregistered.mgz imageDump.mgz' );


cmd=[FSpath '/mri_robust_register --mov imageDump.mgz  --dst ' targetRegFileNameCropped ...
    ' -lta trash.lta --mapmovhdr imageDump_coregistered.mgz  --sat 50'];
system(cmd);
% system([cmd ' >/dev/null']);
system('mv imageDump_coregistered.mgz imageDump.mgz' );

cmd=[FSpath '/mri_robust_register --mov imageDump.mgz  --dst ' targetRegFileNameCropped ...
    ' -lta trash.lta --mapmovhdr imageDump_coregistered.mgz --affine --sat 50'];
system(cmd);
% system([cmd ' >/dev/null']);
system('mv imageDump_coregistered.mgz imageDump.mgz' );



% Now, the idea is to refine the affine transform based on the hippo
% First, we prepare a modifided ASEG that we'll segment

% There's a bunch of labels in the ASEG don't have in our atlas...
ASEGbackup=ASEG;
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
ASEG.vol(ASEG.vol==0)=1;



% Write to disk
myMRIwrite(ASEG,'asegMod.mgz','float',tempdir);


% Eugenio July 2017
% We now merge hippo, amygdala and cortex in cheating image
ASEGcha=ASEG;
ASEGcha.vol(ASEGcha.vol==17)=3;
ASEGcha.vol(ASEGcha.vol==18)=3;
myMRIwrite(ASEGcha,'asegModCHA.mgz','float',tempdir);



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
% Eugenio July 2017
% [ synIm, transform ] = kvlReadCroppedImage( 'asegMod.mgz', boundingFileName );
[ synIm, transform ] = kvlReadCroppedImage( 'asegModCHA.mgz', boundingFileName );

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


% In this first step, we're cheating in that we're going to segment the "fake"
% synthetic image obtained from the ASEG. The means are equal to the ASEG labels

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

% Create mask for valid region of cuboid

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
        subplot( 3, 4, cheatingLabel )
        showImage( priors( :, :, :, cheatingLabel ) )
    end
    title('Priors for segmentation of fake intensity image')
end


cheatingImageBuffer=synImBuffer;
cheatingImageBuffer(~MASK)=0;
cheatingImage = kvlCreateImage( cheatingImageBuffer );

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
    
    % Eugenio November 2017: GEMS2 (note that it uses variances instead of precisions)
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
        disp(['Resolution ' num2str(multiResolutionLevel) ', iteration ' num2str(positionUpdatingIterationNumber)]);
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


% This is very important! We must warp the mesh to T2 space, which ideally would
% be perfectly overlapping, but in practice is often a bit off (motion
% between aquisitions in protocol?)
% First, make sure transform is ras2ras (rather than vox2vox...)
cmd=[FSpath '/mri_concatenate_lta -out_type   1 ' T2transform ' identity.nofile T1toT2.lta'];
system([cmd ' >/dev/null']);
T2transform=[tempdir '/T1toT2.lta'];


% cmd=[FSpath 'kvlApplyTransform imageDump.mgz '];
% fid=fopen(T2transform);
% while 1
%     tline = fgetl(fid);
%     if ~isempty(strfind(tline,'1 4 4')), break, end
% end
% for i=1:3
%     cmd=[cmd ' ' fgetl(fid)];
% end
% fclose(fid);
% system([cmd ' >/dev/null']);
% system('mv imageDump_transformed.mgz imageDump.mgz');

T=zeros(4,4);
T(4,4)=1;
fid=fopen(T2transform);
while 1
    tline = fgetl(fid);
    if ~isempty(strfind(tline,'1 4 4')), break, end
end
for i=1:3
    tline = fgetl(fid);
    [token,remain]=strtok(tline); T(i,1)=str2double(token);
    [token,remain]=strtok(remain); T(i,2)=str2double(token);
    [token,remain]=strtok(remain); T(i,3)=str2double(token);
    T(i,4)=str2double(remain);
end
fclose(fid);

aux=myMRIread('imageDump.mgz',0,tempdir);
aux.vox2ras0=T*aux.vox2ras0;
aux.vox2ras=T*aux.vox2ras;
aux.vox2ras1=T*aux.vox2ras1;
myMRIwrite(aux,'imageDump.mgz','float',tempdir);



% Next step is bias field correcting the T2 image.  To do so, we propagate
% the ASEG from T1 to T2 space, , merge left/right WM, cortex, etc, get rid
% of the labels that are at boundaries (to increase specificity) and use the
% rest  of the segmentation to bias field correct.
cmd=[FSpath '/mri_convert ' subjectDir '/' subjectName ...
    '/mri/aseg.mgz asegT2space.mgz -at ' T2transform ' -rt nearest -odt float'];
system([cmd ' >/dev/null']);
T2correctedFilename='T2corrected.mgz';

if ByPassBF>0
    
    system([FSpath '/mri_convert ' T2volumeFileName ' ' T2correctedFilename]);
    
else
    
    mri=myMRIread('asegT2space.mgz',0,tempdir);
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
    mriT2=myMRIread(T2volumeFileName,0,tempdir);
    biasFieldOrder=4;
    PSI=prepBiasFieldBase(size(mriT2.vol),biasFieldOrder);
    MASK=mri.vol>0 & ~isnan(mri.vol) & ~isnan(mriT2.vol);
    
    
    L=mri.vol(MASK);
    X=log(1+mriT2.vol(MASK));
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
    T2corr=log(1+mriT2.vol);
    for i=1:length(C)
        T2corr=T2corr-C(i)*PSI(:,:,:,i);
    end
    T2corr=exp(T2corr)-1; % back to natural
    T2corr(mriT2.vol==0)=0;
    T2corr(isnan(mriT2.vol))=0;
    aux=mri;
    aux.vol=T2corr;
    myMRIwrite(aux,T2correctedFilename,'float',tempdir);
    
end


% Upsample the T2 volume to work resolution
imageFileName='T2isotropic.mgz';
cmd=[FSpath '/mri_convert  ' T2correctedFilename ' ' imageFileName ' -rt  cubic -vs ' num2str(resolution) ' '  num2str(resolution)  ' '  num2str(resolution)];
system([cmd ' >/dev/null']);
% Eugenio: we want to mask out non-brain voxels and also the cerebellum,
% brainstem and 3rd/4th ventricles, which can be annoying later on.
system([FSpath '/mri_binarize --i asegMod.mgz --min 1.5 --dilate 1 --o asegModBinDilated.mgz >/dev/null']);
system([FSpath '/mri_convert asegModBinDilated.mgz asegModBinDilatedResampled.mgz -odt float -rt nearest -rl T2isotropic.mgz -at ' T2transform ' >/dev/null']);
system([FSpath '/mri_mask -T 0.5 T2isotropic.mgz  asegModBinDilatedResampled.mgz T2isotropic.mgz >/dev/null']);

% Eugenio: let's try masking anything that is not close to the hippo,
% brainstem and 3rd/4th ventricles, which can be annoying later on.
system([FSpath '/mri_binarize --i asegMod.mgz --min 16.5 --max 18.5 --o hippoMask.mgz >/dev/null']);
system([FSpath '/mri_convert hippoMask.mgz hippoMaskResampled.mgz -rt nearest -odt float -rl T2isotropic.mgz -at ' T2transform ' >/dev/null']);
system([FSpath '/mri_binarize --i  hippoMaskResampled.mgz --min 0.5 --dilate ' num2str(round(3/resolution)) '  --o hippoMaskResampledDilated.mgz  >/dev/null']);
system([FSpath '/mri_mask -T 0.5 T2isotropic.mgz  hippoMaskResampledDilated.mgz T2isotropic.mgz >/dev/null']);


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


% Eugenio July 2017

FreeSurferLabelGroups=[];

FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex','Left-Hippocampus','Left-Amygdala','subiculum-head','subiculum-body','Hippocampal_tail','GC-ML-DG-head','GC-ML-DG-body','CA4-head','CA4-body','presubiculum-head','presubiculum-body',...
    'CA1-head','CA1-body','parasubiculum','CA3-head','CA3-body','HATA','Lateral-nucleus','Paralaminar-nucleus',...
    'Basal-nucleus','Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
    'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus',...
    'Anterior-amygdaloid-area-AAA'};
FreeSurferLabelGroups{end+1}={'molecular_layer_HP-body','molecular_layer_HP-head'};
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
disp('Computing hyperparameters for estimation of Gaussians of T2 data')

if UseWholeBrainInHyperPar>0
    cmd=[FSpath '/mri_convert ' subjectDir '/' subjectName '/mri/aseg.mgz wmparcT2space.mgz -at ' T2transform ' -rt nearest -odt float'];
else
    cmd=[FSpath '/mri_convert ' subjectDir '/' subjectName '/mri/wmparc.mgz wmparcT2space.mgz -at ' T2transform ' -rt nearest -odt float'];
end
system([cmd '  >/dev/null']);

T2DATA=myMRIread(T2correctedFilename,0,tempdir);
WMPARC=myMRIread('wmparcT2space.mgz',0,tempdir);
T2DATA.vol(~imdilate(WMPARC.vol>0,createSphericalStrel(1)))=0;

nHyper=zeros(length(sameGaussianParameters),1);
meanHyper=zeros(length(sameGaussianParameters),1);
BGindex=0;
for g=1:length(sameGaussianParameters)
    labels=sameGaussianParameters{g};
    if any(labels==0)
        BGindex=g;
    end
    
    if UseWholeBrainInHyperPar>0
        
        if any(labels==3 | labels==17 | labels==18 | labels > 7000 | labels==226)  % gray matter + alveus/ML
            listMask=[3 42 17 53 18 54];
        elseif any(labels==2)  % white matter
            listMask=[2 41];
        elseif any(labels==26)    % accumbens area
            listMask=[26 58];
        elseif any(labels==4)  % CSF
            listMask=[4 43 14 15] ;
        elseif any(labels==0) % Background
            listMask=0;
        elseif any(labels==13)
            listMask=[13 52];
        elseif any(labels==12)
            listMask=[12 51];
        elseif any(labels==11)
            listMask=[11 50];
        elseif any(labels==10)
            listMask=[10 49];
        elseif any(labels==31)
            listMask=[31 63];
        elseif any(labels==28)
            listMask=[28 60];
        else
            listMask=[];
        end
        
    else
        
        if any(labels==3 | labels==17 | labels==18 | labels > 7000 | labels==226)  % gray matter + alveus/ML
            if strcmp(side,'left')>0, listMask=17; else, listMask=53; end
            
        elseif any(labels==2)  % white matter
            if strcmp(side,'left')>0, listMask=[3006 3007 3016]; else, listMask=[4006 4007 4016]; end
            
        elseif any(labels==26)    % accumbens area
            if strcmp(side,'left')>0, listMask=26; else, listMask=58; end
            
        elseif any(labels==4)  % CSF
            if strcmp(side,'left')>0, listMask=4; else, listMask=43; end
            
        elseif any(labels==0) % Background
            listMask=0;
        elseif any(labels==13)
            if strcmp(side,'left')>0, listMask=13; else, listMask=52; end
            
        elseif any(labels==12)
            if strcmp(side,'left')>0, listMask=12; else, listMask=51; end
            
        elseif any(labels==11)
            if strcmp(side,'left')>0, listMask=11; else, listMask=50; end
            
        elseif any(labels==10)
            if strcmp(side,'left')>0, listMask=10; else, listMask=49; end
            
        elseif any(labels==31)
            if strcmp(side,'left')>0, listMask=31; else, listMask=63; end
            
        elseif any(labels==28)
            if strcmp(side,'left')>0, listMask=28; else, listMask=60; end
            
        else
            listMask=[];
        end
    end
    if length(listMask)>0
        MASK=zeros(size(T2DATA.vol));
        for l=1:length(listMask)
            MASK=MASK | WMPARC.vol==listMask(l);
        end
        MASK=imerode(MASK,createSphericalStrel(1,T2DATA.volres));
        data=T2DATA.vol(MASK & T2DATA.vol>0);
        meanHyper(g)=median(data);
        % nHyper(g)=length(data)/2; % because voxels are 5 times smaller, this is like a /10
        % nHyper(g)=prod(T2DATA.volres)*length(data)/resolution^3/10; % /10 is out of thin air
        nHyper(g)=10+prod(T2DATA.volres)*length(data)/resolution^3;
    end
end
% if any nan, replace by background or mean of other labels
ind=find(isnan(meanHyper));
if ~isnan(meanHyper(BGindex)) && meanHyper(BGindex)~=0
    meanHyper(ind)=meanHyper(BGindex);
else
    aux=meanHyper(:);
    aux=aux(aux~=0);
    aux=aux(~isnan(aux));
    meanHyper(ind)=mean(aux);
end
nHyper(ind)=10;




% Here's the part where we simulate partial voluming!
disp('Estimating typical intensities of molecular layer and alveus')
WMind=[];
GMind=[];
MLind=[];
ALind=[];
CSFind=[];
FISSind=[];
for g=1:length(sameGaussianParameters)
    labels=sameGaussianParameters{g};
    if any(labels==2)
        WMind=g;
    end
    if any(labels==3)
        GMind=g;
    end
    if any(labels==245)  % Eugenio July 2017 (changed 214 by 245)
        MLind=g;
    end
    if any(labels==201)
        ALind=g;
    end
    if any(labels==4)
        CSFind=g;
    end
    if any(labels==215)
        FISSind=g;
    end
end


% Eugenio July 2017: again, rasterize priors one at the time

% priors = kvlRasterizeAtlasMesh( mesh, imageSize );
% priors=double(priors)/65535;
% suma=sum(priors,4);
% maskPriors=suma>.97;
% priors=priors./(eps+repmat(suma,[1 1 1 numberOfClasses]));
% priors=cumsum(priors,4);
% priorsCS=priors;
% clear priors;
% aux=rand(size(maskPriors));
% L=zeros(size(aux));
% for l=numberOfClasses:-1:1
%     L(aux<=priorsCS(:,:,:,l))=l;
% end
ss=[];
ss.Type='twister';
ss.Seed=58059898;
ss.State=uint32([3410437968 493000191 3937088844 2855859649 2055488253 4097603571 1434441578 3002745787 2261240414 2987560210 638430850 3764563185 4024560358 2845954880 1312077213 99046256 724219742 1509071479 2238342804 1354188848 562641151 3308086931 1111208124 3097183433 564927124 1520449620 929509102 859043821 3494335064 30251652 3274900264 633992030 12663028 2368225825 2382998372 4270645069 3706517678 1240495411 4109670805 1688366278 1985578320 3694255828 1103116839 3623746828 2370170267 3356379455 3812986031 3962487973 1117093995 1832124552 730164038 2688373380 1304038691 1133279105 1552678960 3485771784 2067203224 2672129996 2962006980 3811688264 2717503765 4253343185 531441539 1854060530 1488760948 1950262546 3699662348 2197811004 257371645 2932604825 432317016 1665242145 1342657666 510145158 4159693820 3167963749 2582698127 646844076 4224883034 3794533477 1714572503 2220104276 1818580560 2569089347 2486634147 149318545 1710902182 2913078558 1833836366 1230877566 2203621072 510333940 4218667072 2003694149 3899147439 3809686201 2691577592 2278093501 2776020136 3049065324 902390906 4248357917 2091592233 2920425164 2147920508 3195775060 1291700508 379824627 3244040512 281860832 1476677340 1862926681 2143256170 1511616358 3986711235 1472605582 3153067138 1406970208 1041278349 2464407436 2019913624 113009395 3926758644 3581575268 3854073575 840774617 2288814655 998711104 704513897 3965592645 581434087 644907376 264697324 1793526872 3075839507 2835430148 1107661823 1218930003 712025076 3593284091 1484105012 3932223704 2855162513 3994338726 784103664 814417492 1391485041 1602528562 2553133262 2346940117 732289055 1460085355 682492132 965862542 592114344 531690530 2964784340 1540752198 509261826 1465593050 1368541007 1119754638 3933145889 3409194796 1151629249 2543096396 1162952254 120686071 1072796031 1635344410 806990539 1949240388 2401163536 193410941 572929815 1826115942 1709368906 3823569594 3697467021 2606037820 3460208313 3093265774 1210164444 2422290984 1298978291 1597052602 2965991289 3328760101 169935670 3323761361 1034344567 3516510303 1016628606 4066423238 4067124643 3930184562 4104336116 1608506064 2207906500 812669751 1229667590 662572360 433350103 1284930189 953736903 3740888124 1279257809 2785226672 2413884704 1146523319 306943635 522178216 1644263167 4178778721 2407862653 3744063202 714083124 2677063966 4084664585 1903195336 1111708054 1750521315 1651027149 2095369713 96224378 1230846433 2255416769 283958389 848108856 1921761083 2263229615 3918442926 2278922043 2153996993 3053108311 3239184297 860769316 3265577964 1430072418 3870511250 205440931 106415127 2234205764 3974853631 1538962068 2749692830 1026086325 3882157666 2248910349 1703623791 4083086019 69088996 2089223477 3734553120 2288528134 2432850993 3775495635 2584140198 2597747462 886813745 3456528552 4126479183 2665494661 2056877935 2985474899 3334840693 1281194334 2805210650 2541614153 1608904543 1418149132 434633073 3450315291 2544246669 2337360302 3129734298 2591899202 956047506 571391348 4229678133 3339475591 2949518232 1541551447 2820896354 2634046742 3181523285 2013225811 1377073713 1421226137 428907431 3800857064 2742201395 2748902434 1645914435 611390336 113587038 1354467128 1554559976 3466792955 1787621901 1598573513 3157490081 1626593369 826436489 3671892643 3296344051 2643159537 2353297949 2087399669 108248550 3793060546 1667342514 3373407241 2948473492 2295503311 3095880945 648763624 634116976 4006007422 3125045892 2732167107 850138491 2594167551 1505554032 1687909106 2207354331 3265559100 2233507139 730839813 2829260419 3347277264 638525895 1736140020 3934299063 2335726919 1325176590 1863169739 4168059729 2705917166 2321630844 2657034283 961648448 381928833 2265082133 2083664274 2995096459 2848411249 2316887691 2591369222 1788138968 2311732152 2700483846 2576503411 3163343589 682564134 2055292476 824192588 3910582528 892349198 2130497059 1750018378 4175911463 1649173829 3717354654 1351865356 121532870 3453639903 1111023320 70709822 3953205824 1437739351 98337808 1615459686 159735467 2842564047 3273421125 1219401764 3904863290 876248348 3895207109 399482665 2638589731 4224607413 2288927228 1979741565 3083924839 3957274918 1296279006 3108756163 3081916387 89252673 3512991406 2866816427 2795366356 2516904582 3189375678 1380240981 2933642160 1251619690 42023279 263914470 436372455 1597753674 3365407821 1162660455 3135443575 2878706972 840746981 712563146 3586177285 4070568555 3488466883 1530585718 2961550898 769570043 2125673190 3205355305 3401304424 1724485905 693603705 223043625 3173316491 1308457270 2711024482 2804374120 3421132247 2728580540 3263618421 1102717826 3173417011 2724946155 1871187383 2079031531 769248695 4253236639 3311615520 3892062642 770082126 1357814258 1869469601 197863657 2900487501 3255380578 1458605047 3588007173 2501664808 150350842 3021307640 3976191881 2244533893 3946113114 2737382653 2623823297 147068974 347282931 965676538 2275342135 3512707344 1842292017 953709884 1882522525 3606021990 2306000931 826732670 1125648729 3616888784 3556707343 2345214961 1921304782 959778465 795295762 2991903222 1449749849 28980688 2946684790 4282850618 1450192187 1475726809 3938786168 1795023466 410052089 1589871039 3326940022 3818994470 2013723215 3044376647 958995322 3248307235 956083958 2072225296 3448052119 2447195658 2588071956 526006921 1983505568 1674635144 3264951886 2643183237 1318057102 3101414955 2151502608 3755605801 382470011 2878244912 1611623580 706493943 2025110342 1058251689 10454356 1861420073 3013200906 3729580155 738956171 3570317972 4203775383 1334062118 4257296000 3231017523 2818007598 2704948055 2817133652 2642265792 4104267226 1126945493 2596352021 1308859143 252889897 2689131845 344501177 1675665621 3259490044 4094531007 2795650781 1055448013 1035228504 1893552699 2769041259 3983042469 3425284573 3771000695 2741194467 3648153682 211324206 1165361053 1836003831 3002789777 1495668652 2831482182 2519660847 3363669882 2504184968 4141605892 2660518191 3741135717 1154098513 4194132382 2963865630 2547715443 2622625150 538064975 1338423991 2519859283 3951446572 149485389 3462323962 1975074798 4238312472 1812337747 4293462789 2706554319 150102680 2693831106 4288739420 3294196988 3330354565 4256435469 400255261 119744491 4224275315 4021312149 3322490141 2623776852 3053382684 3712861343 815139994 2482144449 3173066022 3288937807 560434293 3486803424 1917193992 1050249372 1487181270 2306492417 3814024376 903907259 3699798550 2757533368 2990388180 1902458037 2611644959 610140678 3632137931 269691837 548141251 2921620252 3314105880 3511432536 4205862093 1004944964 4042395795 3997534221 3892804844 2293573507 1720224801 1691692800 3886465211 1427362579 1214103038 1996640154 2630509595 833738701 1269966945 2011651581 2]);
rng(ss);
aux=rand(imageSize);
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
        M=prior>=(aux*65535);
        L=double(M);
        SET=M;
    else
        prior = kvlRasterizeAtlasMesh( mesh, imageSize, l-1 );
        sumpriors=sumpriors+prior;
        M=sumpriors>=(aux*65535) & SET==0;
        L(M)=l;
        SET(M)=1;
    end
end
suma=single(sumpriors)/65535;
maskPriors=suma>.97;
L(SET==0 & maskPriors)=size(reducedAlphas,2); % in case..


I=zeros(size(L));
for l=1:numberOfClasses
    if l==ALind || l==MLind
        I(L==l)=meanHyper(WMind);
    elseif l==FISSind
        I(L==l)=meanHyper(CSFind);
    else
        I(L==l)=meanHyper(l);
    end
end
I(~maskPriors | L==BGindex)=0;

semiW=(T2DATA.volres/resolution-1)/2; % in voxels
semiW(semiW<0)=0;
I_PV=GaussFilt3d(I,1.25*semiW);

if ~isempty(MLind)
    data=I_PV(L==MLind);
    meanHyper(MLind)=median(data);
    % nHyper(MLind)=300/resolution^3; % yes, I pull it out of thin air. But that's what hyperparameters are for :-D
    % nHyper(MLALind)=(nHyper(WMind)+nHyper(GMind))/2/100*2;
    nHyper(MLind)=(nHyper(WMind)+nHyper(GMind))/2;
end

if ~isempty(ALind)
    data=I_PV(L==ALind);
    meanHyper(ALind)=median(data);
    % nHyper(MLind)=300/resolution^3; % yes, I pull it out of thin air. But that's what hyperparameters are for :-D
    % nHyper(MLALind)=(nHyper(WMind)+nHyper(GMind))/2/100*2;
    nHyper(ALind)=(nHyper(WMind)+nHyper(GMind))/2;
end

if ~isempty(FISSind)
    data=I_PV(L==FISSind);
    meanHyper(FISSind)=median(data);
    nHyper(FISSind)=(nHyper(CSFind)+nHyper(GMind))/2;
end

%
% Multi-resolution scheme
%
% Specify here the size of the standard deviation of the Gaussian kernel used to smooth the priors/mesh. Use
% if you don't want to use multi-resolution
meshSmoothingSigmas = [ 1.5 .75 0 ]';
imageSmoothingSigmas = [0 0 0]';
maxItNos=[7 5 3];  % each iteration has 20 deformation steps

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
    
    
    
    
    % Iterations
    % Eugenio November 2017: we now do 30 instead of 20, since it's a bit faster
    maximumNumberOfIterations = maxItNos(multiResolutionLevel);  % Maximum number of iterations (includes one imaging model parameter estimation and
    positionUpdatingMaximumNumberOfIterations = 30;
    
    
    
    if FAST>0, maximumNumberOfIterations = 2; end  % in case we just wanna cruise throught it :-)
    
    
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
        %         % Ignore everything that's has zero intensity
        %         priors = priors( maskIndices, : );
        
        data = data( maskIndices );
        
        % Start EM iterations. Initialize the parameters if this is the
        % first time ever you run this
        posteriors = double( priors ) / 65535;
        if ( ( multiResolutionLevel == 1) & ( iterationNumber == 1 ) )
            
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
                    variances( classNumber ) = median(variances);
                    
                end
                
            end
            variances(variances==0)=median(variances); % added by Eugenio, prevents nans...
            
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
            
            minLogLikelihood = minLogLikelihood - sum( log( normalizer ) ) % This is what we're optimizing with EM
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
                    variances( classNumber ) = median(variances);
                    
                end
                
            end
            variances(variances==0)=median(variances); % added by Eugenio, prevents nans...
            
            
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

disp(['Fitting mesh to image data mask took ' num2str(etime(clock,time_ref_optimization)) ' seconds']);

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
system([FSpath '/mri_convert asegMod.mgz asmr1.mgz -rl T2isotropic.mgz -rt nearest -at '  T2transform ' -odt float >/dev/null']);
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

% reorient image.mgz
aux=zeros(size(tmp2.vol)+shiftNeg);
aux2=myMRIread('image.mgz',0,tempdir);
aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=aux2.vol;
aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
tmp3.vol=aux;
myMRIwrite(tmp3,'image.mgz','float',tempdir);

% Compute posteriors and volumes

% July 2017: make this memory efficient, and add new substructures. We also
% compute the segmentation along the way. Also, we do hippo and amygdala
% together. Ah! And we also do volumes for head and body
fidHP=fopen([tempdir '/volumesHippo.txt'],'w');
strOfInterestHP={'subiculum-body','subiculum-head','Hippocampal_tail','molecular_layer_HP-body','molecular_layer_HP-head','hippocampal-fissure','GC-ML-DG-body','GC-ML-DG-head','CA4-body','CA4-head','presubiculum-body','presubiculum-head','CA1-body','CA1-head','parasubiculum','fimbria','CA3-body','CA3-head','HATA'};
totVolHP=0;
foundHP=zeros([1,numberOfClasses]);
HPbodyList={'subiculum-body','CA1-body','presubiculum-body','molecular_layer_HP-body','CA3-body','GC-ML-DG-body','CA4-body','fimbria'};
HPheadList={'subiculum-head','presubiculum-head','CA1-head','parasubiculum','molecular_layer_HP-head','GC-ML-DG-head','CA4-head','CA3-head','HATA'};
totVolHPbody=0;
totVolHPhead=0;

fidAM=fopen([tempdir '/volumesAmygdala.txt'],'w');
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
        vol=resolution^3*(sum(double(post(:))/65535));
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
            myMRIwrite(tmp3,['posterior_' side '_' strtrim(lower(names(i,:))) '_T1_' suffix '_' suffixUser '.mgz'],'float',tempdir);
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



% MAP estimates
% Eugenio July 2017
% [~,inds]=max(posteriorsFull,[],4);
inds=L;

kk1=FreeSurferLabels(inds);
kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
aux=zeros(size(tmp2.vol)+shiftNeg);
aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
tmp3.vol=aux;
myMRIwrite(tmp3,'discreteLabels_all.mgz','float',tempdir);
tmp3.vol(tmp3.vol<200)=0;

% Eugenio July 2017
% tmp3.vol(tmp3.vol>226 & tmp3.vol<7000)=0;
tmp3.vol(tmp3.vol>246 & tmp3.vol<7000)=0;

tmp3.vol(tmp3.vol==201)=0; % alveus
tmp3Mask=getLargestCC(tmp3.vol>0);
tmp3.vol(~tmp3Mask)=0;
myMRIwrite(tmp3,'discreteLabels.mgz','float',tempdir);


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
myMRIwrite(tmp4,'discreteLabelsWholeBodyHead.mgz','float',tempdir);

% Second: head+body of each subfield
tmp4=tmp3;
tmp4.vol(tmp3.vol==233 | tmp3.vol==234)=204; % presubiculum
tmp4.vol(tmp3.vol==235 | tmp3.vol==236)=205; % subiculum
tmp4.vol(tmp3.vol==237 | tmp3.vol==238)=206; % CA1
tmp4.vol(tmp3.vol==239 | tmp3.vol==240)=208; % CA3
tmp4.vol(tmp3.vol==241 | tmp3.vol==242)=209; % CA4
tmp4.vol(tmp3.vol==243 | tmp3.vol==244)=210; % GC-DG
tmp4.vol(tmp3.vol==245 | tmp3.vol==246)=214; % ML

myMRIwrite(tmp4,'discreteLabelsMergedBodyHead.mgz','float',tempdir);

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
        seg=label*ones(size(mini));
    else
        dist=dmap(mask);
        m=dist<mini;
        mini(m)=dist(m);
        seg(m)=label;
    end
end
VOL(mask)=seg;
tmp4.vol(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6))=VOL;
myMRIwrite(tmp4,'discreteLabelsMergedBodyHeadNoMLorGCDG.mgz','float',tempdir);




% Eugenio July 2017
% I disabled this for now ...
if  MRFconstant>0
    
    % Eugenio July 2017
    error('MRF smoothing disabled for now');
    
    %     disp('Markov random field label "smoothing"...')
    %     unaryTermWeight=MRFconstant;
    %     [~,inds]=max(posteriorsFull,[],4);
    %     tmp=FreeSurferLabels(inds);
    %     kk=zeros(size(tmp)); kk(maskIndices)=1; tmp=tmp.*kk;
    %     tmp(tmp<200)=0; tmp(tmp>226 & tmp<7000)=0;
    %     [~,cropping]=cropLabelVol(tmp);
    %     Ct=zeros([cropping(4)-cropping(1)+1,cropping(5)-cropping(2)+1,cropping(6)-cropping(3)+1,numberOfClasses]);
    %     for c=1:numberOfClasses
    %         Ct(:,:,:,c)=-log(1e-12+(double(posteriorsFull(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6),c))/65535)*(1-1e-12))*unaryTermWeight;
    %     end
    %     varParas = [size(Ct,1); size(Ct,2); size(Ct,3); size(Ct,4); 100; 1e-3; 0.25; 0.11];
    %     penalty = ones([size(Ct,1),size(Ct,2),size(Ct,3)]);
    %     u = CMF3D_ML_mex(single(penalty), single(Ct), single(varParas));
    %     [~,ind] = max(u, [], 4);
    %     SEG=FreeSurferLabels(ind);
    %     SEG(SEG>226 & SEG<7000)=0; SEG(SEG<200)=0;
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

% Convert segmentation to 1 mm FreeSurfer Space (with and withoutresampling)
% Also, convert T2 scans (withoutresampling)
% Convert to 1 mm FreeSurfer Space
% Eugenio November 2017: smarter/smoother resampling
if SMOOTH_LABEL_RESAMPLE>0
    refFile=[subjectDir '/' subjectName '/mri/norm.mgz'];
    applyLTAsmoothLabels('discreteLabels.mgz',T2transform,'discreteLabelsResampledT1.mgz',refFile,1,FSpath,tempdir);
    applyLTAsmoothLabels('discreteLabelsWholeBodyHead.mgz',T2transform,'discreteLabelsWholeBodyHeadResampledT1.mgz',refFile,1,FSpath,tempdir);
    applyLTAsmoothLabels('discreteLabelsMergedBodyHead.mgz',T2transform,'discreteLabelsMergedBodyHeadResampledT1.mgz',refFile,1,FSpath,tempdir);
    applyLTAsmoothLabels('discreteLabelsMergedBodyHeadNoMLorGCDG.mgz',T2transform,'discreteLabelsMergedBodyHeadNoMLorGCDGResampledT1.mgz',refFile,1,FSpath,tempdir);
else
    system([FSpath '/mri_convert  discreteLabels.mgz  discreteLabelsResampledT1.mgz -rt nearest -odt float ' ...
        ' -rl ' subjectDir '/' subjectName '/mri/norm.mgz -ait ' T2transform]);
    system([FSpath '/mri_convert  discreteLabelsWholeBodyHead.mgz  discreteLabelsWholeBodyHeadResampledT1.mgz -rt nearest -odt float ' ...
        ' -rl ' subjectDir '/' subjectName '/mri/norm.mgz -ait ' T2transform]);
    system([FSpath '/mri_convert  discreteLabelsMergedBodyHead.mgz  discreteLabelsMergedBodyHeadResampledT1.mgz -rt nearest -odt float ' ...
        ' -rl ' subjectDir '/' subjectName '/mri/norm.mgz -ait ' T2transform]);
    system([FSpath '/mri_convert  discreteLabelsMergedBodyHeadNoMLorGCDG.mgz  discreteLabelsMergedBodyHeadNoMLorGCDGResampledT1.mgz -rt nearest -odt float ' ...
        ' -rl ' subjectDir '/' subjectName '/mri/norm.mgz -ait ' T2transform]);
end

system([FSpath '/mri_vol2vol --mov discreteLabels.mgz  --o discreteLabelsT1space.mgz --no-resample ' ...
    ' --targ ' subjectDir '/' subjectName '/mri/norm.mgz --lta-inv ' T2transform]);
system([FSpath '/mri_vol2vol --mov discreteLabelsWholeBodyHead.mgz  --o discreteLabelsWholeBodyHeadT1space.mgz --no-resample ' ...
    ' --targ ' subjectDir '/' subjectName '/mri/norm.mgz --lta-inv ' T2transform]);
system([FSpath '/mri_vol2vol --mov discreteLabelsMergedBodyHead.mgz  --o discreteLabelsMergedBodyHeadT1space.mgz --no-resample ' ...
    ' --targ ' subjectDir '/' subjectName '/mri/norm.mgz --lta-inv ' T2transform]);
system([FSpath '/mri_vol2vol --mov discreteLabelsMergedBodyHeadNoMLorGCDG.mgz  --o discreteLabelsMergedBodyHeadNoMLorGCDGT1space.mgz --no-resample ' ...
    ' --targ ' subjectDir '/' subjectName '/mri/norm.mgz --lta-inv ' T2transform]);

system([FSpath '/mri_vol2vol --mov ' T2volumeFileName '  --o T2inT1space.mgz --no-resample ' ...
    ' --targ ' subjectDir '/' subjectName '/mri/norm.mgz --lta-inv ' T2transform]);



if WRITE_POSTERIORS>0
    d=dir('posterior_*.mgz');
    for i=1:length(d)
        system([FSpath '/mri_vol2vol --mov ' d(i).name '  --o ' d(i).name(1:end-3) 'T1space.mgz --no-resample ' ...
            ' --targ ' subjectDir '/' subjectName '/mri/norm.mgz --lta-inv ' T2transform]);
    end
end

% Move to MRI directory
% Eugenio July 2017: add new segmentation maps, and simplified code
% (left/right)
system(['mv T2inT1space.mgz ' subjectDir '/' subjectName '/mri/' suffixUser '.FSspace.mgz']);

system(['mv discreteLabelsT1space.mgz ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygLabels-' suffixUser '.' suffix '.mgz']);
system(['mv discreteLabelsResampledT1.mgz ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygLabels-' suffixUser '.' suffix '.FSvoxelSpace.mgz']);
system(['mv discreteLabelsWholeBodyHeadT1space.mgz ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygLabels-' suffixUser '.' suffix '.HBT.mgz']);
system(['mv discreteLabelsWholeBodyHeadResampledT1.mgz ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygLabels-' suffixUser '.' suffix '.HBT.FSvoxelSpace.mgz']);
system(['mv discreteLabelsMergedBodyHeadT1space.mgz ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygLabels-' suffixUser '.' suffix '.FS60.mgz']);
system(['mv discreteLabelsMergedBodyHeadResampledT1.mgz ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygLabels-' suffixUser '.' suffix '.FS60.FSvoxelSpace.mgz']);
system(['mv discreteLabelsMergedBodyHeadNoMLorGCDGT1space.mgz ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygLabels-' suffixUser '.' suffix '.CA.mgz']);
system(['mv discreteLabelsMergedBodyHeadNoMLorGCDGResampledT1.mgz ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygLabels-' suffixUser '.' suffix '.CA.FSvoxelSpace.mgz']);

system(['mv volumesHippo.txt ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoSfVolumes-' suffixUser '.' suffix '.txt']);
system(['mv volumesAmygdala.txt ' subjectDir '/' subjectName '/mri/' side(1) 'h.amygNucVolumes-' suffixUser '.' suffix '.txt']);


if WRITE_POSTERIORS>0
    d=dir('posterior_*T1space.mgz');
    for i=1:length(d)
        if isempty(strfind(lower(d(i).name),'left-amygdala'))
            system(['mv ' d(i).name ' ' subjectDir '/' subjectName '/mri/' d(i).name(1:end-12) '.mgz']);
        end
    end
end
% if ByPassBF==0 &&  exist([subjectDir '/' subjectName '/T2corrected.mgz'],'file')==0
%     system(['mv T2corrected.mgz ' subjectDir '/' subjectName]);
% end


if WRITE_MESHES>0
    system(['mv warpedMesh.txt.gz   ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygMesh-' suffixUser '.' suffix '.txt.gz']);
    system(['mv warpedMeshNoAffine.txt.gz   ' subjectDir '/' subjectName '/mri/' side(1) 'h.hippoAmygMeshAtlasSpace-' suffixUser '.' suffix '.txt.gz']);
    system(['mv image2.mgz   ' subjectDir '/' subjectName '/mri/' side(1) 'h.imageForMesh-' suffixUser '.' suffix '.mgz']);
    fid=fopen([subjectDir '/' subjectName '/mri/' side(1) 'h.affineTransformMesh-' suffixUser '.' suffix '.txt'],'w');
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

