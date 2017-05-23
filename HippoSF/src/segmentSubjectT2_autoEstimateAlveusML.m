
% Segments the subfields from a T2 scan; requires the FS recon from the original MPRAGE

%
% segmentSubjectT2(subjectName,subjectDir,T2volumeFileName,resolution,atlasMeshFileName,atlasDumpFileName,compressionLUTfileName,K,side,optimizerType,suffix,MRFconstant,ByPassBF,UseWholeBrainInHyperPar)
%
% - subjectName: FreeSurfer subject name
% - subjectDir: FreeSurfer subject directory
% - T2volumeFileName: T2 image
% - resolution: voxel size at which we want to work (in mm).
% - atlasMeshFileName: the atlas to segment the data
% - atlasDumpFileName: corresponding imageDump.mgz (name *must* be imageDump.mgz)
% - compressionLUTfileName: corresponding compressionLUT.txt
% - K: stiffness of the mesh in the segmentation (around 0.05)
% - side: 'left' or 'right'
% - optimizerType: must be 'LM' or 'ConjGrad'
% - suffix: for output directory, e.g. 'T1based_GGAWLnoSimil','v10',...
% - suffixUser: for the user to identify outputs from different T2s
% - FSpath: path to FreeSurfer executables
% - MRFconstant (optional): make it >0 for MRF cleanup (5 is reasonable, larger is smoother)
%           It does NOT affect volumes, which are computed from soft posteriors anyway


function segmentSubjectT2_autoEstimateAlveusML(subjectName,subjectDir,T2volumeFileName,resolution,atlasMeshFileName,atlasDumpFileName,compressionLUTfileName,K,side,optimizerType,suffix,suffixUser,FSpath,MRFconstant,ByPassBF,UseWholeBrainInHyperPar)


% clear
% subjectName='HCA6018857_V1_A_highres_recon';
% subjectDir='/tmp/';
% T2volumeFileName='/autofs/space/tiamat_003/users/subjects/HCP/HCA6018857_V1_A_highres_recon/mri/TSE.FSspace.mgz';
% resolution=0.4;
% atlasMeshFileName='/usr/local/freesurfer/dev/average/HippoSF/atlas/AtlasMesh.gz';
% atlasDumpFileName='/usr/local/freesurfer/dev/average/HippoSF/atlas/AtlasDump.mgz';
% compressionLUTfileName='/usr/local/freesurfer/dev/average/HippoSF/atlas/compressionLookupTable.txt';
% K=0.05;
% side='right';
% % optimizerType='LM';
% optimizerType='ConjGrad';
% suffix='T2';
% suffixUser='Eugenio';
% FSpath='/usr/local/freesurfer/dev/bin/';
% MRFconstant=0;
% ByPassBF=1;
% UseWholeBrainInHyperPar=0;



% clear
% subjectName='subject1';
% subjectDir='/autofs/space/panamint_005/users/iglesias/data/WinterburnHippocampalAtlas/FSdirConformed/';
% T2volumeFileName='/autofs/space/panamint_005/users/iglesias/data/WinterburnHippocampalAtlas/subject1/T2.0.6.nii.gz';
% resolution=0.4;
% atlasMeshFileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/output/CurrentMeshCollection30.gz';
% atlasDumpFileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/output/imageDumpWithAmygdala.mgz';
% compressionLUTfileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/compressionLookupTable.txt';
% K=0.05;
% side='right';
% % optimizerType='LM';
% optimizerType='ConjGrad';
% suffix='T2basedHRtest';
% suffixUser='T2';
% FSpath='/usr/local/freesurfer/dev/bin/';
% MRFconstant=0;
% ByPassBF=1;
% UseWholeBrainInHyperPar=0;


% clear
% % subjectName='quicknick';
% % subjectDir='/autofs/space/terrier_002/users/nicks/subjects/nicks_recons/';
% % T2volumeFileName=[subjectDir '/' subjectName '/mri/T2.mgz'];
% subjectName='036_S_2378_17-Jul-2013';
% subjectDir='/autofs/space/panamint_005/users/iglesias/data//ADNI_FS_subjDirFS53/';
% T2volumeFileName='/autofs/space/panamint_005/users/iglesias/data//ADNI_FS_subjDirFS53/036_S_2378_17-Jul-2013.T2.nii.gz'; % good case
% resolution=0.4;
% %
% % atlasMeshFileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_131115_CJ_GD_allBuckner/output/CurrentMeshCollection30_fixed_splittedClasses.gz';
% % atlasDumpFileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_131115_CJ_GD_allBuckner/output/imageDumpWithAmygdala.mgz';
% % compressionLUTfileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_131115_CJ_GD_allBuckner/compressionLookupTableSplitted.txt';
% atlasMeshFileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/output/CurrentMeshCollection30.gz';
% atlasDumpFileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/output/imageDumpWithAmygdala.mgz';
% compressionLUTfileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/compressionLookupTable.txt';
% K=0.05;
% side='left';
% % optimizerType='LM';
% optimizerType='ConjGrad';
% suffix='T2based_atlas_131115_autoMLAL_justTesting';
% suffixUser='T2test';
% FSpath='/usr/local/freesurfer/stable5_3_0/bin/';
% MRFconstant=0;
% ByPassBF=1;
% UseWholeBrainInHyperPar=0;



% subjectName='subject3';
% subjectDir='/autofs/space/panamint_005/users/iglesias/data/WinterburnHippocampalAtlas/FSdir/';
% T2volumeFileName=['/autofs/space/panamint_005/users/iglesias/data/WinterburnHippocampalAtlas/' subjectName '/T2.nii'];
% resolution=0.3;
% atlasMeshFileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/output/CurrentMeshCollection30.gz';
% atlasDumpFileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/output/imageDumpWithAmygdala.mgz';
% compressionLUTfileName='/autofs/space/panamint_005/users/iglesias/atlases/atlasHippoBuckner_140408_CJ_GD_allBuckner/compressionLookupTable.txt';
% K=0.05;
% side='left';
% % optimizerType='LM';
% optimizerType='ConjGrad';
% suffix='T2based_atlas_140408_T2based';
% FSpath='/autofs/space/panamint_005/users/iglesias/software/Atlas3D_131025-Release/bin/';
% MRFconstant=5;

DEBUG=0;
FAST=0; % set it to one to optimize just a bit (go through code fast)
WRITE_POSTERIORS=0;
aux=getenv('WRITE_POSTERIORS');
if ~isempty(aux)
    if str2double(aux)>0
        WRITE_POSTERIORS=1;
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
elseif strcmp(optimizerType,'LM')==0 && strcmp(optimizerType,'ConjGrad')==0
    error('Optimizer type must be ''LM'' or ''ConjGrad''');
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
    if isunix
        addpath('/cluster/koen/eugenio/GEMS-Release-linux/bin')
    elseif ismac
        addpath('/cluster/koen/eugenio/GEMS-Release-mac/bin')
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
[ synIm, transform ] = kvlReadCroppedImage( 'asegMod.mgz', boundingFileName );
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
% The variances we just set to a small value.
% We also build compressed classes (which group together the hippocampus and the amygdala)
FreeSurferLabelGroups=[];
FreeSurferLabelGroups{end+1}={'Left-Amygdala','Lateral-nucleus','Paralaminar-nucleus','Basal-nucleus','Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
    'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus','Anterior-amygdaloid-area-AAA'};
FreeSurferLabelGroups{end+1}={'Left-Hippocampus','alveus','subiculum','Hippocampal_tail','molecular_layer_HP','GC-ML-DG','CA4','CA1','CA3','HATA','fimbria','presubiculum','parasubiculum','hippocampal-fissure','Left-hippocampus-intensity-abnormality'};
FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex'};
FreeSurferLabelGroups{end+1}={'Left-Cerebral-White-Matter'};
FreeSurferLabelGroups{end+1}={'Left-Lateral-Ventricle'};
FreeSurferLabelGroups{end+1}={'Left-choroid-plexus'};
FreeSurferLabelGroups{end+1}={'Background','Background-CSF','Background-vessels','Background-tissue','Unknown'};
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


cheatingMeans=zeros(length( sameGaussianParameters),1);
cheatingVariances=0.01*ones(length( sameGaussianParameters),1);
for l=1:length(sameGaussianParameters)
    label= sameGaussianParameters{l}(1);
    if (label>=200 && label<=226) || label==106,  cheatingMeans(l)=17; % HIPPO SF or hippo abnormality -> HIPPO
    elseif label>=7000,  cheatingMeans(l)=18; % AMYGDALOID SUBNUCLEI -> AMYGDALA
    elseif label==0, cheatingMeans(l)=1; % BACKGROUND is 1 instead of 0
    elseif label==999, cheatingMeans(l)=55; cheatingVariances(l)=55^2; % This is the generic, "suspicious" label ...
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


priors = kvlRasterizeAtlasMesh( mesh, synSize );
if ~isdeployed && DEBUG>0
    figure
    for cheatingLabel = 1 : size( reducedAlphas, 2 )
        subplot( 3, 4, cheatingLabel )
        showImage( priors( :, :, :, cheatingLabel ) )
    end
    title('Priors for segmentation of fake intensity image')
end

MASK=imerode(sum(double(priors)/65535,4)>0.99,createSphericalStrel(3));
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
if strcmp(optimizerType,'LM')>0
    maxIt=[100,50];
else
    maxIt=[200,100];
end


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
    
    % Set up the black box optimizer for the mesh nodes
    if ( exist( 'optimizer', 'var' ) == 1 )
        % The optimizer is very memory hungry when run in multithreaded mode.
        % Let's clear any old ones we may have lying around
        kvlClear( optimizer );
    end
    
    % Now the optimization per-se
    % [ cost gradient ] = kvlEvaluateMeshPosition( mesh,cheatingImage, transform, cheatingMeans, cheatingVariances );
    if strcmp(optimizerType,'LM')>0
        cheatingOptimizer = kvlGetLevenbergMarquardtOptimizer( mesh, cheatingImage, transform );
        maximalDeformationStopCriterion = 0.001;
        relativeChangeInCostStopCriterion = 1e-7;
    else
        cheatingOptimizer = kvlGetConjugateGradientOptimizer( mesh, cheatingImage, transform );
        maximalDeformationStopCriterion = 1e-10;  % worth it to be precise, this is fast compared to the latter optimization anyway...
        relativeChangeInCostStopCriterion = 1e-10;
    end
    maxpuin=maxIt(multiResolutionLevel);
    kvlSetOptimizerProperties( cheatingOptimizer, cheatingMeans, 1./reshape(cheatingVariances,[1 1 length(cheatingVariances)]) );
    
    if FAST>0
        maxpuin=20;
    end
    
    for positionUpdatingIterationNumber = 1 : maxpuin
        disp(['Resolution ' num2str(multiResolutionLevel) ', iteration ' num2str(positionUpdatingIterationNumber)]);
        % Calculate a good step. The first one is very slow because of various set-up issues
        tic
        [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStep( cheatingOptimizer , 1);
        elapsedTime = toc;
        disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
        minLogLikelihoodTimesPrior
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
    MASK=mri.vol>0;
    
    
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
priors = kvlRasterizeAtlasMesh( mesh, imageSize ); % Without specifying a specific label, will rasterize all simultaneously
if ~isdeployed && DEBUG>0
    colorCodedPriors = kvlColorCodeProbabilityImages( priors, colors );
    figure
    showImage( colorCodedPriors );
    title('Color-coded priors')
    clear colorCodedPriors
end


% We're not interested in image areas that fall outside our cuboid ROI where our atlas is defined. Therefore,
% generate a mask of what's inside the ROI. Also, by convention we're skipping all voxels whose intensities
% is exactly zero (which we do to remove image areas far away for FreeSurfer's ASEG hippo segmentation) -
% also include that in the mask
mask = imerode(single( sum( priors, 4 ) / 65535 ) > 0.99,createSphericalStrel(5));
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
FreeSurferLabelGroups{end+1}={'Left-Cerebral-Cortex','Left-Hippocampus','Left-Amygdala','subiculum','Hippocampal_tail','GC-ML-DG','CA4','presubiculum',...
    'CA1','parasubiculum','CA3','HATA','Lateral-nucleus','Paralaminar-nucleus',...
    'Basal-nucleus','Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
    'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus','Anterior-amygdaloid-area-AAA'};
FreeSurferLabelGroups{end+1}={'Left-Cerebral-White-Matter','fimbria'};
FreeSurferLabelGroups{end+1}={'molecular_layer_HP'};
FreeSurferLabelGroups{end+1}={'alveus'};
FreeSurferLabelGroups{end+1}={'hippocampal-fissure'};
FreeSurferLabelGroups{end+1}={'Left-Lateral-Ventricle','Background-CSF','SUSPICIOUS','Left-hippocampus-intensity-abnormality'};
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
    if any(labels==214)
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

priors = kvlRasterizeAtlasMesh( mesh, imageSize );
priors=double(priors)/63535;
suma=sum(priors,4);
maskPriors=suma>.97;
priors=priors./(eps+repmat(suma,[1 1 1 numberOfClasses]));
priors=cumsum(priors,4);
priorsCS=priors;
clear priors;
aux=rand(size(maskPriors));
L=zeros(size(aux));
for l=numberOfClasses:-1:1
    L(aux<=priorsCS(:,:,:,l))=l;
end
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
    
    % Set up the black box optimizer for the mesh nodes
    if ( exist( 'optimizer', 'var' ) == 1 )
        % The optimizer is very memory hungry when run in multithreaded mode.
        % Let's clear any old ones we may have lying around
        try
            kvlClear( optimizer );
        catch ME
        end
    end
    
    
    % Now the optimization per-se
    if strcmp(optimizerType,'LM')>0
        optimizer = kvlGetLevenbergMarquardtOptimizer( mesh, image, transform );
        maximalDeformationStopCriterion = 0.01;
        positionUpdatingMaximumNumberOfIterations = 20;
    else
        optimizer = kvlGetConjugateGradientOptimizer( mesh, image, transform );
        maximalDeformationStopCriterion = 1e-10;
        positionUpdatingMaximumNumberOfIterations = 20;
    end
    
    
    maximumNumberOfIterations = maxItNos(multiResolutionLevel);  % Maximum number of iterations (includes one imaging model parameter estimation and
    % one deformation optimization; the latter always does 20 steps).
    
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
        priors = kvlRasterizeAtlasMesh( mesh, imageSize );
        priors = reshape( priors, [ prod( imageSize ) numberOfClasses ] ); % Easier to work with vector notation in the computations
        
        % Ignore everything that's has zero intensity
        priors = priors( maskIndices, : );
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
            posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfClasses ] );
            
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
        haveMoved = false; % Keep track if we've ever moved or not
        kvlSetOptimizerProperties( optimizer, means', reshape(1./variances,[1 1 length(variances)]));
        for positionUpdatingIterationNumber = 1 : positionUpdatingMaximumNumberOfIterations
            % Calculate a good step. The first one is very slow because of various set-up issues
            disp(['Resolution level ' num2str(multiResolutionLevel) ' iteration ' num2str(iterationNumber) ' deformation iterations ' num2str(positionUpdatingIterationNumber)]);
            tic
            [ minLogLikelihoodTimesPrior, maximalDeformation ] = kvlDeformOneStep( optimizer, 1 );
            elapsedTime = toc;
            disp( [ 'Did one deformation step of max. ' num2str( maximalDeformation )  ' voxels in ' num2str( elapsedTime ) ' seconds' ] )
            minLogLikelihoodTimesPrior
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
kvlClear( optimizer )

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

% Get the priors as dictated by the current mesh position
data = double( reshape( imageBuffer, [ prod( imageSize ) 1 ] ) ); % Easier to work with vector notation in the computations
priors = kvlRasterizeAtlasMesh( mesh, imageSize );
priors = reshape( priors, [ prod( imageSize ) numberOfClasses ] );  % Easier to work with vector notation in the computations

% Ignore everything that's has zero intensity
priors = priors( maskIndices, : );
data = data( maskIndices );

% Calculate the posteriors
posteriors = zeros( size( priors ), 'double' );
for classNumber = 1 : numberOfClasses
    % Get the parameters from the correct Gaussian
    mu = means( reducingLookupTable( classNumber ) );
    variance = variances( reducingLookupTable( classNumber ) );
    prior = single( priors( :, classNumber ) ) / 65535;
    
    posteriors( :, classNumber ) = ( exp( -( data - mu ).^2 / 2 / variance ) .* prior ) ...
        / sqrt( 2 * pi * variance);
end
normalizer = sum( posteriors, 2 ) + eps;
posteriors = posteriors ./ repmat( normalizer, [ 1 numberOfClasses ] );
posteriors = uint16( round( posteriors * 65535 ) );

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
priorsFull = kvlRasterizeAtlasMesh( mesh, imageSize );
posteriorsFull=priorsFull;

fid=fopen([tempdir '/volumesHippo.txt'],'w');
% strOfInterest={'alveus','subiculum','Hippocampal_tail','molecular_layer_HP','hippocampal-fissure','GC-ML-DG','CA4','presubiculum','CA1','parasubiculum','fimbria','CA3','HATA'};
% Ignore alveus... good for guiding, but volumes are so unreliable :-S
strOfInterest={'subiculum','Hippocampal_tail','molecular_layer_HP','hippocampal-fissure','GC-ML-DG','CA4','presubiculum','CA1','parasubiculum','fimbria','CA3','HATA'};
totVol=0;
found=zeros(1,size(priorsFull,4));
for i=1:size(priorsFull,4)
    tmp=posteriorsFull(:,:,:,i);
    tmp(maskIndices)=posteriors(:,i);
    posteriorsFull(:,:,:,i)=tmp;
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
        vol=resolution^3*(sum(sum(sum(double(posteriorsFull(:,:,:,i))/65535))));
        fprintf(fid,'%s %f\n',str,vol);
        if isempty(strfind(lower(names(i,:)),'hippocampal-fissure'))  % don't count the fissure towards the total volume
            totVol=totVol+vol;
            
            if WRITE_POSTERIORS>0
                kk1=double(posteriorsFull(:,:,:,i))/65535;
                kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
                aux=zeros(size(tmp2.vol)+shiftNeg);
                aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
                aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
                tmp3.vol=aux;
                myMRIwrite(tmp3,['posterior_' side '_' strtrim(lower(names(i,:))) '_' suffixUser '_' suffix '.mgz'],'float',tempdir);
            end
            
        end
    end
end
if sum(found>0)>1
    fprintf(fid,'Whole_hippocampus %f\n',totVol);
    fclose(fid);
else
    fclose(fid);
    delete([tempdir '/volumesHippo.txt']);
end


fid=fopen([tempdir '/volumesAmygdala.txt'],'w');
strOfInterest={'Left-Amygdala','Lateral-nucleus','Paralaminar-nucleus',...
    'Basal-nucleus','Hippocampal-amygdala-transition-HATA','Accessory-Basal-nucleus','Amygdala-background',...
    'Corticoamygdaloid-transitio','Central-nucleus','Cortical-nucleus','Medial-nucleus','Anterior-amygdaloid-area-AAA'};
totVol=0;
found=zeros(1,size(priorsFull,4));
for i=1:size(priorsFull,4)
    tmp=posteriorsFull(:,:,:,i);
    tmp(maskIndices)=posteriors(:,i);
    posteriorsFull(:,:,:,i)=tmp;
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
        vol=resolution^3*(sum(sum(sum(double(posteriorsFull(:,:,:,i))/65535))));
        fprintf(fid,'%s %f\n',str,vol);
        if isempty(strfind(lower(names(i,:)),'hippocampal-fissure'))  % don't count the fissure towards the total volume
            totVol=totVol+vol;
            
            if WRITE_POSTERIORS>0
                kk1=double(posteriorsFull(:,:,:,i))/65535;
                kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
                aux=zeros(size(tmp2.vol)+shiftNeg);
                aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
                aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
                tmp3.vol=aux;
                myMRIwrite(tmp3,['posterior_' side '_' strtrim(lower(names(i,:))) '_' suffixUser '_' suffix '.mgz'],'float',tempdir);
            end
            
        end
    end
end
if sum(found>0)>1
    fprintf(fid,'Whole_amygdala %f\n',totVol);
    fclose(fid);
else
    fclose(fid);
    delete([tempdir '/volumesAmygdala.txt']);
end

% MAP estimates
[~,inds]=max(posteriorsFull,[],4);
kk1=FreeSurferLabels(inds);
kk2=zeros(size(kk1)); kk2(maskIndices)=1; kk1=kk1.*kk2;
aux=zeros(size(tmp2.vol)+shiftNeg);
aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(kk1,[2 1 3]);
aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
tmp3.vol=aux;
myMRIwrite(tmp3,'discreteLabels_all.mgz','float',tempdir);
tmp3.vol(tmp3.vol<200)=0;
tmp3.vol(tmp3.vol>226 & tmp3.vol<7000)=0;
tmp3.vol(tmp3.vol==201)=0; % alveus
tmp3Mask=getLargestCC(tmp3.vol>0);
tmp3.vol(~tmp3Mask)=0;
myMRIwrite(tmp3,'discreteLabels.mgz','float',tempdir);



if  MRFconstant>0
    disp('Markov random field label "smoothing"...')
    unaryTermWeight=MRFconstant;
    [~,inds]=max(posteriorsFull,[],4);
    tmp=FreeSurferLabels(inds);
    kk=zeros(size(tmp)); kk(maskIndices)=1; tmp=tmp.*kk;
    tmp(tmp<200)=0; tmp(tmp>226 & tmp<7000)=0;
    [~,cropping]=cropLabelVol(tmp);
    Ct=zeros([cropping(4)-cropping(1)+1,cropping(5)-cropping(2)+1,cropping(6)-cropping(3)+1,numberOfClasses]);
    for c=1:numberOfClasses
        Ct(:,:,:,c)=-log(1e-12+(double(posteriorsFull(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6),c))/65535)*(1-1e-12))*unaryTermWeight;
    end
    varParas = [size(Ct,1); size(Ct,2); size(Ct,3); size(Ct,4); 100; 1e-3; 0.25; 0.11];
    penalty = ones([size(Ct,1),size(Ct,2),size(Ct,3)]);
    u = CMF3D_ML_mex(single(penalty), single(Ct), single(varParas));
    [~,ind] = max(u, [], 4);
    SEG=FreeSurferLabels(ind);
    SEG(SEG>226 & SEG<7000)=0; SEG(SEG<200)=0;
    
    data=zeros(size(inds));
    data(cropping(1):cropping(4),cropping(2):cropping(5),cropping(3):cropping(6))=SEG;
    aux=zeros(size(tmp2.vol)+shiftNeg);
    aux(1+shiftNeg(1):shiftNeg(1)+size(tmp2.vol,1),1+shiftNeg(2):shiftNeg(2)+size(tmp2.vol,2),1+shiftNeg(3):shiftNeg(3)+size(tmp2.vol,3))=permute(data,[2 1 3]);
    aux=aux(1+shiftPos(1):end,1+shiftPos(2):end,1+shiftPos(3):end);
    tmp3.vol=aux;
    tmp3Mask=getLargestCC(tmp3.vol>0);
    tmp3.vol(~tmp3Mask)=0;
    myMRIwrite(tmp3,'discreteLabels_MRF.mgz','float',tempdir);
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
system([FSpath '/mri_convert  discreteLabels.mgz  discreteLabelsResampledT1.mgz -rt nearest -odt float ' ...
    ' -rl ' subjectDir '/' subjectName '/mri/norm.mgz -ait ' T2transform]);

system([FSpath '/mri_vol2vol --mov discreteLabels.mgz  --o discreteLabelsT1space.mgz --no-resample ' ...
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
system(['mv T2inT1space.mgz ' subjectDir '/' subjectName '/mri/' suffixUser '.FSspace.mgz']);
if strcmp(side,'right')>0
    system(['mv discreteLabelsT1space.mgz ' subjectDir '/' subjectName '/mri/rh.hippoSfLabels-' suffixUser '.' suffix '.mgz']);
    system(['mv discreteLabelsResampledT1.mgz ' subjectDir '/' subjectName '/mri/rh.hippoSfLabels-' suffixUser '.' suffix '.FSvoxelSpace.mgz']);
    system(['mv volumesHippo.txt ' subjectDir '/' subjectName '/mri/rh.hippoSfVolumes-' suffixUser '.' suffix '.txt']);
    system(['mv volumesAmygdala.txt ' subjectDir '/' subjectName '/mri/rh.hippoSfVolumes-' suffixUser '.amygdala.' suffix '.txt']);
else
    system(['mv discreteLabelsT1space.mgz ' subjectDir '/' subjectName '/mri/lh.hippoSfLabels-' suffixUser '.' suffix '.mgz']);
    system(['mv discreteLabelsResampledT1.mgz ' subjectDir '/' subjectName '/mri/lh.hippoSfLabels-' suffixUser '.' suffix '.FSvoxelSpace.mgz']);
    system(['mv volumesHippo.txt ' subjectDir '/' subjectName '/mri/lh.hippoSfVolumes-' suffixUser '.' suffix '.txt']);
    system(['mv volumesAmygdala.txt ' subjectDir '/' subjectName '/mri/lh.hippoSfVolumes-' suffixUser '.amygdala.' suffix '.txt']);
end
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

