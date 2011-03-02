function [dr,dl]=cc_cut_adf(subject, name,th_pval)
%           For each subject "subject":
% Computes the Dice coefficients measuring the overlap
% of the WM volume in right and left hemispheres to check  
%    if the corpus_callosum is correctly located. 
%
% Uses .lta transform and p values
%


%
% cc_cut_afd.m
%
% Original Author: Laurence Wastiaux
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.4 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%




if (nargin<1 | nargin>3)
    msg=sprintf('USAGE: [Dr,Dl]=cc_cut_adf(Subject)');
    disp(msg)
end

%%  Load the .lta Talairach transform matrix  %%
%% -----------------------------------------  %%
matname=strcat(subject, '/mri/transforms/talairach.lta');
fid=fopen(matname, 'r');
if fid == -1
    message=sprintf('Cannot open the file %s',matname);
    disp(message);
else
    while feof(fid) == 0
        linef=fgetl(fid);
        if(findstr(linef, '1 4 4'))
            pos=ftell(fid);
            break
        end    
    end
    A=(fscanf(fid, '%g',12))';
    A=(reshape(A, [4 3]))';
    status=fclose(fid);
end
Tal_transform=A;

%%  Load the wm & aseg volume    %%
%% ----------------------------  %%

SubjectDir=strcat(subject,'/');
CorDir1=strcat(subject,'/mri/wm/');
CorDir2=strcat(subject,'/mri/aseg/');
d1=dir(CorDir1);
d2=dir(CorDir2);

if (length(d1)<3)
    mess1=sprintf('Cannot find the coronals of the volume wm for the subject %s',subject);
    error(mess1)
elseif (length(d2)<3)
    mess2=sprintf('Cannot find the coronals of the volume aseg for the subject %s',subject);
    error(mess2)
else
    
    [vol1 M1]=load_cor2(SubjectDir,'wm');
    [vol2 M2]=load_cor2(SubjectDir,'aseg');
    
end

%% Read position of the cutting plane  %%
%% Assumes that a log file is present  %%
%% ----------------------------------  %%

%logfilename=strcat(subject,'/mri/ponscc.log'); % Location where the log file can normally be found
logfilename=strcat('/autofs/homes/001/wastiaux/adf/cutting_planes/logfiles_neo2/',name,'/ponscc.log'); %temporary location for the subjects of the training set

fp=fopen(logfilename, 'r');
if (fp==-1)
    msg=sprintf('Cannot open file %s, try to generate a new log file', logfilename);
    disp(msg)
    subject_wm=strcat(subject,'/mri/wm');
    cmd=sprintf('~/work/dev/mri_fill/mri_fill -l ~/tmp.log %s ~/filled.mgz',subject_wm);
    s=unix(cmd);
    logfilename2=sprintf('~/tmp.log');
    fp=fopen(logfilename2, 'r');
    if (fp==-1)
        mess=sprintf('Cannot open file %s', logfilename2);
        error(mess)
    end
end
if (fp~=-1)
    while feof(fp) == 0
        s=fscanf(fp,'%s',1);
        if (strcmp(s, 'CC:')>0 | strcmp(s, 'CC-CRS:')>0);
            x_cc=fscanf(fp,'%d',1); % only x_cc is needed since the cutting plane is defined in a sagittal slice in the talairach space
            fscanf(fp,'%s',1);
            y_cc=fscanf(fp,'%d',1);
            fscanf(fp,'%s',1);
            z_cc=fscanf(fp,'%d',1);
            break
        end    
    end
end
V_cc=[x_cc ; y_cc ; z_cc ];

%%    Compute statistics   %%
%% ----------------------- %%
rh_wm_1=0; lh_wm_1=0;
rh_wm_2=0; lh_wm_2=0;
rh_wm_12=0; lh_wm_12=0;

if ( (size(vol1)) ~= (size(vol2)) )
    msg=sprintf('Wm and aseg volumes do not have the same size');
    error(msg);
else
    sz=size(vol1);
    for i=1:sz(1)
        for j=1:sz(2)
            for k=1:sz(3)
                if (vol2(i,j,k)==41 | vol2(i,j,k)==46) % the cerrebellum wm is taken into account 
                    rh_wm_2 = rh_wm_2 +1;
                end
                if (vol2(i,j,k)==2 | vol2(i,j,k)==7)
                    lh_wm_2 = lh_wm_2 +1;
                end
                if (vol1(i,j,k)~=0)
                    v=[j i k 1]';
                    v_tal=Tal_transform*v;
                    if (v_tal(1) < V_cc(1))
                        rh_wm_1 = rh_wm_1 +1;
                    else
                        lh_wm_1 = lh_wm_1 +1;
                    end
                end
                if (vol1(i,j,k)~=0 && (vol2(i,j,k)==41 | vol2(i,j,k)==46))
                    v=[j i k 1]';
                    v_tal=Tal_transform*v;
                    if (v_tal(1) < V_cc(1))
                        rh_wm_12 = rh_wm_12 +1;
                    end
                end
                if (vol1(i,j,k)~=0 && (vol2(i,j,k)==2 | vol2(i,j,k)==7))
                    v=[j i k 1]';
                    v_tal=Tal_transform*v;
                    if (v_tal(1) > V_cc(1))
                        lh_wm_12 = lh_wm_12 +1;
                    end
                end
                    
            end
        end
    end
end

%% Compute Dice coefficients  %%
%% ------------------------- %%

dr= 2*rh_wm_12/(rh_wm_1+rh_wm_2);
dl= 2*lh_wm_12/(lh_wm_1+lh_wm_2);
[rpval lpval]=compute_pval(dr, dl);

if (rpval < th_pval)
    mess1=sprintf('The Dice coefficient for the right hemisphere is too low (Dr=%3g, rpval=%3g)', dr, rpval);
    disp(mess1)
end
if (lpval < th_pval)
    mess2=sprintf('The Dice coefficient for the left hemisphere is too low (Dl=%3g, lpval=%.3g)', dl, lpval);
    disp(mess2)
end
if( (rpval>=th_pval) & (lpval >= th_pval))
    mess=sprintf('Vertical cutting plane location OK (Dl=%3g, lpval=%.3g,Dr=%3g,rpval=%3g)',dl, lpval, dr, rpval);
    disp(mess)
end

% subfunction compute_pval() %
function [rpinf,lpinf]=compute_pval(rval, lval)
%load('/space/okapi/3/data/laurence/ADF/cutting_planes/Dice_cc_lrh_lta.mat'); %loads Dr Dl
%rh_stat_file='/space/okapi/3/data/laurence/ADF/cutting_planes/rh.CorpusCallosumCutDice.adf';
%lh_stat_file='/space/okapi/3/data/laurence/ADF/cutting_planes/lh.CorpusCallosumCutDice.adf';
%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end
rh_stat_file=strcat(fsafdDir, '/rh.CorpusCallosumCutDice.adf');
lh_stat_file=strcat(fsafdDir, '/lh.CorpusCallosumCutDice.adf');
fidrh=fopen(rh_stat_file);
fidlh=fopen(lh_stat_file);
if(fidrh==-1 | fidlh==-1)
    mess=sprintf('Could not find %s or %s', rh_stat_file, lh_stat_file);
    error(mess)
end
while(strfind(fgetl(fidrh), '#'))
    pos=ftell(fidrh);
end
fseek(fidrh, pos, 'bof');
Dr=fscanf(fidrh, '%g');
while(strfind(fgetl(fidlh), '#'))
    pos=ftell(fidlh);
end
fseek(fidlh, pos, 'bof');
Dl=fscanf(fidlh, '%g');
%Distribution %
pas=0.05;
x=0:pas:1;
[hr] = hist(Dr,x);
[hl] = hist(Dl,x);
pr = hr/sum(hr);
pl = hl/sum(hl);
dinfr=find(x<=rval);
dinfl=find(x<=lval);
xinfr=x(dinfr);
xinfl=x(dinfl);
pinfr=pr(1:length(xinfr));
pinfl=pl(1:length(xinfl));
if (rval>=0 & length(xinfr) >1 )
    rpinf=trapz(xinfr,pinfr)/pas;
elseif (rval>=0 & (length(xinfr)<2))
    pas2=pas/10;
    x2=0:pas2:1;
    [hr2] = hist(Dr,x2);
    pr2 = hr2/sum(hr2);
    dinfr2=find(x2<=rval);
    xinfr2=x2(dinfr2);
    pinfr2=pr2(1:length(xinfr2));
   if(length(xinfr2)>1)
        rpinf=trapz(xinfr2,pinfr2)/pas2;
    else
        rpinf=0;
    end 
else
    rpinf=0;
end
if (lval>=0 & length(xinfl) >1 )
    lpinf=trapz(xinfl,pinfl)/pas;
elseif (lval>=0 & (length(xinfl)<2))
    pas2=pas/10;
    x2=0:pas2:1;
    [hl2] = hist(Dl,x2);
    pl2 = hl2/sum(hl2);
    dinfl2=find(x2<=lval);
    xinfl2=x2(dinfl2);
    pinfl2=pl2(1:length(xinfl2));
   if(length(xinfl2)>1)
        lpinf=trapz(xinfl2,pinfl2)/pas2;
    else
        lpinf=0;
    end 
else
    lpinf=0;
end




    
    

        
    
  
   
