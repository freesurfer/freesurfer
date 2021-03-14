 function [ddr, ddl, Isubj]=cc_cut_dir_adf(dirname, th_pval)
% For each subject in the directory "dirname":
% Computes the Dice coefficients measuring the overlap
% of the WM volume in right and left hemispheres to check  
% if the corpus_callosum is correctly located. 
%
% Uses the p values
% 


%
% cc_cut_dir_afd.m
%
% Original Author: Laurence Wastiaux
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


if (nargin<1 | nargin>2)
    msg=sprintf('USAGE: [Dr,Dl]=cc_cut_dir_adf(dirname, th_pval)');
    disp(msg)
end

%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end

% Load Buckner stats %
% ------------------ %
%rh_stat_file='/space/okapi/3/data/laurence/ADF/cutting_planes/rh.CorpusCallosumCutDice.adf';
%lh_stat_file='/space/okapi/3/data/laurence/ADF/cutting_planes/lh.CorpusCallosumCutDice.adf';
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
DDr=fscanf(fidrh, '%g');
while(strfind(fgetl(fidlh), '#'))
    pos=ftell(fidlh);
end
fseek(fidlh, pos, 'bof');
DDl=fscanf(fidlh, '%g');

files=dir(dirname);
nf=files;
I=[];
flag=1;
R=0;
L=0;
for i=5:300
    flag=1;
    SubjectDir=strcat(dirname,'/',files(i).name);
    subject_name=files(i).name;
    
    %% Load Talairach transform matrix %%
    %%---------------------------------%%
    
    ttfile=strcat(SubjectDir,'/mri/transforms/talairach.lta');
    fid=fopen(ttfile, 'r');
    if (fid ~= -1) 
        while feof(fid) == 0
            linef=fgetl(fid);
            if (strfind(linef, '1 4 4'))
                pos=ftell(fid);
                break
            end    
        end
        A=(fscanf(fid, '%g',12))';
        Tal_transform =(reshape(A, [4 3]))';
    else
        flag=0;
    end
    %%  Load the wm & aseg volume    %%
    %% ----------------------------  %%

    SDir=strcat(SubjectDir,'/');
    CorDir1=strcat(SubjectDir,'/mri/wm/');
    CorDir2=strcat(SubjectDir,'/mri/aseg/');
    d1=dir(CorDir1);
    d2=dir(CorDir2);

    if (length(d1)<3 | (strcmp(files(i).name,'010611_vc7044')==1) | ( length(strfind(files(i).name,'0'))==0 &&(length(strfind(files(i).name,'1'))==0 )))
        mess1=sprintf('Cannot find the coronals of the volume wm for the subject %s',SubjectDir);
        disp(mess1)
        flag=0;
    elseif (length(d2)<3)
        mess2=sprintf('Cannot find the coronals of the volume aseg for the subject %s',SubjectDir);
        disp(mess2)
        flag=0;
    end
    if (flag~=0)
    
        [vol1 M1 t1]=load_cor2(SDir,'wm');
        [vol2 M2 t2]=load_cor2(SDir,'aseg');
        if ((t1~=1) | (t2~=1))
            flag=0;
        else   

            %% Read position of the cutting plane  %%
            %% Assumes that a log file is present  %%
            %% ----------------------------------  %%

            logfilename=strcat('/autofs/homes/001/wastiaux/adf/cutting_planes/logfiles_neo2/',files(i).name,'/ponscc.log'); % temporary location for the patients of the Buckner the training set...
            %logfilename=strcat(SDir, '/mri/ponscc.log'); % normal location
            fp=fopen(logfilename, 'r'); 
            if (fp==-1)
                msg=sprintf('Cannot open file %s', logfilename);
                disp(msg)
                flag=0;
            else
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
                V_cc=[x_cc ; y_cc ; z_cc; 1];
                V_cc_tal= Tal_transform*V_cc;
            end
        end
    end
    
    if(flag~=0)
        %%    Compute statistics   %%
        %% ----------------------- %%
        rh_wm_1=0; lh_wm_1=0;
        rh_wm_2=0; lh_wm_2=0;
        rh_wm_12=0; lh_wm_12=0;

        if ( (size(vol1)) ~= (size(vol2)) )
            msg=sprintf('Wm and aseg volumes do not have the same size');
            disp(msg);
            flag=0;
        else
            sz=size(vol1);
            for u=1:sz(1)
                for j=1:sz(2)
                    for k=1:sz(3)
                        if (vol2(u,j,k)==41 | vol2(u,j,k)==46) % the cerrebellum wm is taken into account 
                            rh_wm_2 = rh_wm_2 +1;
                        end
                        if (vol2(u,j,k)==2 | vol2(u,j,k)==7)
                            lh_wm_2 = lh_wm_2 +1;
                        end
                        if (vol1(u,j,k)~=0)
                            v=[j u k 1]';
                            v_tal=Tal_transform*v;
                            if (v_tal(1) < V_cc(1))
                                rh_wm_1 = rh_wm_1 +1;
                            else
                                lh_wm_1 = lh_wm_1 +1;
                            end
                        end
                        if (vol1(u,j,k)~=0 && (vol2(u,j,k)==41 | vol2(u,j,k)==46))
                            vr=[j u k 1]';
                            vr_tal=Tal_transform*vr;
                            if (vr_tal(1) < V_cc(1))
                                rh_wm_12 = rh_wm_12 +1;
                            end
                        end
                        if (vol1(u,j,k)~=0 && (vol2(u,j,k)==2 | vol2(u,j,k)==7))
                            vl=[j u k 1]';
                            vl_tal=Tal_transform*vl;
                            if (vl_tal(1) > V_cc(1))
                                lh_wm_12 = lh_wm_12 +1;
                            end
                        end
                    end 
                end
            end
            
            %% Compute Dice coefficients  %%
            %% -------------------------  %%
           
            dr= 2*rh_wm_12/(rh_wm_1+rh_wm_2);
            R=[R dr];
            nf(length(R)-1).name=subject_name;
            dl= 2*lh_wm_12/(lh_wm_1+lh_wm_2);
            L=[L dl];
            Isubj=[Isubj i];
            [rpval lpval]=compute_pval(dr ,dl, DDr, DDl);
            if(rpval<th_pval | lpval<th_pval)
                mess=sprintf('Wrong cc location for %s, (dr=%g, rpval=%.3g, dl=%g, lpval=%.3g)', files(i).name,dr, rpval, dl, lpval);
                disp(mess)
            else
                mess2=sprintf('cc location OK for %s, (dr=%g, rpval=%.3g, dl=%g, lpval=%.3g)', files(i).name,dr, rpval, dl, lpval);
                disp(mess2)
            end
        end
        ddr=R;
        ddl=L;
    else
        i=i+1;
    end 
end
ddr=ddr(2:end);
ddl=ddl(2:end);

% subfunction compute_pval() %
function [rpinf,lpinf]=compute_pval(rval, lval, Dr, Dl)
%load('/space/okapi/3/data/laurence/ADF/cutting_planes/Dice_cc_rh_lh.mat'); %loads Dr Dl
%load('/space/okapi/3/data/laurence/ADF/cutting_planes/Dice_cc_lrh_tal.mat'); %loads R L
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

        
    
  
   
