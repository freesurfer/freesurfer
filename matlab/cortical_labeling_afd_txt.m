function [A_lh, A_rh]=cortical_label_adf(subject, p_val)
% Computes the area of the different cortical labels 
% and compare them to the normal range
% Uses p_value to detect the abnormal areas
% Uses the lh/rh.parc.txt files
%


%
% cortical_labeling_afd_txt.m
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


if (nargin<2 | nargin>2)
    msg=sprintf('USAGE: [L, R]=cortical_label_adf(subject, p_val)');
    disp(msg)
end

A_lh=zeros(1,84);
A_rh=zeros(1,84);
namerh=[];
namelh=[];

labmapfile=('/space/lyon/1/fsdev/freesurfer_dev/Simple_surface_labels2002.txt');
[label name val1 val2 val3 val4]=textread(labmapfile,'%d %s %d %d %d %d',85);


%% a- Check if the left and right stats files exist %%
LabelDir=strcat(subject, '/label');
dlabel=dir(LabelDir);
if(length(dlabel)<=2)
    mess=sprintf('Directory %s is empty', LabelDir);
    error(mess)
end
if(length(dlabel)>2)
    if ((exist(strcat(LabelDir, '/rh.parc.txt'))==0) | (exist(strcat(LabelDir, '/lh.parc.txt'))==0))
        mess=sprintf('Cannot find the stats files for subject %s', subject);
        error(mess)
    end
end

%% b- Read areas for each label & for both hemispheres using the lh/rh.parc.txt files %%
sn=strread(subject, '%s', 'delimiter', '/');
subject_name=char(sn(length(sn)));
rh_LabelFile=strcat(subject, '/label/rh.parc.txt');
lh_LabelFile=strcat(subject, '/label/lh.parc.txt');

lh_fid=fopen(lh_LabelFile);
rh_fid=fopen(rh_LabelFile);
if(lh_fid == -1)
    msg=sprintf('Cannot open file %s', lh_LabelFile);
    error(msg)
elseif(rh_fid ==-1)
    msg=sprintf('Cannot open file %s', rh_LabelFile);
    error(msg)
else
    while(feof(lh_fid)==0)
        s=fgetl(lh_fid);
        if(strfind(s,'structure name')>1)
            s=fgetl(lh_fid);
            lhpt=ftell(lh_fid);
            break;
        end
    end
    fseek(lh_fid, lhpt,'bof');
    while(feof(lh_fid)==0)
        s=fgetl(lh_fid);
        a=strread(s,'%d',2);  %area
        aa=strread(s,'%s',11); %name
        llab=0;
        for test=1:length(name)
            if(strfind(char(name(test)),char(aa(11)))>=1)
                llab=label(test);
                break;
            end  
        end 
        if(llab==0)
            disp(sprintf('impossible to find the label')); 
        elseif(A_lh(llab)~=0)
            disp(sprintf('already attributed'));
        else
            A_lh(llab)=a(2);
        end 
    end
    while(feof(rh_fid)==0)
        s=fgetl(rh_fid);
        if(strfind(s,'structure name')>1)
            s=fgetl(rh_fid);
            rhpt=ftell(rh_fid);
            break;
        end
    end
    fseek(rh_fid, rhpt,'bof');
    while(feof(rh_fid)==0)
        s=fgetl(rh_fid);
        a=strread(s,'%d',2);
        aa=strread(s,'%s',11);
        llab=0;
        for test=1:length(name)
            if(strfind(char(name(test)),char(aa(11)))>=1)
                llab=label(test);
                break;
            end  
        end 
        if(llab==0)
            disp(sprintf('impossible to find the label')); 
        elseif(A_rh(llab)~=0)
            disp(sprintf('already attributed'));
        else
            A_rh(llab)=a(2);
        end 
    end    
end
nl=length(name)-1;
if(length(A_lh)< nl)
    for lleft=(length(A_lh)+1):nl
        A_lh(lleft)=0;
    end
end
if(length(A_rh)< nl)
    for lleft=(length(A_rh)+1):nl
        A_rh(lleft)=0;
    end
end
if (sum(A_lh)==0 | sum(A_rh)==0)
    msg=sprintf('Total area is null');
    error(msg)
else
    total_lharea=sum(A_lh);
    total_rharea=sum(A_rh);
end
for i=1:nl
    A_lh(i)= 100 * A_lh(i) / total_lharea;
    A_rh(i)= 100 * A_rh(i) / total_rharea;
end
% Compute p_values %

count=0;
for k =1:length(A_lh)
    %disp(k)
    
    [Rpval_inf,Rpval_sup, Lpval_inf, Lpval_sup]=compute_pval(k,A_rh(k),A_lh(k));
    if ((Lpval_inf<p_val | Lpval_sup <p_val) & (A_lh(k)~=0))
        fprintf('%s, left hemisphere: percent area of %s (%d) = %g (pval_inf=%g, pval_sup=%g)\n' ,subject_name, char(name(k+1)), k , A_lh(k),Lpval_inf, Lpval_sup);
        count=1;
    end
    if ((Rpval_inf<p_val | Rpval_sup <p_val)& (A_rh(k)~=0))
        fprintf('%s, right hemisphere: percent area of %s (%d) = %g (pval_inf=%g, pval_sup=%g)\n' ,subject_name, char(name(k+1)),k, A_rh(k),Rpval_inf,Rpval_sup);
        count=1;
    end
   
end
if(count==0)
    fprintf('%s: subject is normal\n', subject_name);
end


% subfunction compute_pval() %
function [pinf_rh, psup_rh, pinf_lh, psup_lh]=compute_pval(llabel,R,L)
%load('/space/okapi/3/data/laurence/ADF/cortical_labeling/PercentArea_labels.mat'); % loads D_lh and D_rh
% load('/space/okapi/3/data/laurence/ADF/cortical_labeling/PercentArea_labels_parctxt.mat');
% D_lh=D2_lh;
% D_rh=D2_rh;
%rh_stat_file='/space/okapi/3/data/laurence/ADF/cortical_labeling/rh.CorticalLabelingPercentArea_txt.adf';
%lh_stat_file='/space/okapi/3/data/laurence/ADF/cortical_labeling/lh.CorticalLabelingPercentArea_txt.adf';
%%% Get the table's directory %%%
if(getenv('FREESURFER_HOME'))
    fsh=getenv('FREESURFER_HOME');
    fsafdDir=strcat(fsh, '/fsafd');
else
    error(sprintf('Impossible to find FREESURFER_HOME\n'));
end
rh_stat_file=strcat(fsafdDir, '/rh.CorticalLabelingPercentArea_txt.adf');
lh_stat_file=strcat(fsafdDir, '/lh.CorticalLabelingPercentArea_txt.adf');
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
D_rhtmp=fscanf(fidrh, '%g');
nrowr=length(D_rhtmp)/84;
D_rh=(reshape(D_rhtmp, [84, nrowr]))';
while(strfind(fgetl(fidlh), '#'))
    pos=ftell(fidlh);
end
fseek(fidlh, pos, 'bof');
D_lhtmp=fscanf(fidlh, '%g');
nrow=length(D_lhtmp)/84;
D_lh=(reshape(D_lhtmp, [84,nrow]))';
%distribution%
pas=0.01;
x=0:pas:18;
[h1] = hist(D_lh(:,llabel),x);
[h2] = hist(D_rh(:,llabel),x);
%[h1] = hist(D2_lh(:,llabel),x);
%[h2] = hist(D2_rh(:,llabel),x);
p1 = h1/sum(h1);
p2 = h2/sum(h2);
d1inf=find(x<=L);
d2inf=find(x<=R);
d1sup=find(x>L);
d2sup=find(x>R);
x1inf=x(d1inf);
x1sup=x(d1sup);
x2inf=x(d2inf);
x2sup=x(d2sup);
p1inf=p1(1:length(x1inf));
p2inf=p2(1:length(x2inf));
p1sup=p1(length(x)-length(x1sup)+1:end);
p2sup=p2(length(x)-length(x2sup)+1:end);
if( L>=0 & length(x1inf)>1 & length(x1sup)>1)
    pinf_lh=trapz(x1inf,p1inf)/pas;
    psup_lh=trapz(x1sup,p1sup)/pas;
else   % No subdivisions of the last intervals
    pinf_lh=0;
    psup_lh=0;
end
if(R>=0 & length(x2inf)>1 & length(x2sup)>1)
    pinf_rh=trapz(x2inf,p2inf)/pas;
    psup_rh=trapz(x2sup,p2sup)/pas;
else
    pinf_rh=0;
    psup_rh=0;
end



  
