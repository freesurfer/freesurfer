function [proba, pinf]=talmat(filename, th_pval, DirTable)
%
% Computes the probability of the Talairach transform matrix
%       'filename/mri/transforms/talairach.xfm'.  
%  Uses the mean vector and covariance matrix obtained with talairachin_table.m from 
%     the data set (default data set: /space/neo/2/recon/buckner)
%  Uses th_pval as a threshold for the p-values to detect the unlikely transform matrices
%


%
% talairaching_afd.m
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



if (nargin<2 | nargin>3)
    msg=sprintf('USAGE: [proba, pinf]=talmat(Subject, th_pval, <DirTable>)');
    disp(msg)
end

if(nargin==3)
    fsafdDir=DirTable;
else
    %%% Get the tables'directory %%%
    if(getenv('FREESURFER_HOME'))
        fsh=getenv('FREESURFER_HOME');
        fsafdDir=strcat(fsh, '/fsafd');
    else
        error(sprintf('Impossible to find FREESURFER_HOME\n'));
    end
end

%%% 1x9 mean vector obtained from the training set %%%
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_mean.mat'); %Loads mu
%file_mu='/space/okapi/3/data/laurence/ADF/talairaching/TalairachingMean.adf'; %Loads mu
file_mu=strcat(fsafdDir, '/TalairachingMean.adf');
fi=fopen(file_mu);
pos=0;
if(fi==-1)
    mess=sprintf('Could not find %s', file_mu);
    error(mess)
else
    while(strfind(fgetl(fi), '#'))  % skip the header
        pos=ftell(fi);
    end
    fseek(fi, pos, 'bof');
    mu=(fscanf(fi, '%g'))';
    fclose(fi);
end

%%% 9x9 covariance matrix obtained from the training set %%%
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_regularizedCov2.mat'); %loads sigma
%sigma_file='/space/okapi/3/data/laurence/ADF/talairaching/TalairachingCovariance.adf';
sigma_file=strcat(fsafdDir, '/TalairachingCovariance.adf');
fis=fopen(sigma_file);
pos=0;
if(fis==-1)
    mess=sprintf('Could not find %s', sigma_file);
    error(mess)
else
    while(strfind(fgetl(fis), '#'))  % skip the header
        pos=ftell(fis);
    end
    fseek(fis, pos, 'bof');
    sig=fscanf(fis, '%g');
    sigma=reshape(sig, [9,9]);
    fclose(fis);
end

matname=strcat(filename, '/mri/transforms/talairach.xfm');
fid=fopen(matname, 'r');
if fid == -1
    message=sprintf('Cannot open the file %s',matname);
    disp(message);
else
    while feof(fid) == 0
        linef=fgetl(fid);
        nb=findstr(linef, 'Linear_Transform');
        if nb == 1
            pos=ftell(fid);
            break
        end    
    end
    A=(fscanf(fid, '%g',12))';
    A=[A(1) A(2) A(3) A(5) A(6) A(7) A(9) A(10) A(11)]; % supression of the translation parameters
    status=fclose(fid);
end
proba=mvnpdf(A,mu,sigma);
[pinf]=compute_pval(proba, fsafdDir);
if ( (pinf < th_pval) )
    mess1=sprintf('Talairach Transform: failed (p=%g, pval=%g)', proba, pinf);
    %mess2=sprintf('Talairach Transform: failed');
    disp(mess1)
else
    mess2=sprintf('Talairach Transform: OK (p=%g, pval=%g)', proba, pinf);
    %mess2=sprintf('Talairach Transform: OK');
    disp(mess2)
end


function [p_inf]=compute_pval(val, tableDir)
%load('/space/okapi/3/data/laurence/ADF/talairaching/transfo_param_probas.mat'); %loads y
%stat_file='/space/okapi/3/data/laurence/ADF/talairaching/TalairachingProbas.adf';
stat_file=strcat(tableDir, '/TalairachingProbas.adf');
fid=fopen(stat_file);
pos=0;
if(fid==-1)
    mess=sprintf('Could not find %s', stat_file);
    error(mess)
else
    while(strfind(fgetl(fid), '#'))
        pos=ftell(fid);
    end
    fseek(fid, pos, 'bof');
    y=fscanf(fid, '%g');
end
pas=0.05;
x=0:pas:1;
[h] = hist(y,x);
p = h/sum(h);
dinf=find(x<=val);
xinf=x(dinf);
pinf=p(1:length(xinf));
if (val>=0 & length(xinf) >1 )
    p_inf=trapz(xinf,pinf)/pas;
elseif (val>=0 & (length(xinf)<2 ))
    pas2=pas/5;
    x2=0:pas2:1;
    [h2] = hist(y,x2);
    p2 = h2/sum(h2);
    dinf2=find(x2<=val);
    xinf2=x2(dinf2);
    pinf2=p2(1:length(xinf2));
    p_inf=trapz(xinf2,pinf2)/pas2;
    p_sup=trapz(xsup2,psup2)/pas2;
else
    p_inf=0;
end
