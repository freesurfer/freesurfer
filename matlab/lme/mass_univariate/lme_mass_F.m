function fstats = lme_mass_F(stats,CM,prs)
% fstats = lme_mass_F(stats,CM,prs) 
% 
% Inference for the fixed-effects in the linear mixed-effects model. Allows 
% to test a different model at each voxel/vertex (Depends on the Statistics
% toolbox).
%
% Input
% stats: Structure array containing statistiscs for every voxel/vertex (see
% lme_FSfit for more details on these statistics).
% CM: Structure array containing a contrast matrix CM.C for each and every 
% voxel/vertex. If length(CM)==1 the same contrast matrix CM.C is tested at
% each and every voxel/vertex.
%
% Output
% fstats.F: F-Statistic.
% fstats.pval: P-value of the F-Statistic.
% fstats.sgn: Sign of the contrast.
% fstats.df: Degrees of freedom of the F-Statistic.
%
% $Revision: 1.1.2.2 $  $Date: 2013/02/23 21:08:10 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:08:10 $
%    $Revision: 1.1.2.2 $
% References: Bernal-Rusiel J.L., Greve D.N., Reuter M., Fischl B., Sabuncu
% M.R., 2012. Statistical Analysis of Longitudinal Neuroimage Data with Linear 
% Mixed Effects Models, NeuroImage, doi:10.1016/j.neuroimage.2012.10.065.
%
tic;
if nargin < 2
    error('Too few inputs');
elseif nargin < 3
    prs = 8;
end;
nv = length(stats);
nCM = length(CM);
if nCM ~= nv
    if nCM == 1
        CM(1:nv) = struct('C',CM.C);
    else
        error(['The number of elements (contrasts) in CM must be equal to the'...
            ' number of locations (length(stats))']);
    end
end;
for i=1:nv
    if (stats(i).lreml ~= -10^10)  &&  (size(CM(i).C,2) ~= length(stats(i).Bhat))
        error(['The number of colums in the contrast matrix CM('...
            num2str(i) ').C is different from the corresponding number'...
            ' of fixed effects in stats(' num2str(i) ').Bhat']);
    end;
end;
F = zeros(1,nv);
pval = ones(1,nv);
sgn = zeros(1,nv);
df = zeros(2,nv);
if (prs==1) || (matlabpool('size') ~= prs)
    if (matlabpool('size') > 0)
        matlabpool close;
    end;
    if (prs>1)
        matlabpool(prs);
    end;
end;
display(' ');
display('Computing F statistics, degrees of freedom and P-values ...');
parfor i=1:nv
    if stats(i).lreml ~= -10^10
        Bhat = stats(i).Bhat;
        C = CM(i).C;
        Zcols = stats(i).Zcols;
        q = length(Zcols);
        nth = q*(q+1)/2+1;
        invEI = stats(i).invEI;
        Pth = stats(i).Pth;
        Qthth = stats(i).Qthth;
        %Estimation of the bias in the covariance matrix and computation of the
        %F-statistic and the degrees of freedom of the test.
        Bias = 0;
        CBhat = stats(i).CovBhat;
        OM = C'*(C*CBhat*C')^-1*C;
        A1 = 0; A2 = 0;
        Term1 = OM*CBhat;
        Term2 = CBhat*OM*CBhat;
        for k=1:nth
            Pk = squeeze(Pth(k,:,:));
            for j=1:nth
                Qkj = squeeze(Qthth(k,j,:,:));
                Pj = squeeze(Pth(j,:,:));
                Bias = Bias + invEI(k,j)*(Qkj-Pk*CBhat*Pj);
                A1 = A1+invEI(k,j)*trace(Term1*Pk*CBhat)*trace(Term1*Pj*CBhat);
                A2 = A2+invEI(k,j)*trace(Term1*Pk*Term2*Pj*CBhat);
            end;
        end;
        szC = size(C,1);
        Bdf = (A1+6*A2)/(2*szC);
        g = ((szC+1)*A1-(szC+4)*A2)/((szC+2)*A2);
        d = 3*szC+2*(1-g);
        c1 = g/d;
        c2 = (szC-g)/d;
        c3 = (szC+2-g)/d;
        EF = (1-A2/szC)^-1;
        VF = (2/szC)*(1+c1*Bdf)/((1-c2*Bdf)^2*(1-c3*Bdf));
        ro = VF/(2*EF^2);
        m = 4+(szC+2)/(szC*ro-1);
        l = m/(EF*(m-2));
        Bias = CBhat*Bias*CBhat;
        CovBhat = CBhat + 2*Bias;
        Fstat = l*Bhat'*C'*(C*CovBhat*C')^-1*C*Bhat/szC;
        if Fstat<0
           Fstat = 0; 
        end;
        F(i) = Fstat;
        pval(i) = 1-fcdf(Fstat,szC,m);
        contrast = C*Bhat;
        sgn(i) = sign(contrast(1));
        df(:,i) = [szC;m];
    end;
end;
fstats.F = F;
fstats.pval = pval;
fstats.sgn = sgn;
fstats.df = df;
if (matlabpool('size') > 0)
    matlabpool close;
end;
toc;