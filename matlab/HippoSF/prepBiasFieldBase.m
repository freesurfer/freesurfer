function PSI=prepBiasFieldBase(siz,biasFieldOrder)

% Prepare basis functions for the bias field, and initialize weights
% disp('Preparing basis functions for bias field')
[X Y Z]=ndgrid(1:siz(1),1:siz(2),1:siz(3));
X=X-min(X(:)); X=X/max(X(:)); X=2*X-1;
Y=Y-min(Y(:)); Y=Y/max(Y(:)); Y=2*Y-1;
Z=Z-min(Z(:)); Z=Z/max(Z(:)); Z=2*Z-1;
orders=[];
for o=0:biasFieldOrder
    for x=0:o
        for y=0:o
            for z=0:o
                if x+y+z==o
                    orders{end+1}=[x y z];
                end
            end
        end
    end
end
PSI=zeros([siz length(orders)],'single');
for d=1:length(orders)
    psi=ones(siz);
    for x=1:orders{d}(1)
        psi=psi.*X;
    end
    for y=1:orders{d}(2)
        psi=psi.*Y;
    end
    for z=1:orders{d}(3)
        psi=psi.*Z;
    end
    PSI(:,:,:,d)=psi;
end