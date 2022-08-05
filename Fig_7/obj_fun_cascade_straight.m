function [f,fL,fpr] = obj_fun_cascade_straight(Hi,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq)
    [M,K,L1]=size(tRL);
    [N,~]=size(phi);
    %%
    Hi=permute(Hi,[2,1,3]);  %% N M K
    fpr=0;
    for m0=1:M
        for k0=1:K
            h=Hi(:,m0,k0);
            Cirs=iC_IRS(:,:,m0,k0);
            fpr=fpr-h'*Cirs*h;
        end
    end
    fpr=real(fpr);
    Hi=permute(Hi,[2,1,3]); %% M N K
    %% prior
    %% likelihood
    fL=0;
    for l0=1:L1
        tmp=0;
        theta=phi(:,l0);
        for k0=1:K
            tmpi=tRL(:,k0,l0)-Hi(:,:,k0)*theta;
            tmp=tmp-sum(abs(tmpi(:)).^2);
        end
        if l0==1
            tmp=tmp/sigma1_sq;
        else
            tmp=tmp/sigmaL_sq;
        end
        fL=fL+tmp;
    end
    A=0;
    for k0=1:K
        A=A+Hi(:,:,k0).*bar_x(k0);
    end
    for l0=(L1+1):N
        bar_l0=l0-L1;
        theta=phi(:,l0);
        tmp=bar_r(:,bar_l0)-A*theta;
        tmp=-sum(abs(tmp(:)).^2)/bar_sigmaL_sq;
        fL=fL+tmp;
    end
    f=fL+fpr;
end