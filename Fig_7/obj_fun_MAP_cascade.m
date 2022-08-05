function [f,fL,fpr] = obj_fun_MAP_cascade(Hu,Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq)
    [M,K,L1]=size(tRL);
%     [~,L2]=size(bar_r);
%     [N,~,~]=size(iFCu);
    [N,~]=size(phi);
    %%
    MAP_Hi=zeros(M,N,K);
    for k0=1:K
        MAP_Hi(:,:,k0)=Hg.'*diag(Hu(:,k0));
    end
    MAP_Hi=permute(MAP_Hi,[2,1,3]);
    fpr=0;
    for m0=1:M
        for k0=1:K
            h=MAP_Hi(:,m0,k0);
            Cirs=iC_IRS(:,:,m0,k0);
            fpr=fpr-h'*Cirs*h;
        end
    end
    fpr=real(fpr);
    %% prior
%     fpr=0;
%     for m0=1:M
%         fpr=fpr-Hg(:,m0)'*iFCg(:,:,m0)*Hg(:,m0);
%     end
%     for k0=1:K
%         fpr=fpr-Hu(:,k0)'*iFCu(:,:,k0)*Hu(:,k0);
%     end
%     fpr=real(fpr);
    %% likelihood
    fL=0;
    for l0=1:L1
        tmp=0;
        theta=phi(:,l0);
        for k0=1:K
            tmpi=tRL(:,k0,l0)-Hg.'*diag(theta)*Hu(:,k0);
            tmp=tmp-sum(abs(tmpi(:)).^2);
        end
        if l0==1
            tmp=tmp/sigma1_sq;
        else
            tmp=tmp/sigmaL_sq;
        end
        fL=fL+tmp;
    end
    for l0=(L1+1):N
        bar_l0=l0-L1;
        theta=phi(:,l0);
        tmp=bar_r(:,bar_l0)-Hg.'*diag(theta)*Hu*bar_x;
%         for k0=1:K
%             tmp=tmp-bar_x(k0)*Hg.'*diag(theta)*Hu(:,k0);
%         end
        tmp=-sum(abs(tmp(:)).^2)/bar_sigmaL_sq;
        fL=fL+tmp;
    end
    f=fL+fpr;
end