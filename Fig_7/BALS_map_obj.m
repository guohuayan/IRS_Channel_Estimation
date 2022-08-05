function f = BALS_map_obj(R,Hg,Hu,phi_t,iC_IRS,sigma_h)
    [M,K,tL]=size(R);
    [N,~]=size(Hu);
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
    %%
    f=0;
    for l0=1:tL
        t_theta=phi_t(:,l0);
        tmp=R(:,:,l0)-Hg.'*diag(t_theta)*Hu;
        err=sum(abs(tmp(:)).^2);
        f=f+err;
    end
    f=-f+fpr*sigma_h;
end

