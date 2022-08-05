function Hu = BALS_Hu_map(R,Hg,phi_t,iC_IRS,sigma_h)
    [M,K,tL]=size(R);
    [N,~]=size(Hg);
    %%
    %%
    iFCu=zeros(N,N,K);
    for k0=1:K
        tmp=zeros(N,N);
        for m0=1:M
            hg=diag(Hg(:,m0));
            Cirs=iC_IRS(:,:,m0,k0);
            tmp=tmp+hg'*Cirs*hg;
        end
        iFCu(:,:,k0)=tmp;
    end
    iFCu=iFCu.*sigma_h;
    %%
    D_w=zeros(M,N,tL);
    for l0=1:tL
        t_theta=phi_t(:,l0);
        D_w(:,:,l0)=Hg.'*diag(t_theta);
    end
    A=zeros(N,N);
    B=zeros(N,K);
    for l0=1:tL
        A=A+D_w(:,:,l0)'*D_w(:,:,l0);
        B=B+D_w(:,:,l0)'*R(:,:,l0);
    end
    Hu=zeros(N,K);
    for k0=1:K
        Eta=A+iFCu(:,:,k0);
        Hu(:,k0)=Eta\B(:,k0);
    end
end

