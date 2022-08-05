function Hg = BALS_Hg_map(R,Hu,phi_t,iC_IRS,sigma_h)
    [M,K,tL]=size(R);
    [N,~]=size(Hu);
    %%
    iFCg=zeros(N,N,M);
    for m0=1:M
        tmp=zeros(N,N);
        for k0=1:K
            hu=diag(Hu(:,k0));
            Cirs=iC_IRS(:,:,m0,k0);
            tmp=tmp+hu'*Cirs*hu;
        end
        iFCg(:,:,m0)=tmp;
    end
    iFCg=iFCg.*sigma_h;
    %%
    D_w=zeros(K,N,tL);
    for l0=1:tL
        t_theta=phi_t(:,l0);
        D_w(:,:,l0)=Hu.'*diag(t_theta);
    end
    A=zeros(N,N);
    B=zeros(N,M);
    for l0=1:tL
        A=A+D_w(:,:,l0)'*D_w(:,:,l0);
        B=B+D_w(:,:,l0)'*R(:,:,l0).';
    end
    Hg=zeros(N,M);
    for m0=1:M
        Eta=A+iFCg(:,:,m0);
        Hg(:,m0)=Eta\B(:,m0);
    end
end

