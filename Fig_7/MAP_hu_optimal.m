function Hu = MAP_hu_optimal(Hg,tRL,bar_r,bar_x,phi,iC_IRS,sigma1_sq,sigmaL_sq,bar_sigmaL_sq)
    [M,K,L1]=size(tRL);
%     [~,L2]=size(bar_r);
    [N,~]=size(phi);
    %%
    DL=zeros(M,N,N);
    for n0=1:N
        theta_t=phi(:,n0);
        DL(:,:,n0)=Hg.'*diag(theta_t);
    end
    A1=1/sigma1_sq*DL(:,:,1)'*DL(:,:,1);
    for l0=2:L1
        A1=A1+1/sigmaL_sq*DL(:,:,l0)'*DL(:,:,l0);
    end
    B1=0;
    for l0=(L1+1):N
        B1=B1+1/bar_sigmaL_sq*DL(:,:,l0)'*DL(:,:,l0);
    end
    C=bar_x*bar_x';
    nu=1/sigma1_sq*DL(:,:,1)'*tRL(:,:,1);
    for l0=2:L1
        nu=nu+1/sigmaL_sq*DL(:,:,l0)'*tRL(:,:,l0);
    end
    for l0=(L1+1):N
        bar_l0=l0-L1;
        nu=nu+1/bar_sigmaL_sq*DL(:,:,l0)'*bar_r(:,bar_l0)*bar_x';
    end
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
    IFC=cell(K,1);
    for k0=1:K
        IFC{k0}=iFCu(:,:,k0);
    end
    IFC=blkdiag(IFC{:});
    IFC=IFC+kron( eye(K),A1 )+kron(C.',B1); 
    IFC=IFC-diag(IFC)+real(diag(IFC));
    IFC=(IFC+IFC')/2;
%         issymmetric(IFC)
%         tic
    x=IFC\nu(:);
%         toc
    Hu = reshape(x,[N,K]);
end