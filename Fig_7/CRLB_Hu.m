function Fu = CRLB_Hu(bar_x,phi,iECu,ECg,sigma1_sq,sigmaL_sq,bar_sigmaL_sq,K,L1)
    [N,~]=size(phi);
    %%
    S_Cg=sum(ECg,3);
    E_DL=zeros(N,N,N);
    for n0=1:N
        theta_t=phi(:,n0);
        E_DL(:,:,n0)=diag(theta_t)'*S_Cg*diag(theta_t);
    end
%     DL=zeros(M,N,N);
%     for n0=1:N
%         theta_t=phi(:,n0);
%         DL(:,:,n0)=Hg.'*diag(theta_t);
%     end
    A1=1/sigma1_sq*E_DL(:,:,n0);
    for l0=2:L1
        A1=A1+1/sigmaL_sq*E_DL(:,:,n0);
    end
    B1=0;
    for l0=(L1+1):N
        B1=B1+1/bar_sigmaL_sq*E_DL(:,:,n0);
    end
    C=bar_x*bar_x';
    %%
    iFCu=iECu;
    %%
    IFC=cell(K,1);
    for k0=1:K
        IFC{k0}=iFCu(:,:,k0);
    end
    IFC=blkdiag(IFC{:});
    IFC=IFC+kron( eye(K),A1 )+kron(C.',B1); 
    IFC=IFC-diag(IFC)+real(diag(IFC));
    Fu=(IFC+IFC')/2;
end