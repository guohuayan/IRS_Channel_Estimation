function Fg = CRLB_Hg(phi,iECg,ECu,sigma1_sq,sigmaL_sq,bar_sigmaL_sq,M,K,L1)
%     [M,K,L1]=size(tRL);
    [N,~]=size(phi);
    %%
    C_B=ECu;
%     B=zeros(N,N,K);
%     for k0=1:K
%         for l0=1:N
%             theta=phi(:,l0);
%             B(:,l0,k0)=diag(theta)*Hu(:,k0);
%         end
%     end
    tmp_eta=0;
    for k0=1:K
        for l0=1:L1
            if l0==1
                tmp_eta=tmp_eta+1/sigma1_sq*C_B(:,:,k0);
            else
                tmp_eta=tmp_eta+1/sigmaL_sq*C_B(:,:,k0);
            end
        end
    end
    tmp_eta2=0;
    for l0=(1+L1):N
        tmp=sum(C_B,3);
        tmp_eta2=tmp_eta2+1/bar_sigmaL_sq*tmp;
    end
    iFCg=iECg;
    %%
    Fg=zeros(N,N,M);
    for m0=1:M
        Eta=iFCg(:,:,m0)+tmp_eta+tmp_eta2;
        Fg(:,:,m0)=Eta;
    end
end

