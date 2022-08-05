function Hg = ML_hg_propose(Hu,Hg,tRL,bar_r,bar_x,phi,sigma1_sq,sigmaL_sq,bar_sigmaL_sq)
    [M,K,L1]=size(tRL);
    [~,L2]=size(bar_r);
    [N,~]=size(phi);
    %%
    B=zeros(N,N,K);
    for k0=1:K
        for l0=1:N
            theta=phi(:,l0);
            B(:,l0,k0)=diag(theta)*Hu(:,k0);
        end
    end
    tmp_eta=0;
    for k0=1:K
        for l0=1:L1
            if l0==1
                tmp_eta=tmp_eta+1/sigma1_sq*conj(B(:,l0,k0))*B(:,l0,k0).';
            else
                tmp_eta=tmp_eta+1/sigmaL_sq*conj(B(:,l0,k0))*B(:,l0,k0).';
            end
        end
    end
    tmp_eta2=0;
    for l0=(1+L1):N
        tmp=0;
        for k0=1:K
            tmp=tmp+bar_x(k0)*B(:,l0,k0).';
        end
        tmp_eta2=tmp_eta2+1/bar_sigmaL_sq*(tmp'*tmp);
    end
    for m0=1:M
        nu=0;
        for k0=1:K
            for l0=1:L1
                if l0==1
                    nu=nu+1/sigma1_sq*tRL(m0,k0,l0)*conj(B(:,l0,k0));
                else
                    nu=nu+1/sigmaL_sq*tRL(m0,k0,l0)*conj(B(:,l0,k0));
                end
            end
        end
        for l0=(L1+1):N
            bar_l0=l0-L1;
            tmp=0;
            for k0=1:K
                tmp=tmp+conj(bar_x(k0)*B(:,l0,k0));
            end
            nu=nu+1/bar_sigmaL_sq*bar_r(m0,bar_l0)*tmp;
        end
        Eta=tmp_eta+tmp_eta2;
        Hg(:,m0)=Eta\nu;
    end
end

