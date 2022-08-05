function Hu = MAP_hu_propose(Hu,Hg,tRL,bar_r,bar_x,phi,phi1,phi2,iFCg,iFCu,sigma1_sq,sigmaL_sq,bar_sigmaL_sq)
    [M,K,L1]=size(tRL);
    [~,L2]=size(bar_r);
    [N,~]=size(phi);
    %%
    DL=zeros(M,N,N);
    for n0=1:N
        theta_t=phi(:,n0);
        DL(:,:,n0)=Hg.'*diag(theta_t);
    end
    Eta_temp1=1/sigma1_sq*DL(:,:,1)'*DL(:,:,1);
    for l0=2:L1
        Eta_temp1=Eta_temp1+1/sigmaL_sq*DL(:,:,l0)'*DL(:,:,l0);
    end
    Eta_temp2=0;
    for l0=(L1+1):N
        Eta_temp2=Eta_temp2+1/bar_sigmaL_sq*DL(:,:,l0)'*DL(:,:,l0);
    end
    Eta=zeros(N,N,K);
    for k0=1:K
        Etmp=iFCu(:,:,k0)+Eta_temp1+Eta_temp2*abs(bar_x(1))^2; 
        Eta(:,:,k0)=Etmp\eye(N);
    end
    f0=obj_likelihood(Hu,DL,tRL,bar_r,bar_x,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
    while(1)
        Hu_old=Hu;
        for k0=1:K
            nu=1/sigma1_sq*DL(:,:,1)'*tRL(:,k0,1);
            for l0=2:L1
                nu=nu+1/sigmaL_sq*DL(:,:,l0)'*tRL(:,k0,l0);
            end
            for l0=(L1+1):N
                bar_l0=l0-L1;
                tmp1=0;
                Dt=DL(:,:,l0);
                for ii0=1:K
                    if ii0~=k0
                        tmp1=tmp1+bar_x(ii0)*Dt*Hu(:,ii0);
                    end
                end
                nu=nu+1/bar_sigmaL_sq*(bar_x(k0)*Dt)'*(bar_r(:,bar_l0)-tmp1);
            end
    %         Eta=iFCu(:,:,k0)+Eta_temp1+Eta_temp2*abs(bar_x(k0))^2;       
    %         Hu(:,k0)=Eta\nu;
%             Hu(:,k0)=Eta*nu;
            Hu(:,k0)=Eta(:,:,k0)*nu;
        end
        f1=obj_likelihood(Hu,DL,tRL,bar_r,bar_x,sigma1_sq,sigmaL_sq,bar_sigmaL_sq);
        if f1<f0
            Hu=Hu_old;
            fprintf('err HU\n');
            break
        elseif abs(log(abs(f1-f0)/abs(f0)))>3  %abs(f1-f0)<1e-6
            break
        end
        f0=f1;
    end
end

function [fL] = obj_likelihood(Hu,DL,tRL,bar_r,bar_x,sigma1_sq,sigmaL_sq,bar_sigmaL_sq)
    [~,K,L1]=size(tRL);
    [N,~]=size(Hu);
    %% prior
    %% likelihood
    fL=0;
    for l0=1:L1
        tmp=0;
        for k0=1:K
            tmpi=tRL(:,k0,l0)-DL(:,:,l0)*Hu(:,k0);
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
        tmp=bar_r(:,bar_l0)-DL(:,:,l0)*Hu*bar_x;
        tmp=-sum(abs(tmp(:)).^2)/bar_sigmaL_sq;
        fL=fL+tmp;
    end
end

