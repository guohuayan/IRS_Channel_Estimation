function tmpH = MAP_hi_propose_cascade(Hi,Aw,B1,B0,bar_x,iC_IRS,phi,L1,bar_sigmaL_sq,k0,xi)
    [M,N,K]=size(Hi);
    %%
    temp_par=0;
    for i0=1:K
        temp_par=temp_par+Hi(:,:,i0).*bar_x(i0);
    end
    B2=zeros(M,N);
    B3=zeros(N,N);
    for l0=(L1+1):N
        theta=phi(:,l0);
        tmp=-temp_par*(theta*theta');
        tmp=tmp/bar_sigmaL_sq;
        B2=B2+tmp;
        tmp=theta*theta'/bar_sigmaL_sq;
        B3=B3+tmp;
    end 
    %%
%     Ck0=iC_IRS(:,:,:,k0);
    B_tmp=B1(:,:,k0)+B0+B2*conj(bar_x(k0))+Hi(:,:,k0)*B3+xi.*Hi(:,:,k0);
%     A=A.';
    B_tmp=B_tmp.';
    %%
    tmpH=zeros(N,M);
    for m0=1:M
%         Ckm=Ck0(:,:,m0);
%         tmpH(:,m0)=(A+Ckm)\B_tmp(:,m0);
        Atmp=Aw(:,:,m0,k0);
        tmpH(:,m0)=Atmp*B_tmp(:,m0);
    end
    tmpH=tmpH.';
    %%
end