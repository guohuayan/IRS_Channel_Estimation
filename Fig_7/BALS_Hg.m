function Hg = BALS_Hg(R,Hu,phi_t,Hg_old)
    [M,K,tL]=size(R);
    [N,~]=size(Hu);
    D_w=zeros(N,K,tL);
    for l0=1:tL
        t_theta=phi_t(:,l0);
        D_w(:,:,l0)=diag(t_theta)*Hu;
    end
    A=zeros(N,N);
    B=zeros(M,N);
    for l0=1:tL
        A=A+D_w(:,:,l0)*D_w(:,:,l0)';
        B=B+R(:,:,l0)*D_w(:,:,l0)';
    end
    if rank(A)==N
        Hg=B/A;
    else
        Hg=Hg_old;
    end
end

