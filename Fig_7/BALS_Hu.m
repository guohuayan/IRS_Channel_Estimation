function Hu = BALS_Hu(R,Hg,phi_t,Hu_old)
    [M,K,tL]=size(R);
    [~,N]=size(Hg);
    A_w=zeros(M,N,tL);
    for l0=1:tL
        t_theta=phi_t(:,l0);
        A_w(:,:,l0)=Hg*diag(t_theta);
    end
    A=zeros(N,N);
    B=zeros(N,K);
    for l0=1:tL
        A=A+A_w(:,:,l0)'*A_w(:,:,l0);
        B=B+A_w(:,:,l0)'*R(:,:,l0)';
    end
    if rank(A)==N
        Hu=A\B;
    else
        Hu=Hu_old;
    end
end

