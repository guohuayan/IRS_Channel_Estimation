function B = Khatri_Rao(G,Phi,N,M)
    tau=N/M;
    B=zeros(N,N);
    for n0=1:tau
        active=(1:M).'+(n0-1)*M;
        theta=Phi(n0,:);
        B(active,:)=G*diag(theta);
    end
end

