function f = opt_phi(theta,Ew)
    [~,~,K]=size(Ew);
    ft=zeros(1,K);
    for k0=1:K
        E=Ew(:,:,k0);
        ft(k0)=real(theta'*E*theta);
    end
    f=min(ft);
end

