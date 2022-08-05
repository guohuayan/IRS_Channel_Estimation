function f = BALS_obj(R,Hg,Hu,phi_t)
    [M,K,tL]=size(R);
    f=0;
    for l0=1:tL
        t_theta=phi_t(:,l0);
        tmp=R(:,:,l0)-Hg*diag(t_theta)*Hu;
        err=sum(abs(tmp(:)).^2);
        f=f+err;
    end
end

