function [Pnn, Ncd] = cluster_elimination(Pn,Thr)
    N=length(Pn);
    Pmax=max(Pn);
    Pthr=Pmax*Thr;
    tmp=[];
    Ncd=N;
    for n0=1:N
        if Pn(n0)<Pthr
            tmp=[tmp,n0];
            Ncd=Ncd-1;
        end
    end
    Pnn=Pn;
    Pnn(tmp)=[];
end

