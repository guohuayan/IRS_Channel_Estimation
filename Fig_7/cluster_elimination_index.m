function [idx] = cluster_elimination_index(Pn,Thr)
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
    idx=tmp;
%     Pnn=Pn;
%     Pnn(tmp)=[];
end