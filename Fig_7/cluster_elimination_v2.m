function [Pnn] = cluster_elimination_v2(Pn,Thr)
%% just assign the eliminated cluster by p=0
    N=length(Pn);
    Pmax=max(Pn);
    Pthr=Pmax*Thr;
    Ncd=N;
    for n0=1:N
        if Pn(n0)<Pthr
            Pn(n0)=0;
        end
    end
    Pnn=Pn;
end