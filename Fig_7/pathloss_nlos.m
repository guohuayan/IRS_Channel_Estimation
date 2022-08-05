function [path] = pathloss_nlos(d,fc)
%     path1=32.4+21*log10(d)+20*log10(fc);
%     path2=22.4+35.3*log10(d)+21.3*log10(fc);
    path1=32.4+17.3*log10(d)+20*log10(fc);
    path2=22.4+31.9*log10(d)+20*log10(fc);
    path=max(path1,path2);
end

